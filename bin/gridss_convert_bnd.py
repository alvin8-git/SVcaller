#!/usr/bin/env python3
"""Convert GRIDSS BND VCF to simple SV VCF (DEL/DUP/INV).

GRIDSS outputs all variants as BND (breakend) pairs. This script pairs
records by their EVENT tag, classifies same-chromosome pairs as DEL/DUP/INV
based on ALT bracket orientation, applies QUAL thresholds, and emits a
standard SV VCF that Jasmine can merge correctly.

ALT bracket orientation rules (VCF 4.2 spec):
  DEL: seq[chr:pos2[ at p1  AND  ]chr:pos1]seq at p2  (seq before/after bracket)
  DUP: ]chr:pos2]seq at p1  AND  seq[chr:pos1[ at p2  (seq after/before bracket)
  INV-h: seq[chr:pos2[ at p1 AND seq[chr:pos1[ at p2  (both seq-before)
  INV-t: ]chr:pos2]seq at p1 AND ]chr:pos1]seq at p2  (both seq-after)

Only same-chromosome pairs with SVLEN >= min_svlen pass; inter-chr pairs
and sub-threshold calls are discarded.
"""
import argparse, gzip, re, sys
from collections import defaultdict


def open_vcf(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def bracket_orient(alt: str):
    """Return (leading, bracket_char) for the sequence portion of a BND ALT.

    leading=True  → bases precede the first bracket  e.g. `seq[chr:pos[`
    leading=False → bases follow the last bracket     e.g. `]chr:pos]seq`
    """
    # Find leftmost bracket character
    for i, c in enumerate(alt):
        if c in '[]':
            return i > 0, c
    return False, None


def classify_pair(alt1: str, alt2: str):
    """Classify a BND pair into DEL/DUP/INV or None."""
    lead1, br1 = bracket_orient(alt1)
    lead2, br2 = bracket_orient(alt2)
    if br1 is None or br2 is None:
        return None
    # Trailing2: does alt2 have bases AFTER its last bracket?
    last_br2 = max(alt2.rfind('['), alt2.rfind(']'))
    trailing2 = last_br2 < len(alt2) - 1

    if lead1 and trailing2:
        return 'DEL'
    if not lead1 and not trailing2:
        return 'DUP'
    # INV: both leading or both trailing
    if lead1 and not trailing2:
        return 'INV'
    if not lead1 and trailing2:
        return 'INV'
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('vcf', help='GRIDSS BND VCF (gz or plain)')
    ap.add_argument('--out', required=True, help='Output VCF path')
    ap.add_argument('--min-qual-del', type=float, default=500,
                    help='Min QUAL for DEL calls (default 500)')
    ap.add_argument('--min-qual-dup', type=float, default=1000,
                    help='Min QUAL for DUP calls (default 1000)')
    ap.add_argument('--min-qual-inv', type=float, default=1000,
                    help='Min QUAL for INV calls (default 1000)')
    ap.add_argument('--min-svlen', type=int, default=50,
                    help='Min SV length bp (default 50)')
    ap.add_argument('--max-svlen', type=int, default=100_000_000,
                    help='Max SV length bp (default 100 Mb)')
    args = ap.parse_args()

    qual_thresh = {'DEL': args.min_qual_del, 'DUP': args.min_qual_dup,
                   'INV': args.min_qual_inv}

    # First pass: collect PASS records keyed by EVENT
    events = defaultdict(list)
    header_lines = []

    with open_vcf(args.vcf) as fh:
        for line in fh:
            if line.startswith('#'):
                header_lines.append(line.rstrip('\n'))
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8:
                continue
            chrom, pos_s, vid, ref, alt, qual_s, filt, info = fields[:8]
            if filt != 'PASS':
                continue
            info_d = {}
            for tok in info.split(';'):
                if '=' in tok:
                    k, v = tok.split('=', 1)
                    info_d[k] = v
                else:
                    info_d[tok] = True
            event = info_d.get('EVENT', '')
            if not event:
                continue
            try:
                qual = float(qual_s)
                pos  = int(pos_s)
            except ValueError:
                continue
            events[event].append({
                'chrom': chrom, 'pos': pos, 'alt': alt,
                'qual': qual, 'ref': ref, 'vid': vid,
            })

    # Second pass: convert paired BND events → simple SV records
    out_records = []
    seen_events = set()

    for event, recs in events.items():
        if len(recs) != 2 or event in seen_events:
            continue
        seen_events.add(event)

        a, b = sorted(recs, key=lambda r: (r['chrom'], r['pos']))
        if a['chrom'] != b['chrom']:
            continue

        svlen = b['pos'] - a['pos']
        if svlen < args.min_svlen or svlen > args.max_svlen:
            continue

        svtype = classify_pair(a['alt'], b['alt'])
        if svtype is None:
            continue

        qual = max(a['qual'], b['qual'])
        if qual < qual_thresh[svtype]:
            continue

        chrom = a['chrom']
        start = a['pos']
        end   = b['pos']

        info_parts = [
            f'SVTYPE={svtype}',
            f'END={end}',
            f'SVLEN={-svlen if svtype == "DEL" else svlen}',
            'IMPRECISE',
        ]

        out_records.append((
            chrom, start, f'GRIDSS_{event}_{svtype}',
            a['ref'], f'<{svtype}>',
            f'{qual:.1f}', 'PASS',
            ';'.join(info_parts),
            'GT', '0/1',
        ))

    # Write output VCF
    with open(args.out, 'w') as out:
        # Emit minimal header from original + add needed INFO/FORMAT lines
        emitted = False
        for hline in header_lines:
            if hline.startswith('#CHROM') and not emitted:
                out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
                out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n')
                out.write('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT">\n')
                out.write('##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">\n')
                out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                out.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
                emitted = True
            # Skip duplicate INFO/FORMAT/FILTER lines that GRIDSS already emits
            if re.match(r'^##(INFO|FORMAT|FILTER)=', hline):
                continue
            out.write(hline + '\n')

        out_records.sort(key=lambda r: (r[0], r[1]))
        for rec in out_records:
            out.write('\t'.join(str(x) for x in rec) + '\n')

    print(f'gridss_convert_bnd: wrote {len(out_records)} SVs to {args.out}', file=sys.stderr)
    # Summary by type
    from collections import Counter
    type_counts = Counter(r[4].strip('<>') for r in out_records)
    for sv, n in sorted(type_counts.items()):
        print(f'  {sv}: {n}', file=sys.stderr)


if __name__ == '__main__':
    main()
