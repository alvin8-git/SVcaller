#!/usr/bin/env python3
"""Recompute cross-caller support for translocations (TRA/BND) in a Jasmine-merged VCF.

Jasmine does not co-cluster interchromosomal breakends across callers: Manta and
Delly each emit the same translocation as separate single-caller records
(SUPP_VEC has one bit set), so every TRA lands at SUPP=1 and is wiped out by the
downstream SUPP>=2 consensus gate. On COLO829 this gave 0/26 translocation recall
despite the callers finding 26/26.

This step re-clusters TRA records by breakend proximity (regardless of which
caller emitted them), then for each cluster emits ONE representative whose SUPP /
SUPP_VEC reflect the number of distinct callers that agreed. Non-TRA records pass
through untouched. Mate records (A->B and B->A) and cross-caller duplicates fold
into the same cluster.

ponytail: O(n^2) within chrom-pair buckets; fine for ~10k TRA. If TRA counts ever
explode, bucket by binned coordinate instead.
"""
import argparse
import gzip
import re
import sys

BR = re.compile(r'[\[\]]([^\[\]:]+):(\d+)[\[\]]')


def open_maybe_gz(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)


def is_tra(alt, info):
    if 'SVTYPE=TRA' in info or 'SVTYPE=BND' in info:
        return True
    return '[' in alt or ']' in alt


def parse_vec(info):
    m = re.search(r'SUPP_VEC=([0-9]+)', info)
    return m.group(1) if m else None


def vec_bits(vec):
    return {i for i, ch in enumerate(vec) if ch == '1'}


def qual(q):
    try:
        return float(q)
    except (ValueError, TypeError):
        return -1.0


def breakend(chrom, pos, alt):
    m = BR.search(alt)
    if not m:
        return None
    return (chrom, int(pos), m.group(1), int(m.group(2)))


def same_event(a, b, w):
    """a,b are (c1,p1,c2,p2). Match if same unordered chrom-pair and both
    breakends within window w (in either mate orientation)."""
    c1, p1, c2, p2 = a
    d1, q1, d2, q2 = b
    if {c1, c2} != {d1, d2}:
        return False
    if c1 == d1 and abs(p1 - q1) <= w and abs(p2 - q2) <= w:
        return True
    if c1 == d2 and abs(p1 - q2) <= w and abs(p2 - q1) <= w:
        return True
    return False


def cluster(tras, w):
    """Greedy union by proximity, bucketed on the unordered chrom-pair so we only
    compare records that can possibly match. Returns list of index-lists."""
    buckets = {}
    for i, t in enumerate(tras):
        key = tuple(sorted((t['be'][0], t['be'][2])))
        buckets.setdefault(key, []).append(i)
    clusters = []
    for idxs in buckets.values():
        used = [False] * len(idxs)
        for a in range(len(idxs)):
            if used[a]:
                continue
            group = [idxs[a]]
            used[a] = True
            # transitive single-link: keep absorbing matches into the group
            changed = True
            while changed:
                changed = False
                for b in range(len(idxs)):
                    if used[b]:
                        continue
                    if any(same_event(tras[g]['be'], tras[idxs[b]]['be'], w) for g in group):
                        group.append(idxs[b])
                        used[b] = True
                        changed = True
            clusters.append(group)
    return clusters


def set_info(info, key, value):
    parts = [p for p in info.split(';') if p and p != '.']
    out, seen = [], False
    for p in parts:
        if p.startswith(key + '='):
            out.append(f'{key}={value}')
            seen = True
        else:
            out.append(p)
    if not seen:
        out.append(f'{key}={value}')
    return ';'.join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('vcf')
    ap.add_argument('--out', required=True)
    ap.add_argument('--window', type=int, default=1000)
    args = ap.parse_args()

    header, others, tras = [], [], []
    with open_maybe_gz(args.vcf) as fh:
        for line in fh:
            if line.startswith('#'):
                header.append(line)
                continue
            f = line.rstrip('\n').split('\t')
            alt, info = f[4], f[7]
            if is_tra(alt, info):
                be = breakend(f[0], f[1], alt)
                vec = parse_vec(info)
                if be is None or vec is None:
                    others.append(line)          # can't cluster -> leave as-is
                    continue
                tras.append({'f': f, 'be': be, 'vec': vec, 'q': qual(f[5])})
            else:
                others.append(line)

    collapsed = []
    for group in cluster(tras, args.window):
        members = [tras[i] for i in group]
        bits = set()
        veclen = max(len(m['vec']) for m in members)
        for m in members:
            bits |= vec_bits(m['vec'])
        supp = len(bits)
        new_vec = ''.join('1' if i in bits else '0' for i in range(veclen))
        rep = max(members, key=lambda m: m['q'])      # highest-QUAL representative
        info = set_info(rep['f'][7], 'SUPP', str(supp))
        info = set_info(info, 'SUPP_VEC', new_vec)
        rep['f'][7] = info
        collapsed.append('\t'.join(rep['f']) + '\n')

    with open(args.out, 'w') as out:
        out.writelines(header)
        out.writelines(others)
        out.writelines(collapsed)

    sys.stderr.write(
        f'tra_consensus: {len(tras)} TRA records -> {len(collapsed)} clusters '
        f'(window={args.window})\n')


def _selfcheck():
    # Manta(bit0) and Delly(bit1) emit the same translocation ~200bp apart -> SUPP=2.
    # A Manta-only TRA stays SUPP=1. A mate record folds into its event.
    global BR
    rows = [
        # chrom pos id ref alt qual filter info
        ['1', '1000', 'a', 'T', 'T]5:2000]', '50', 'PASS', 'SVTYPE=TRA;SUPP=1;SUPP_VEC=10'],
        ['1', '1150', 'b', 'T', 'T]5:2100]', '60', 'PASS', 'SVTYPE=TRA;SUPP=1;SUPP_VEC=01'],
        ['5', '2050', 'b2', 'T', ']1:1120]T', '40', 'PASS', 'SVTYPE=TRA;SUPP=1;SUPP_VEC=01'],  # Delly mate
        ['9', '7000', 'c', 'T', 'T]12:8000]', '30', 'PASS', 'SVTYPE=TRA;SUPP=1;SUPP_VEC=10'],  # Manta-only
    ]
    tras = []
    for f in rows:
        tras.append({'f': f, 'be': breakend(f[0], f[1], f[4]),
                     'vec': parse_vec(f[7]), 'q': qual(f[5])})
    clusters = cluster(tras, 1000)
    supps = sorted(len(set().union(*[vec_bits(tras[i]['vec']) for i in g])) for g in clusters)
    assert supps == [1, 2], f'expected [1,2] got {supps}'
    assert set_info('SVTYPE=TRA;SUPP=1;SUPP_VEC=10', 'SUPP', '2') == 'SVTYPE=TRA;SUPP=2;SUPP_VEC=10'
    assert set_info('.', 'SUPP', '3') == 'SUPP=3'
    print('tra_consensus self-check OK')


if __name__ == '__main__':
    if len(sys.argv) == 2 and sys.argv[1] == '--selfcheck':
        _selfcheck()
    else:
        main()
