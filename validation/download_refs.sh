#!/usr/bin/env bash
# Download GRCh38 reference, GIAB truth sets, AnnotSV db, EH catalog, cytobands.
# Run from the SVcaller project root.
set -euo pipefail

REF_DIR="/data/alvin/ref/GRCh38"
GIAB_DIR="/data/alvin/ref/GIAB"
mkdir -p "${REF_DIR}" "${GIAB_DIR}"

echo "=== Downloading GRCh38 reference ==="
if [ ! -f "${REF_DIR}/GRCh38.fasta" ]; then
    wget -q -O "${REF_DIR}/GRCh38.fasta.gz" \
        "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    bgzip -d "${REF_DIR}/GRCh38.fasta.gz"
    samtools faidx "${REF_DIR}/GRCh38.fasta"
    samtools dict  "${REF_DIR}/GRCh38.fasta" > "${REF_DIR}/GRCh38.dict"
    echo "Reference downloaded and indexed."
else
    echo "Reference already exists, skipping."
fi

echo "=== Building BWA-MEM2 index (via Docker) ==="
if [ ! -f "${REF_DIR}/GRCh38.fasta.0123" ]; then
    docker run --rm \
      -v "${REF_DIR}:/ref" \
      quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_1 \
      bwa-mem2 index /ref/GRCh38.fasta
    echo "BWA-MEM2 index built."
else
    echo "BWA-MEM2 index already exists, skipping."
fi

echo "=== Downloading GRCh38 cytobands ==="
if [ ! -f "/data/alvin/SVcaller/assets/GRCh38_cytobands.txt" ] || \
   [ "$(wc -l < /data/alvin/SVcaller/assets/GRCh38_cytobands.txt)" -lt 100 ]; then
    wget -q -O /data/alvin/tmp/cytobands.txt.gz \
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"
    bgzip -d -c /data/alvin/tmp/cytobands.txt.gz \
        > /data/alvin/SVcaller/assets/GRCh38_cytobands.txt
    echo "Cytobands downloaded ($(wc -l < /data/alvin/SVcaller/assets/GRCh38_cytobands.txt) bands)."
else
    echo "Cytobands already present."
fi

echo "=== Downloading GIAB SV truth sets ==="

# --- v0.6 (legacy; deletion-biased; ~12K SVs) ---
GIAB_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6"
for FNAME in "HG002_SVs_Tier1_v0.6.vcf.gz" "HG002_SVs_Tier1_v0.6.vcf.gz.tbi" "HG002_SVs_Tier1_v0.6.bed"; do
    DEST="${GIAB_DIR}/$(echo "${FNAME}" | sed 's/HG002_SVs_Tier1_v0.6/HG002_SV_v0.6/')"
    if [ ! -f "${DEST}" ]; then
        wget -q -O "${DEST}" "${GIAB_BASE}/${FNAME}"
        echo "Downloaded: ${DEST}"
    fi
done

# --- v1.0 (preferred; multi-platform HiFi+ONT+short-read; ~75K SVs incl. insertions) ---
# Check the latest release path at:
#   https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/
# Update GIAB_V1_BASE and filenames below once the path is confirmed.
GIAB_V1_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v1.0"
V1_VCF="HG002_SVs_Tier1_v1.0.vcf.gz"
if [ ! -f "${GIAB_DIR}/HG002_SV_v1.0.vcf.gz" ]; then
    if wget -q --spider "${GIAB_V1_BASE}/${V1_VCF}" 2>/dev/null; then
        wget -q -O "${GIAB_DIR}/HG002_SV_v1.0.vcf.gz"     "${GIAB_V1_BASE}/${V1_VCF}"
        wget -q -O "${GIAB_DIR}/HG002_SV_v1.0.vcf.gz.tbi" "${GIAB_V1_BASE}/${V1_VCF}.tbi"
        echo "Downloaded GIAB SV v1.0: ${GIAB_DIR}/HG002_SV_v1.0.vcf.gz"
    else
        echo "WARN: GIAB SV v1.0 URL not yet reachable — verify path at NCBI FTP and rerun."
        echo "      Falling back to v0.6 for benchmarking."
        echo "      Pass --giab_truth ${GIAB_DIR}/HG002_SV_v0.6.vcf.gz until v1.0 is available."
    fi
else
    echo "GIAB SV v1.0 already present."
fi

echo "=== Downloading ExpansionHunter catalog ==="
EH_URL="https://github.com/Illumina/RepeatCatalogs/raw/main/hg38/variant_catalog.json"
if [ ! -s "/data/alvin/SVcaller/assets/eh_catalog.json" ] || \
   grep -q '"note": "Stub catalog"' /data/alvin/SVcaller/assets/eh_catalog.json 2>/dev/null; then
    wget -q -O /data/alvin/SVcaller/assets/eh_catalog.json "${EH_URL}"
    echo "EH catalog downloaded."
else
    echo "EH catalog already present."
fi

echo ""
echo "=== All reference data ready ==="
echo "  Reference: ${REF_DIR}/GRCh38.fasta"
echo "  GIAB truth: ${GIAB_DIR}/HG002_SV_v0.6.vcf.gz"
echo "  EH catalog: /data/alvin/SVcaller/assets/eh_catalog.json"
