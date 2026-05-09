FROM python:3.11-slim

RUN apt-get update && apt-get install -y \
    git samtools tabix bgzip \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    numpy==1.26.4 \
    scipy==1.13.0 \
    pysam==0.22.0

RUN git clone --depth=1 --branch v1.1 \
    https://github.com/Illumina/SMNCopyNumberCaller.git /opt/smn && \
    cd /opt/smn && pip install --no-cache-dir -e .

ENV PATH="/opt/smn/smn_copy_number_caller:${PATH}"

LABEL maintainer="svcaller" \
      version="1.1" \
      description="SMNCopyNumberCaller v1.1 for SMN1/SMN2 copy number calling"
