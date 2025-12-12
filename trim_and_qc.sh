#!/usr/bin/env bash
set -euo pipefail

# Simple trimming + QC pipeline for paired-end scEM-seq reads
# Step 1: fixed-length trimming with cutadapt
# Step 2: quality control and filtering with fastp
#
# Usage:
#   bash trim_and_qc.sh sample_R1.fastq.gz sample_R2.fastq.gz
#
# Output:
#   sample_R1.cut.fastq.gz
#   sample_R2.cut.fastq.gz
#   sample_R1.cut.clean.fastq.gz
#   sample_R2.cut.clean.fastq.gz
#   sample.fastp.html
#   sample.fastp.json

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <R1.fastq.gz> <R2.fastq.gz>" >&2
    exit 1
fi

R1="$1"
R2="$2"

# Basic checks
if [ ! -f "$R1" ]; then
    echo "ERROR: R1 file not found: $R1" >&2
    exit 1
fi
if [ ! -f "$R2" ]; then
    echo "ERROR: R2 file not found: $R2" >&2
    exit 1
fi

# Example: sample_R1.fastq.gz  -> base_R1 = sample_R1
base_R1=$(basename "$R1" .fastq.gz)
base_R2=$(basename "$R2" .fastq.gz)

# Common prefix: sample
prefix=${base_R1%_R1*}

echo "[INFO] R1   : $R1"
echo "[INFO] R2   : $R2"
echo "[INFO] base_R1 : $base_R1"
echo "[INFO] prefix  : $prefix"

# ------------ 1) cutadapt fixed-length trimming ------------

cut_R1="${base_R1}.cut.fastq.gz"
cut_R2="${base_R2}.cut.fastq.gz"
cut_log="${prefix}.cutadapt.log"

echo "[INFO] Running cutadapt..."
cutadapt \
    -u 85 \               # trim 85 bases from 5' of R1
    -U 20 \               # trim 20 bases from 5' of R2
    --cores 6 \
    -o "${cut_R1}" \
    -p "${cut_R2}" \
    "$R1" \
    "$R2" \
    > "${cut_log}" 2>&1

echo "[INFO] cutadapt finished. Output:"
echo "       ${cut_R1}"
echo "       ${cut_R2}"
echo "       log: ${cut_log}"

# ------------ 2) fastp quality control ------------

clean_R1="${base_R1}.cut.clean.fastq.gz"
clean_R2="${base_R2}.cut.clean.fastq.gz"
fastp_html="${prefix}.fastp.html"
fastp_json="${prefix}.fastp.json"
fastp_log="${prefix}.fastp.log"

echo "[INFO] Running fastp..."
fastp \
    -i "${cut_R1}" \
    -I "${cut_R2}" \
    -o "${clean_R1}" \
    -O "${clean_R2}" \
    -h "${fastp_html}" \
    -j "${fastp_json}" \
    --thread 8 \
    > "${fastp_log}" 2>&1

echo "[INFO] fastp finished. Output:"
echo "       ${clean_R1}"
echo "       ${clean_R2}"
echo "       report: ${fastp_html}, ${fastp_json}"
echo "       log   : ${fastp_log}"

echo "[DONE] Trimming + QC pipeline completed."
