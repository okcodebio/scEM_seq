#!/usr/bin/env bash
# Simple Bismark + BAM merge + sorting pipeline
# Author: YOUR_NAME
# Description:
#   1) Run Bismark mapping on all FASTQ files in the current directory
#   2) Merge BAMs from R1 and R2 folders
#   3) Sort merged BAM by read name and by genomic coordinate

set -euo pipefail

#############################################
# User-configurable parameters
#############################################

# Path to Bismark executable (or just "bismark" if it is in your $PATH)
BISMARK_BIN="bismark"

# Path to Bowtie2 executables
BOWTIE2_PATH="/path/to/bowtie2/bin"

# Path to Bismark genome folder (prepared with `bismark_genome_preparation`)
GENOME_DIR="/path/to/bismark/genome/hg38"

# Threads and memory for samtools sorting
THREADS=10
SORT_MEMORY="10G"

# Sample ID used for naming merged BAM files
SAMPLE_ID="${1:-sample}"   # can be given as first argument; defaults to "sample"

# Directories containing BAM files from R1 and R2 (relative to current workdir)
R1_BAM_DIR="../R1_fq"
R2_BAM_DIR="../R2_fq"

#############################################
# 1. Run Bismark mapping
#############################################

WORKDIR="$(pwd)"

# Collect FASTQ files in the current directory
# (you can adjust the patterns to *.fq, *.fastq.gz, etc. as needed)
FASTQ_FILES=( *.fastq *.fastq.gz )
if [ "${#FASTQ_FILES[@]}" -eq 0 ]; then
  echo "[ERROR] No FASTQ files (*.fastq or *.fastq.gz) found in ${WORKDIR}" >&2
  exit 1
fi

echo "[INFO] Found ${#FASTQ_FILES[@]} FASTQ file(s): ${FASTQ_FILES[*]}"
echo "[INFO] Running Bismark..."

"${BISMARK_BIN}" -q \
  --path_to_bowtie "${BOWTIE2_PATH}" \
  -N 1 \
  --parallel 4 \
  --non_directional \
  --bam \
  --temp_dir "${WORKDIR}" \
  -un "${WORKDIR}" \
  --bowtie2 \
  --genome "${GENOME_DIR}" \
  -o "${WORKDIR}" \
  "${FASTQ_FILES[@]}"

echo "[INFO] Bismark mapping finished."

#############################################
# 2. Merge BAM files from R1 and R2
#############################################

echo "[INFO] Collecting BAM files for merging..."

# Adjust the patterns according to your real filenames
R1_BAMS=( "${R1_BAM_DIR}"/*/*cut.clean_bismark_bt2.bam )
R2_BAMS=( "${R2_BAM_DIR}"/*/*cut.clean_bismark_bt2.bam )

if [ "${#R1_BAMS[@]}" -eq 0 ] || [ "${#R2_BAMS[@]}" -eq 0 ]; then
  echo "[ERROR] No BAM files found in ${R1_BAM_DIR} or ${R2_BAM_DIR}" >&2
  exit 1
fi

MERGED_BAM="${SAMPLE_ID}_R1_R2.bam"

echo "[INFO] Merging BAM files into ${MERGED_BAM}..."

bamtools merge \
  $(printf -- "-in %s " "${R1_BAMS[@]}") \
  $(printf -- "-in %s " "${R2_BAMS[@]}") \
  -out "${MERGED_BAM}" \
  > bamtools_merge.out 2>&1

echo "[INFO] BAM merge finished. Output: ${MERGED_BAM}"

#############################################
# 3. BAM sorting
#############################################

# 3.1 Sort by read name (facilitates downstream processing by barcode/read-pair)
echo "[INFO] Sorting BAM by read name..."

NAME_SORTED_BAM="${SAMPLE_ID}_name_sorted.bam"

samtools sort -n \
  -m "${SORT_MEMORY}" \
  --threads "${THREADS}" \
  -o "${NAME_SORTED_BAM}" \
  "${MERGED_BAM}" \
  > sort_by_name.out 2>&1

echo "[INFO] Name-sorted BAM: ${NAME_SORTED_BAM}"

# 3.2 Sort by genomic coordinate
#     (convenient for detecting overlapping reads, duplicate removal, etc.)
echo "[INFO] Sorting BAM by genomic coordinate..."

COORD_SORTED_BAM="${SAMPLE_ID}_coord_sorted.bam"

samtools sort \
  --threads "${THREADS}" \
  -o "${COORD_SORTED_BAM}" \
  "${MERGED_BAM}" \
  > sort_by_coord.out 2>&1

echo "[INFO] Coordinate-sorted BAM: ${COORD_SORTED_BAM}"

echo "[INFO] Pipeline finished successfully."
