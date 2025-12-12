#!/usr/bin/env bash
# Convert per-cell SAM files to BAM, run Bismark deduplication,
# sort/index BAMs, and prepare methylation extractor commands.
#
# Usage:
#   1) Run this script in the directory containing per-cell SAM files
#      (generated e.g. by split_bam_by_cb.py).
#   2) Adjust the paths and parameters in the "User-configurable parameters"
#      section if needed.

set -euo pipefail

#############################################
# User-configurable parameters
#############################################

# Path to samtools (or simply "samtools" if it is in your $PATH)
SAMTOOLS_BIN="samtools"

# Path to Bismark deduplication script (deduplicate_bismark)
DEDUPLICATE_BISMARK_BIN="deduplicate_bismark"

# Path to Bismark methylation extractor (bismark_methylation_extractor)
BME_BIN="bismark_methylation_extractor"

# Path to Bismark genome folder (prepared with `bismark_genome_preparation`)
METHYLATION_REF="/path/to/bismark/genome/hg38"

# Number of threads for methylation extractor
THREADS=10

# Output file that stores all bismark_methylation_extractor commands
COMMAND_LIST_FILE="bismark_commands.txt"

#############################################
# 1. Convert SAM to BAM and deduplicate
#############################################

echo "[INFO] Step 1: Converting SAM to BAM and running deduplication..."

shopt -s nullglob
sam_files=( *.sam )

if [ "${#sam_files[@]}" -eq 0 ]; then
  echo "[WARN] No SAM files found in current directory."
else
  for sam in "${sam_files[@]}"; do
    base="${sam%.sam}"
    bam="${base}.bam"

    echo "[INFO] Converting ${sam} -> ${bam}"
    "${SAMTOOLS_BIN}" view -h -b "${sam}" > "${bam}"

    echo "[INFO] Running deduplicate_bismark on ${bam}"
    # -s: single-end, --barcode: indicate barcode information,
    # --bam: input is BAM
    "${DEDUPLICATE_BISMARK_BIN}" \
      -s \
      --barcode \
      --output_dir ./ \
      --bam \
      -o "${base}" \
      "${bam}"
    # deduplicate_bismark will produce ${base}.deduplicated.bam
  done
fi

#############################################
# 2. Sort and index deduplicated BAM files
#############################################

echo "[INFO] Step 2: Sorting and indexing deduplicated BAM files..."

dedup_bams=( *.deduplicated.bam )

if [ "${#dedup_bams[@]}" -eq 0 ]; then
  echo "[WARN] No *.deduplicated.bam files found."
else
  for id in "${dedup_bams[@]}"; do
    sorted_bam="${id%.bam}.srt.bam"
    echo "[INFO] Sorting ${id} -> ${sorted_bam}"
    "${SAMTOOLS_BIN}" sort -o "${sorted_bam}" "${id}"

    echo "[INFO] Indexing ${sorted_bam}"
    "${SAMTOOLS_BIN}" index "${sorted_bam}"
  done
fi

#############################################
# 3. Prepare Bismark methylation extractor commands
#############################################

echo "[INFO] Step 3: Preparing Bismark methylation extractor commands..."
echo "[INFO] Genome folder: ${METHYLATION_REF}"
echo "[INFO] Threads: ${THREADS}"

# Clear old command list if it exists
: > "${COMMAND_LIST_FILE}"

# You can choose whether to run methylation extractor on deduplicated BAMs
# or on sorted BAMs. Here we use *.deduplicated.bam by default.
for bam in *.deduplicated.bam; do
  [ -e "${bam}" ] || continue

  cmd="${BME_BIN} -s --no_overlap --comprehensive --bedGraph --counts --report --gzip \
--genome_folder ${METHYLATION_REF} \
--parallel ${THREADS} \
-o ./ \
${bam}"

  echo "${cmd}" >> "${COMMAND_LIST_FILE}"
done

echo "[INFO] Commands for Bismark methylation extractor written to: ${COMMAND_LIST_FILE}"
echo "[INFO] You can run them sequentially, or in parallel, for example:"
echo "       bash ${COMMAND_LIST_FILE}"
echo "       # or: cat ${COMMAND_LIST_FILE} | xargs -L 1 -P 4 bash -c"