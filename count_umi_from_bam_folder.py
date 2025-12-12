#!/usr/bin/env python3
"""
Count UMI-like tags from BAM files in a folder.

This script:
    1) Scans a folder for BAM files matching "*srt.bam"
    2) For each BAM file:
       - Parses read names
       - Extracts the last underscore-separated field as a "UMI / tag"
       - Counts how many reads per tag
       - Saves a per-BAM detail table: "<bam_basename>.reads.info.txt"
    3) Summarizes, per BAM file (treated as one barcode / cell), the
       total number of unique tags and total reads.

Assumptions:
    - Read name structure contains underscore-delimited parts, and the
      last part is a tag you want to count (e.g. "R1", "R2" or UMI-like).
    - Barcode (cell ID) can be extracted from the BAM file name.
      By default, this script uses a regex pattern:

          r'_cc3_srtn_(.*?).deduplicated'

      i.e., it expects filenames that look like:
          <prefix>_cc3_srtn_<BARCODE>.deduplicated.srt.bam

      If your naming scheme is different, please adjust the regex
      in `extract_barcode_from_filename()`.

Usage:
    python count_umi_from_bam_folder.py <bam_folder> <output_summary.tsv>

Example:
    python count_umi_from_bam_folder.py ./bam_dir umi_summary.tsv
"""

import os
import re
import sys
import glob
from typing import Tuple, Dict

import pysam
import pandas as pd


def umi_count_bam(inbam: str) -> Tuple[int, int, pd.DataFrame]:
    """
    Count how many reads are associated with each tag (the last
    underscore-separated field in the read name).

    Parameters
    ----------
    inbam : str
        Path to an input BAM file.

    Returns
    -------
    umi_count : int
        Number of unique tags.
    reads_number : int
        Total number of reads (sum of per-tag counts).
    df : pandas.DataFrame
        DataFrame with columns ['tag', 'count'] listing counts per tag.
    """
    tag_count: Dict[str, int] = {}

    with pysam.AlignmentFile(inbam, "rb") as bam_file:
        for read in bam_file.fetch(until_eof=True):
            # Example read name: CB_UMI_001_R1
            # Here we take the last part, e.g. "R1"
            name_parts = read.query_name.split("_")
            tag = name_parts[-1] if name_parts else ""

            if tag not in tag_count:
                tag_count[tag] = 1
            else:
                tag_count[tag] += 1

    df = pd.DataFrame(list(tag_count.items()), columns=["tag", "count"])
    umi_count = df.shape[0]
    reads_number = int(df["count"].sum())

    return umi_count, reads_number, df


def extract_barcode_from_filename(bam_path: str) -> str:
    """
    Extract barcode (or sample ID) from BAM filename.

    Default pattern:
        r'_cc3_srtn_(.*?).deduplicated'

    This expects something like:
        path/to/XYZ_cc3_srtn_BARCODE.deduplicated.srt.bam

    If the pattern does not match, fall back to the BAM basename
    (without extension).

    Parameters
    ----------
    bam_path : str
        Full path to the BAM file.

    Returns
    -------
    barcode : str
        Extracted barcode or a fallback name.
    """
    basename = os.path.basename(bam_path)
    m = re.search(r"_cc3_srtn_(.*?).deduplicated", basename)
    if m:
        return m.group(1)
    # Fallback: use basename without extension
    return os.path.splitext(basename)[0]


def main():
    if len(sys.argv) != 3:
        print(
            "Usage: python count_umi_from_bam_folder.py "
            "<bam_folder> <output_summary.tsv>"
        )
        sys.exit(1)

    folder_path = sys.argv[1]
    out_f = sys.argv[2]

    if not os.path.isdir(folder_path):
        print(f"[ERROR] Folder not found: {folder_path}")
        sys.exit(1)

    # Find BAM files ending with 'srt.bam'
    bam_files = glob.glob(os.path.join(folder_path, "*srt.bam"))

    if not bam_files:
        print(f"[WARN] No '*srt.bam' files found in: {folder_path}")
        sys.exit(0)

    barcode_list = []

    for fbam in sorted(bam_files):
        umi_count, reads_number, temp_df = umi_count_bam(fbam)

        barcode = extract_barcode_from_filename(fbam)

        print(os.path.basename(fbam), barcode, umi_count, reads_number)

        # Collect summary info
        barcode_list.append([barcode, umi_count, reads_number])

        # Save per-BAM detail table
        # (keeping original style "<bam>.reads.infor.txt" if you prefer,
        #  but here we use ".reads.info.txt" as a more standard name)
        detail_out = os.path.basename(fbam) + ".reads.info.txt"
        temp_df.to_csv(detail_out, sep="\t", header=True, index=False)

    # Build final summary table
    umi_df = pd.DataFrame(
        barcode_list,
        columns=["barcode", "umi_count", "reads_count"],
    )

    print(umi_df)
    umi_df.to_csv(out_f, sep="\t", header=True, index=False)
    print(f"[INFO] Summary written to: {out_f}")


if __name__ == "__main__":
    main()
