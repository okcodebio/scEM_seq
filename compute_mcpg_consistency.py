#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute per-barcode mCpG coverage and consistency metrics.
Usage:
    python compute_mcpg_consistency.py <input_dir> <output_summary.tsv>

Example:
    python compute_mcpg_consistency.py ./clean_txt mcpg_summary.tsv
"""

import os
import glob
import re
import sys
from collections import Counter
from typing import List

import pandas as pd


def extract_barcode_from_filename(path: str) -> str:
    """
    Extract barcode from filename using a regex pattern.

    Default pattern:
        r"_srtn_(.*?).deduplicated"

    If the pattern does not match, the basename (without extension)
    is returned as a fallback.

    Parameters
    ----------
    path : str
        Full path to the file.

    Returns
    -------
    barcode : str
        Extracted barcode or fallback string.
    """
    basename = os.path.basename(path)
    m = re.search(r"_srtn_(.*?).deduplicated", basename)
    if m:
        return m.group(1)
    # Fallback: use the basename without extension
    return os.path.splitext(basename)[0]


def main():
    if len(sys.argv) != 2 and len(sys.argv) != 3:
        print(
            "Usage:\n"
            "  python compute_mcpg_consistency.py <input_dir> <output_summary.tsv>"
        )
        sys.exit(1)

    directory_path = sys.argv[1]
    out_file = sys.argv[2]

    if not os.path.isdir(directory_path):
        print(f"[ERROR] Directory not found: {directory_path}")
        sys.exit(1)

    pattern = os.path.join(directory_path, "*cpg.txt")
    mey_files: List[str] = glob.glob(pattern)

    if not mey_files:
        print(f"[WARN] No files matching '*cpg.txt' in: {directory_path}")
        sys.exit(0)

    data_list = []

    for file_path in sorted(mey_files):
        print(f"[INFO] Processing: {file_path}")

        # Each file is expected to have: chr, pos, status (no header)
        temp_df = pd.read_csv(
            file_path,
            sep="\t",
            header=None,
            names=["chr", "pos", "status"],
        )

        # Sort by chromosome and genomic position
        temp_df_sorted = temp_df.sort_values(by=["chr", "pos"], ascending=[True, True])

        # List of all positions
        pos_list = temp_df_sorted["pos"].tolist()

        if len(pos_list) == 0:
            # No positions detected, skip or record zeros
            barcode = extract_barcode_from_filename(file_path)
            print(f"[WARN] No positions found in file: {file_path}")
            data_list.append([barcode, 0, 0, 0])
            continue

        # Count how many times each position is observed
        counts = Counter(pos_list)

        # Positions covered by >= 2 reads
        count_ge_2 = sum(1 for value in counts.values() if value >= 2)

        # Number of unique positions
        detect_pos = len(counts)

        # Count positions with inconsistent status
        status_2 = 0
        unique_positions = set(pos_list)
        for pos in unique_positions:
            status_list = temp_df_sorted[temp_df_sorted["pos"] == pos]["status"].tolist()
            are_all_same = len(set(status_list)) == 1
            if not are_all_same:
                status_2 += 1

        barcode = extract_barcode_from_filename(file_path)

        # Optional print for quick QC
        rep_pro = count_ge_2 / detect_pos if detect_pos > 0 else 0.0
        ins_pro = status_2 / detect_pos if detect_pos > 0 else 0.0
        print(
            f"[INFO] {barcode} | Detect_pos={detect_pos} | "
            f"Rep_pos={count_ge_2} ({rep_pro:.4f}) | "
            f"Ins_pos={status_2} ({ins_pro:.4f})"
        )

        data_list.append([barcode, detect_pos, count_ge_2, status_2])

    # Build summary table
    mcpg_status = pd.DataFrame(
        data_list,
        columns=["Barcode", "Detect_pos", "Rep_pos", "Ins_pos"],
    )

    # Avoid division-by-zero explicitly
    mcpg_status["Rep_pro"] = mcpg_status.apply(
        lambda row: row["Rep_pos"] / row["Detect_pos"] if row["Detect_pos"] > 0 else 0.0,
        axis=1,
    )
    mcpg_status["Ins_pro"] = mcpg_status.apply(
        lambda row: row["Ins_pos"] / row["Detect_pos"] if row["Detect_pos"] > 0 else 0.0,
        axis=1,
    )

    mcpg_status.to_csv(out_file, sep="\t", header=True, index=False)
    print(f"[INFO] Summary written to: {out_file}")


if __name__ == "__main__":
    main()
