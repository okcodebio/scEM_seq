#!/usr/bin/env python3
"""
Split a BAM file into per-cell SAM files based on cell barcode (CB) encoded
in the read name.

Assumption:
    Read name format is:  <CB>_<UMI>_<READNUM>_<R1/R2>
Example:
    AAACCCAAGAAAGTTA_SOMEUMI_001_R1

This script extracts the CB field and writes reads into separate SAM files,
each containing its own SAM header.

Author: YOUR_NAME
Repository: https://github.com/<your_repo>
"""

import os
import sys
import pysam


def extract_cb_from_read_name(read_name: str):
    """
    Extract the cell barcode (CB) from the read name.

    Expected read name structure:
        CB_UMI_readNum_R1/R2

    Returns:
        CB (str) if extraction is successful, otherwise None.
    """
    try:
        parts = read_name.split("_")
        if len(parts) < 4:
            return None
        cb = parts[0]  # first segment is CB
        return cb
    except Exception:
        return None


def main():
    if len(sys.argv) != 2:
        print("Usage: python split_bam_by_cb.py <input.bam>")
        sys.exit(1)

    bam_file_path = sys.argv[1]

    if not os.path.exists(bam_file_path):
        print(f"[ERROR] File not found: {bam_file_path}")
        sys.exit(1)

    basename = os.path.basename(bam_file_path).replace(".bam", "")

    # Directory for output SAM files
    output_dir = "split_bam_files_sam"
    os.makedirs(output_dir, exist_ok=True)

    print(f"[INFO] Reading BAM: {bam_file_path}")
    bamfile = pysam.AlignmentFile(bam_file_path, "rb")

    # Convert header dict → SAM formatted header text
    header_text = str(pysam.AlignmentHeader.from_dict(bamfile.header))

    written_cbs = set()  # record CBs that already have a SAM header written

    for read in bamfile.fetch(until_eof=True):
        cb = extract_cb_from_read_name(read.query_name)

        if cb is None:
            # If read name does not match expected format
            print(f"[WARN] Failed to parse CB from read: {read.query_name}")
            continue

        out_path = os.path.join(output_dir, f"{basename}_{cb}.sam")

        # If CB first appears → write header + first read
        if cb not in written_cbs:
            with open(out_path, "w") as f:
                f.write(header_text)
                f.write(read.to_string() + "\n")
            written_cbs.add(cb)
        else:
            # Append subsequent reads
            with open(out_path, "a") as f:
                f.write(read.to_string() + "\n")

    bamfile.close()
    print(f"[INFO] Splitting finished. Output directory: {output_dir}")


if __name__ == "__main__":
    main()
