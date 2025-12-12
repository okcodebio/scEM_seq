#!/usr/bin/env python3
"""
Add cell barcode (CB) and UMI information to paired-end FASTQ files.

For each read in R1, the first N bases are treated as the cell barcode and the
following M bases as the UMI. Reads whose barcode is not in the provided
whitelist are discarded. For reads that pass this filter, a new FASTQ header
containing CB and UMI is written for both R1 and the corresponding R2 read.

R1/R2 are assumed to be gzipped FASTQ files.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from typing import Dict, Tuple, Set

from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Add cell barcode and UMI to read names in paired-end FASTQ files."
    )
    parser.add_argument(
        "-i", "--in1", required=True, help="Input FASTQ.gz for read 1"
    )
    parser.add_argument(
        "-I", "--in2", required=True, help="Input FASTQ.gz for read 2"
    )
    parser.add_argument(
        "-o", "--out1", required=True, help="Output FASTQ.gz for read 1"
    )
    parser.add_argument(
        "-O", "--out2", required=True, help="Output FASTQ.gz for read 2"
    )
    parser.add_argument(
        "-cb",
        "--cb",
        required=True,
        help="Plain-text file with one valid cell barcode per line",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        required=True,
        help="Prefix for output cell-barcode summary file",
    )
    parser.add_argument(
        "--cb-len",
        type=int,
        default=17,
        help="Length of cell barcode in R1 (default: 17)",
    )
    parser.add_argument(
        "--umi-len",
        type=int,
        default=12,
        help="Length of UMI in R1 immediately following the barcode (default: 12)",
    )
    return parser.parse_args()


def open_gzip_text(path: str | Path, mode: str = "rt"):
    """
    Convenience wrapper for opening gzipped text files.
    """
    return gzip.open(path, mode)


def read_cell_barcodes(path: str | Path) -> Set[str]:
    """Read valid cell barcodes from a text file (one per line)."""
    cb_set: Set[str] = set()
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cb_set.add(line)
    if not cb_set:
        raise ValueError(f"No barcodes found in {path}")
    return cb_set


def process_read1(
    file_path: str | Path,
    out_file: str | Path,
    cb_list: Set[str],
    cb_len: int,
    umi_len: int,
) -> Tuple[Dict[str, str], Dict[str, int], int, int]:
    """
    Process R1: filter by CB whitelist, extract CB+UMI, rewrite read name.

    Returns
    -------
    cb_umi_dict : dict
        Mapping from original read ID -> new FASTQ header line (including '@').
    cb_count : dict
        Mapping from CB sequence -> number of reads passing filter.
    total_reads : int
        Total number of reads seen in R1.
    invalid_cb_reads : int
        Number of reads whose CB was not in the whitelist.
    """
    cb_umi_dict: Dict[str, str] = {}
    cb_count: Dict[str, int] = {}
    total_reads = 0
    invalid_cb_reads = 0

    with open_gzip_text(file_path, "rt") as f_in, open_gzip_text(out_file, "wt") as f_out:
        for record in SeqIO.parse(f_in, "fastq"):
            total_reads += 1

            seq = str(record.seq)
            if len(seq) < cb_len + umi_len:
                # Too short to contain CB + UMI; skip
                invalid_cb_reads += 1
                continue

            cb = seq[:cb_len]
            umi = seq[cb_len : cb_len + umi_len]

            if cb not in cb_list:
                invalid_cb_reads += 1
                continue

            # New read name:
            # Example: @<orig_id>_<CB>_<UMI> 1:N:0:<UMI>
            # (Space between name and Illumina flags, no tabs.)
            new_read_name = f"@{record.id}_{cb}_{umi} 1:N:0:{umi}"

            # Convert qualities back to ASCII 33-based string
            qual = "".join(chr(q + 33) for q in record.letter_annotations["phred_quality"])

            f_out.write(f"{new_read_name}\n{seq}\n+\n{qual}\n")

            cb_umi_dict[record.id] = new_read_name
            cb_count[cb] = cb_count.get(cb, 0) + 1

    return cb_umi_dict, cb_count, total_reads, invalid_cb_reads


def process_read2(
    file_path: str | Path,
    out_file: str | Path,
    cb_umi_dict: Dict[str, str],
) -> int:
    """
    Process R2: only keep reads whose read ID is present in R1 mapping.

    Returns
    -------
    missing_r2 : int
        Number of R2 reads that did not have a matching R1 entry.
    """
    missing_r2 = 0

    with open_gzip_text(file_path, "rt") as f_in, open_gzip_text(out_file, "wt") as f_out:
        for record in SeqIO.parse(f_in, "fastq"):
            read_id = record.id
            if read_id not in cb_umi_dict:
                missing_r2 += 1
                continue

            # Replace " 1:N:0:" with " 2:N:0:" in the stored header
            r1_header = cb_umi_dict[read_id]
            r2_header = r1_header.replace(" 1:N:0:", " 2:N:0:")

            seq = str(record.seq)
            qual = "".join(chr(q + 33) for q in record.letter_annotations["phred_quality"])

            f_out.write(f"{r2_header}\n{seq}\n+\n{qual}\n")

    return missing_r2


def write_cb_summary(
    cb_count: Dict[str, int],
    total_reads: int,
    invalid_cb_reads: int,
    missing_r2: int,
    prefix: str,
) -> None:
    """Write barcode usage summary and print basic stats to stderr."""
    kept_reads = sum(cb_count.values())

    # Print summary to stderr (so it doesn't mix with downstream pipes)
    sys.stderr.write(
        (
            f"[SUMMARY] Total R1 reads: {total_reads}\n"
            f"[SUMMARY] Reads with valid CB: {kept_reads} "
            f"({kept_reads / total_reads:.4f} of total)\n"
            f"[SUMMARY] Reads with invalid/absent CB: {invalid_cb_reads}\n"
            f"[SUMMARY] R2 reads without matching R1 (dropped): {missing_r2}\n"
        )
    )

    # Save per-barcode counts
    summary_path = f"{prefix}_cell_barcode_info.txt"
    sorted_items = sorted(cb_count.items(), key=lambda x: x[1], reverse=True)
    with open(summary_path, "w") as f:
        for cb, count in sorted_items:
            f.write(f"{cb}\t{count}\n")


def main() -> None:
    args = parse_args()

    cb_list = read_cell_barcodes(args.cb)

    # Process R1
    cb_umi_dict, cb_count, total_reads, invalid_cb_reads = process_read1(
        args.in1,
        args.out1,
        cb_list,
        cb_len=args.cb_len,
        umi_len=args.umi_len,
    )

    # Process R2
    missing_r2 = process_read2(
        args.in2,
        args.out2,
        cb_umi_dict,
    )

    # Summary
    write_cb_summary(
        cb_count,
        total_reads,
        invalid_cb_reads,
        missing_r2,
        prefix=args.prefix,
    )


if __name__ == "__main__":
    main()
