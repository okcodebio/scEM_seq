# scEM_seq
## Description

This code provides a comprehensive computational workflow for processing and analyzing scEM-seq data, from raw sequencing reads to downstream cell-type annotation and clustering. Paired-end FASTQ reads are demultiplexed by cell barcode, trimmed, aligned with Bismark, and processed to generate single-cell DNA methylation profiles and a methylation matrix. Downstream analyses perform VMR detection, matrix imputation, Seurat-based clustering, and cell-type annotation using reference methylomes or scMethID.

This code supports both reference-based and reference-free analyses, enabling cluster-level DNA methylation profiling in immune cells as well as complex tissues such as prostate cancer.

<p align="center">
  <img src="https://github.com/okcodebio/scEM_seq/blob/main/scEM_seq_framework_upload.png"
       width="50%"
       height="auto"
       alt="scEM_seq Framework">
</p>

## Key Features

- **An integrated, end-to-end scEM-seq analysis framework** from raw reads to cell-type–resolved methylation profiles  
- **Single-cell–aware demultiplexing** with explicit barcode and UMI modeling  
- **VMR-centric methylation representation**, enabling robust clustering and characterization of epigenetic heterogeneity  
- **Flexible cell-type inference** with or without reference atlases, supporting both immune and solid tumor contexts

## Installation

```bash
git clone https://github.com/okcodebio/scEM_seq.git
```

## Dependencies

- Cutadapt 4.5
- fastp 0.23.3
- Bowtie2 2.5.1
- Bismark 0.24.1
- Python >= 3.8
- NumPy >= 1.21.0
- pandas >= 1.3.0
