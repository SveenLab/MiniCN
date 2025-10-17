# MiniCN
<img src="https://img.shields.io/badge/R%20tested-4.3.x-blue" alt="R version tested"> ![version](https://img.shields.io/badge/version-v0.1.0-blue) ![license](https://img.shields.io/badge/license-MIT-green)

## Copy Number Alteration Caller for Small Targeted Panels

MiniCN is a lightweight pipeline for DNA copy number calling on small, amplicon-based targeted sequencing panels (\< 1 Mb) using GATK-like coverage files. It supports either matched tumor–normal pairs or a pool-of-normals (PON) for coverage normalization.

MiniCN was originally created for a gene panel (\~60 kb) with amplicons covering the full coding regions of \~20 genes, and introns for only a subset of the genes.

## Installation

The easiest way to install the R package miniCN is via the devtools package:

```         
install.packages("devtools")
library(devtools)
devtools::install_github("SveenLab/MiniCN")
library(miniCN)
```
This package was constructed under **R version 4.3.2** (2023-10-31). 

---

## Quick Start

```         
library(miniCN)

# Example inputs bundled with the package
sample_csv <- system.file("extdata", "sample_sheet.csv", package = "MiniCN")
panel_file <- system.file("extdata", "panel_genes.txt", package = "MiniCN")

# Output directory
out <- "outputs"
if (!dir.exists(out)) {
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
}

# Run copy-number calling
if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
  run_caller(sample_csv = sample_csv, outdir = out, panel_path = panel_file)
}
```

---

## Output files (in `outdir`):

-   0_Excluded_Targets.tsv
-   1_Amplicons_CN.tsv
-   2_Samples_QC.tsv
-   3_Gene_CN_Final.tsv
-   \<sample\>\_CN_plots.pdf
  
---

## Input formats:

### Coverage files

These input files are derived from results of the **GATK `DepthOfCoverage`** tool. 
All coverage files names should have suffix `_coverage.sample_interval_summary` . 
They should be **tab-delimited text files** with a header and contain the following columns (case-sensitive):

Column | Description
|:---|:---
`target` (or `interval_id`, `interval`) | format `chr:start-end`
`total_coverage` (or `totalcount`, `totalcounts`, `count`, `counts`, `coverage`, `totalcoverage`) | total read coverage

> *Example GATK’s `DepthOfCoverage` command:*
>
> ```bash
> gatk DepthOfCoverage \
>   -R reference.fasta \
>   -O sample_coverage \
>   -L targets.interval_list \
>   -I sample.bam
> ```
>
> The output file `sample_coverage.sample_interval_summary` can be directly used as MiniCN input.

### Sample sheet (`.csv`)

Column | Description
|:---|:---
`Tumor Folder` | folder path with tumor coverage files
`Tumor File` | tumor coverage file name
`Normal Folder` | folder path with Normal coverage files
`Normal File` | Normal coverage file name

> Example:
>
> ```         
> Tumor Folder,Tumor File Name,Normal Folder,Normal File Name
> examples/,tumorA,examples/,normalA
> ```

**Using a PON**: set `Normal File Name` to `PON_median` and point `Normal Folder` to a directory containing more than 2 files.

> Example:
>
> ```         
> Tumor Folder,Tumor File Name,Normal Folder,Normal File Name
> examples/,tumorA,examples/PON_normals/,PON_median
> ```

### Panel genes

Plain text, one HGNC symbol per line (e.g., `EGFR`).

Symbols are validated via `org.Hs.eg.db` .

---

## Key Parameters

Flag | Description | Default
|:---|:---|:---
`pc` | pseudocount | 0.5
`HIGH_AMP_T` | high-amplification copy number threshold | log2(17/2) or ~3.09
`LOW_AMP_T` | low-amplification copy number threshold | log2(7/2) or ~1.81
`GAIN_T` | gain copy number threshold | log2(3/2) or ~0.58
`LOSS_T` | loss copy number threshold | log2(1/2) or -1.0
`DEEP_DEL_T` | deletion copy number threshold | -2.0
`MIN_AMPLICONS` | minimum number of amplicons supporting respective changes | 2.0
`MIN_GENE_Z` | minimum z-score cutoff for gains and losses | 3.0
`MAX_Q_VALUE` | maximum q-value cutoff for gains and losses | 0.05

## Citation

To cite package ‘MiniCN’ in publications use:

Nunes, L. (2025). MiniCN: Copy Number Alteration Caller for Small Targeted Panels. R package version 0.1.0. SveenLab, Dept. Molecular Oncology, ICR, Oslo University Hospital. URL: https://github.com/SveenLab/miniCN.
