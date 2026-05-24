# TFG-R
# Transcriptomic and genomic analysis of transcription factors (TFs) in glioblastoma multiforme (GBM)

This repository contains a bioinformatics pipeline developed in R to download, process, and analyze RNA-seq and somatic mutation data from GBM patients available at The Cancer Genome Atlas (TCGA).

The workflow focuses on identifying differentially expressed TFs, performing unsupervised molecular subtyping, evaluating clinical survival rates, and predicting potential drug targets.

## Pipeline Overview

The analysis is structured into different phases:

1. Data acquisition (TCGA-GBM): Automated retrieval of raw RNA-seq counts (STAR-Counts) for primary tumor and solid tissue normal samples using `TCGAbiolinks`.
2. Transcription factor filtering: Querying Ensembl using `biomaRt` to isolate genes mapped to the transcriptional regulation GO term (`GO:0003700`).
3. Differential expression analysis: Statistical testing with `DESeq2` (tumor vs. normal) and subsequent patient subtyping using K-means clustering on Principal Component 1 (PC1).
4. Functional enrichment analysis: Gene Ontology (GO-Biological Processes) over-representation analysis using `clusterProfiler` to identify altered pathways within each subgroup.
5. Survival analysis: Kaplan-Meier survival curve estimation using `survival` and `survminer` to assess the clinical impact of the identified molecular subgroups.
6. Cancer genomics and drug repurposing: Somatic mutation interaction profiling (excluding $IDH1$-mutated samples) via `maftools` and drug candidate discovery targeting upregulated TFs through DSigDB via `enrichR`.

## Prerequisites & Installation

To run this pipeline, you need **R (version $\ge$ 4.0$)** along with the following CRAN and Bioconductor packages:

### Required Packages
```R
# Bioconductor Dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "DESeq2", "biomaRt", "clusterProfiler", 
                       "org.Hs.eg.db", "enrichplot", "survminer", "maftools"))

# CRAN Dependencies
install.packages(c("dplyr", "ggplot2", "ggrepel", "pheatmap", "survival", 
                   "data.table", "enrichR"))
