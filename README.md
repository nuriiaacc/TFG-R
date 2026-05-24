# TFG-R
# Transcriptomic and genomic analysis of transcription factors (TFs) in glioblastoma multiforme (GBM)

This repository contains a bioinformatics pipeline developed in R to download, process, and analyze RNA-seq and somatic mutation data from GBM patients available at The Cancer Genome Atlas (TCGA).

The workflow focuses on identifying differentially expressed TFs, performing unsupervised molecular subtyping, evaluating clinical survival rates, and predicting potential drug targets.

## Overview

The analysis is structured into different phases:

1. Data acquisition (TCGA-GBM): Automated retrieval of raw RNA-seq counts (STAR-Counts) for primary tumor and solid tissue normal samples using `TCGAbiolinks`.
2. Transcription factor filtering: Querying Ensembl using `biomaRt` to isolate genes mapped to the transcriptional regulation GO term (`GO:0003700`).
3. Differential expression analysis: Statistical testing with `DESeq2` (tumor vs. normal).
4. Functional enrichment analysis: Gene Ontology (GO-Biological Processes) over-representation analysis using `clusterProfiler` to identify altered pathways within each subgroup.
5. Survival analysis: Kaplan-Meier survival curve estimation using `survival` and `survminer` to assess the clinical impact of the identified molecular subgroups.
6. Cancer genomics and drug repurposing: Somatic mutation interaction profiling (excluding $IDH1$-mutated samples) via `maftools` and drug candidate discovery targeting upregulated TFs through DSigDB via `enrichR`.

## Prerequisites and installation

To run this pipeline, you need **R (version $\ge$ 4.0)** along with the following CRAN and Bioconductor packages:

### Required Packages
```R
# Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "DESeq2", "biomaRt", "clusterProfiler", 
                       "org.Hs.eg.db", "enrichplot", "survminer", "maftools"))

# CRAN dependencies
install.packages(c("dplyr", "ggplot2", "ggrepel", "pheatmap", "survival", 
                   "data.table", "enrichR"))
```
## Usage

1. Clone this repository to your local machine:

```bash
git clone [https://github.com/nuriiaacc/TFG-R.git](https://github.com/nuriiaacc/TFG-R.git)
```

2. Open the project or the directory in RStudio.

> Before running the pipeline, copy and paste this code into your RStudio console to create the required output directory:

if (!dir.exists("output")) dir.create("output")

3. Run the main R script (`TFG-script.R`). 

> The initial execution might take several minutes to download the large TCGA datasets and query the Ensembl BioMart server.

## Repository structure

* `TFG-script.R`: The core R script containing the workflow.
* `GDCdata/`: (Excluded via .gitignore) Local folder where raw TCGA GDC data files are stored.
* `output/`: Directory where all analytic results are saved:
    * `Resultados_DEA_TFs_GBM.csv`: Main spreadsheet containing the differential expression statistics.
    * `Volcano_TFs.png` & `Heatmap_TFs.pdf`: Visualization plots of the TF expression landscape.
    * `PCA_plot.pdf`: Principal Component Analysis showing tumor/normal separation and K-means subtyping.
    * `Supervivencia_GBM_Subgrupos.pdf`: Kaplan-Meier survival curves and risk tables.
    * `GO_General_Tumor_vs_Normal.pdf`: Dotplots illustrating over-represented Biological Processes.
    * `Dianas_Terapeuticas_Enrichr.pdf`: Bar chart highlighting top potential pharmaceutical drug matches.

---

## Author

* **Nuria Cuenca Campoy** - [*GitHub Profile*](https://github.com/nuriiaacc) - [*LinkedIn*](www.linkedin.com/in/nuria-cuenca-campoy-9396a5227)
