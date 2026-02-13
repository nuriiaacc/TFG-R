# Instalar paquetes
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "DESeq2",
  "edgeR",
  "recount3",
  "ggplot2",
  "pheatmap"
))

library(TCGAbiolinks)
library(SummarizedExperiment)
library(data.table)

#Crear directorio para descargar datos GBM
setwd("C:/Users/nuria/Downloads")
getwd()
dir.create("GDCdata", recursive = TRUE, showWarnings = FALSE)

#OBTENER DATOS
#Obtener datos TCGA-GBM
query_tcga <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query_tcga, directory = "GDCdata")

gbm_data <- GDCprepare(query_tcga, directory = "GDCdata")
gbm_counts <- assay(gbm_data, "unstranded")
gbm_metadata <- colData(gbm_data)

cat("Genes:", nrow(gbm_counts), "\n")
cat("Muestras:", ncol(gbm_counts), "\n")

head(gbm_counts[, 1:5])
dim(gbm_counts)

#Obtener datos GTEx
# Cargar el archivo
gtex_expr <- fread(
  "C:/Users/nuria/OneDrive/Escritorio/CUARTO/TFG/Proyecto TFG-R/Proyecto TFG/gtex/gene_reads_v11_brain_cortex.gct",
  skip = 2,
  data.table = FALSE
)

# Extraer nombres de genes (ENSEMBL IDs)
gene_names <- gtex_expr$Name

# Convertir a matriz
gtex_counts <- as.matrix(gtex_expr[, -c(1, 2)])

# Añadir nombres de genes
rownames(gtex_counts) <- gene_names

head(gtex_counts[, 1:5])
dim(gtex_counts)

cat("Genes:", nrow(gtex_counts), "\n")
cat("Muestras:", ncol(gtex_counts), "\n")

#ANÁLISIS DIFERENCIAL

library(DESeq2)
library(ggplot2)
library(pheatmap)

#Genes comunes
cat("TCGA genes:", nrow(gbm_counts), "\n")
cat("GTEx genes:", nrow(gtex_counts), "\n")

common_genes <- intersect(rownames(gbm_counts), rownames(gtex_counts))
cat("Genes comunes:", length(common_genes), "\n")

#Combinar matrices
counts_combined <- cbind(
  gbm_counts[common_genes, ],
  gtex_counts[common_genes, ]
)

cat("Matriz combinada:", dim(counts_combined), "\n")

#Crear condición
condition <- factor(c(
  rep("GBM", ncol(gbm_counts)),
  rep("Normal", ncol(gtex_counts))
), levels = c("Normal", "GBM"))

print(table(condition))

#DESeq2-Expresión diferencial
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_combined),
  colData = data.frame(condition = condition),
  design = ~ condition
)

dds <- DESeq(dds)

#Resultados análisis de expresión diferencial
res <- results(dds, contrast = c("condition", "GBM", "Normal"))
res_df <- as.data.frame(res)
res_ordered <- res[order(res$padj), ]

cat("\nRESULTADOS\n")
summary(res_ordered)

#Genes expresados diferencialmente
deg <- res_ordered[res_ordered$padj < 0.01 & !is.na(res_ordered$padj), ]

cat("\nGenes significativos (FDR < 0.01):", nrow(deg), "\n")
cat("Upregulados en GBM:", sum(deg$log2FoldChange > 0), "\n")
cat("Downregulados en GBM:", sum(deg$log2FoldChange < 0), "\n")

cat("\nTop 20 genes:\n")
print(head(deg, 20))

#Guardar resultados
write.csv(as.data.frame(deg), 
          "genes_diferencialmente_expresados_GBM_vs_Normal.csv")



#VISUALIZACIÓN
library(ggplot2)
library(pheatmap)

#VOLCANO
res_df <- as.data.frame(res_ordered)
res_df$Gene <- rownames(res_df)

pdf("volcanoplot.pdf", width = 10, height = 8)
EnhancedVolcano(
  res_df,
  lab = res_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2,
  labSize = 3,
  max.overlaps = 20  # aumenta o disminuye el número de genes con etiqueta
)
dev.off()

#PCA
vsd <- vst(dds, blind = FALSE)

pdf("pca_plot.pdf", width = 10, height = 8)
plotPCA(vsd, intgroup = "condition") + 
  theme_minimal() +
  labs(title = "PCA: GBM vs Tejido Normal")
dev.off()

#MAPA DE CALOR
top_genes <- rownames(deg[1:50, ])
heatmap_data <- assay(vsd)[top_genes, ]

pdf("heatmap_top50.pdf", width = 12, height = 10)
pheatmap(heatmap_data,
         scale = "row",
         main = "Top 50 genes diferencialmente expresados",
         annotation_col = data.frame(condition = condition),
         show_colnames = FALSE,
         fontsize = 10)
dev.off()
