# Instalar paquetes
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "DESeq2",
  "ggplot2",
  "pheatmap",
  "data.table",
  "EnhancedVolcano",
  "biomaRt"
))

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(data.table)
library(EnhancedVolcano)
library(biomaRt)

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

head(gbm_counts[, 1:5])
cat("Genes:", nrow(gbm_counts), "\n")
cat("Muestras:", ncol(gbm_counts), "\n")


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
cat("Genes:", nrow(gtex_counts), "\n")
cat("Muestras:", ncol(gtex_counts), "\n")

#ANÁLISIS DIFERENCIAL

cat("TCGA genes:", nrow(gbm_counts), "\n")
cat("GTEx genes:", nrow(gtex_counts), "\n")

#Genes comunes
common_genes <- intersect(rownames(gbm_counts), rownames(gtex_counts))
cat("Genes comunes:", length(common_genes), "\n")

#Combinar matrices
counts_combined <- cbind(
  gbm_counts[common_genes, ],
  gtex_counts[common_genes, ]
)

cat("Matriz combinada:", dim(counts_combined), "\n")

#Crear condición para el análisis de expresión diferencial
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

#FILTRAR LOS FACTORES DE TRANSCRIPCIÓN
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Usamos GO:0003700 = DNA-binding transcription factor activity
tfs_info <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                  filters = 'go',
                  values = 'GO:0003700',
                  mart = ensembl)

tf_ensembl_ids <- unique(tfs_info$ensembl_gene_id)
cat("Se han recuperado", length(tf_ensembl_ids), "factores de transcripción.\n")

#Quitar la versión Ensembl ID y añadir el nombre del gen
res_df$ensembl_base <- gsub("\\..*", "", rownames(res_df))
res_df$Symbol <- tfs_info$hgnc_symbol[match(res_df$ensembl_base, tfs_info$ensembl_gene_id)]

#Usar solo factores de transcripción
res_tfs <- res_df[res_df$ensembl_base %in% tf_ensembl_ids, ]

#Ordenar por padj
res_tfs_ordered <- res_tfs[order(res_tfs$padj), ]

cat("\nRESULTADOS (Solo Factores de Transcripción)\n")
summary(res_tfs_ordered)

# Genes expresados diferencialmente (TFs)
deg_tfs <- res_tfs_ordered[res_tfs_ordered$padj < 0.01 & !is.na(res_tfs_ordered$padj), ]

cat("\nTFs significativos (FDR < 0.01):", nrow(deg_tfs), "\n")
cat("TFs Upregulados en GBM:", sum(deg_tfs$log2FoldChange > 0), "\n")
cat("TFs Downregulados en GBM:", sum(deg_tfs$log2FoldChange < 0), "\n")

cat("\nTop 20 TFs:\n")
print(head(deg_tfs, 20))

# Guardar resultados
write.csv(as.data.frame(deg_tfs), 
          "TFs_diferencialmente_expresados_GBM_vs_Normal.csv")


#VISUALIZACIÓN DE RESULTADOS

library(EnhancedVolcano)

#VOLCANO
#Usando nombres Symbol, IDs de Ensembl.
pdf("volcanoplot_TFs.pdf", width = 10, height = 8)
EnhancedVolcano(
  res_tfs_ordered,
  lab = ifelse(is.na(res_tfs_ordered$Symbol), rownames(res_tfs_ordered), res_tfs_ordered$Symbol),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2,
  labSize = 4,
  title = 'Glioblastoma vs Normal (Solo TFs)',
  max.overlaps = 30 
)
dev.off()

# MAPA DE CALOR
# Usando solo el top de TFs, no todos los genes
# Dataframe donde los nombres de las filas son los IDs de las muestras
anno_col <- data.frame(Condicion = condition)
rownames(anno_col) <- colnames(vsd) 

#Definir los colores (evita errores)
anno_colors <- list(Condicion = c(Normal = "#1B9E77", GBM = "#D95F02"))

#Heatmap
top_tfs <- rownames(deg_tfs[1:50, ]) 
vsd <- vst(dds, blind = FALSE)
heatmap_data <- assay(vsd)[top_tfs, ]

#Descargar pdf heatmap
pdf("heatmap_top50_TFs_corregido.pdf", width = 12, height = 10)
pheatmap(heatmap_data,
         scale = "row",
         main = "Top 50 Factores de Transcripción Diferenciales",
         annotation_col = anno_col,       # Usamos el nuevo objeto corregido
         annotation_colors = anno_colors, # Colores explícitos
         show_colnames = FALSE,
         labels_row = deg_tfs$Symbol[1:50], 
         fontsize = 10,
         clustering_method = "ward.D2")  # Agrupamiento más robusto
dev.off()

