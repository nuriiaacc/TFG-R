#Cargar paquetes necesarios
library(TCGAbiolinks)
library(DESeq2)
library(biomaRt)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(survival)
library(survminer)
library(maftools)
library(data.table)
library(enrichR)

#Obtención de datos TCGA

#Tumores primarios
query_tumor <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)
GDCdownload(query_tumor, directory = "GDCdata")
data_tumor <- GDCprepare(query_tumor, directory = "GDCdata")
tumor_counts <- assay(data_tumor, "unstranded")

#Tejido normal
query_normal <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Solid Tissue Normal"
)
GDCdownload(query_normal, directory = "GDCdata")
data_normal <- GDCprepare(query_normal, directory = "GDCdata")
normal_counts <- assay(data_normal, "unstranded")

#Limpiar versiones de Ensembl
rownames(tumor_counts) <- gsub("\\..*", "", rownames(tumor_counts))
rownames(normal_counts) <- gsub("\\..*", "", rownames(normal_counts))

all_counts <- cbind(tumor_counts, normal_counts)

#Obtener TFs
options(timeout = 600)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

tfs_info <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                  filters = 'go',
                  values = 'GO:0003700',
                  mart = ensembl)

tfs_ids <- unique(tfs_info$ensembl_gene_id)

#Análisis de expresión diferencial

#Crear tabla de condiciones
condition <- factor(c(rep("Tumor", ncol(tumor_counts)), 
                      rep("Normal", ncol(normal_counts))), 
                    levels = c("Normal", "Tumor"))

anno_col <- data.frame(Condicion = condition)
rownames(anno_col) <- colnames(all_counts)

#Crear objeto para DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(all_counts),
                              colData = anno_col,
                              design = ~ Condicion)

#Filtrado y ejecución
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

#Resultados generales (no solo TFs)
res <- results(dds, contrast = c("Condicion", "Tumor", "Normal"), pAdjustMethod = "BH")
res_df <- as.data.frame(res)
res_df$ensembl_id <- rownames(res_df)
saveRDS(res_df, file = "resultados_generales_DEA.rds")

#Filtrado de factores de transcripción (TFs)
res_tfs <- res_df %>%
  filter(ensembl_id %in% tfs_ids) %>%
  left_join(tfs_info, by = c("ensembl_id" = "ensembl_gene_id")) %>%
  filter(!is.na(padj)) %>%
  arrange(padj)
res_tfs_ordered <- res_tfs[order(res_tfs$padj), ]

#Guardar resultados
write.csv(res_tfs_ordered, "output/Resultados_DEA_TFs_GBM.csv", row.names = FALSE)

#Visualización de resultados
#Volcano plot
res_tfs_ordered <- res_tfs_ordered %>%
  mutate(Expresión = case_when(
    log2FoldChange > 2 & padj < 0.05 ~ "UP",
    log2FoldChange < -2 & padj < 0.05 ~ "DOWN",
    TRUE ~ "No significativo"
  ))

volcano <- ggplot(res_tfs_ordered, aes(x = log2FoldChange, y = -log10(padj), color = Expresión)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_minimal() +
  geom_text_repel(data = head(res_tfs_ordered, 15), aes(label = hgnc_symbol), 
                  size = 4, color = "black", fontface = "bold") +
  labs(title = "Factores de transcripción diferenciales en GBM",
      x = "Log2 Fold Change", y = "-log10 P-ajustado")

ggsave("output/Volcano_TFs.png", plot = volcano, width = 8, height = 6)


#Mapa de calor (solo el top de TFs)
#Nombres de las filas son los IDs de las muestras
vsd <- vst(dds, blind = FALSE)
anno_col <- data.frame(Condicion = condition)
rownames(anno_col) <- colnames(vsd) 

#Definir los colores
anno_colors <- list(Condicion = c(Normal = "blue", Tumor = "red"))

#Comparar las 5 muestras normales vs. 5 tumorales
set.seed(42)
muestras_grafico <- c(colnames(normal_counts), sample(colnames(tumor_counts), 5))
top50_ids <- res_tfs_ordered$ensembl_id[1:50]
data_equilibrada <- assay(vsd)[top50_ids, muestras_grafico] 
rownames(data_equilibrada) <- res_tfs_ordered$hgnc_symbol[match(rownames(data_equilibrada), res_tfs_ordered$ensembl_id)]

#Crear el mapa de calor
heatmap <- pheatmap(data_equilibrada,
                    cluster_cols = TRUE,
                    scale = "row",
                    main = "Factores de transcripción diferenciales",
                    annotation_col = anno_col[muestras_grafico, , drop=FALSE],       
                    annotation_colors = anno_colors, 
                    show_colnames = FALSE,
                    fontsize = 8,
                    clustering_method = "ward.D2",
                    width = 8, height = 10,
                    filename = "output/Heatmap_TFs.pdf")

#BOXPLOT PARA UN GEN ESPECÍFICO

#Extraer datos de un gen específico
get_gene_data <- function(dds, gene_nm, res_df) {
  #Encontrar ENSEMBL ID a partir de Symbol
  target_id <- res_df$ensembl_id[res_df$hgnc_symbol == gene_nm]
  
  #Extraer conteos normalizados
  d <- plotCounts(dds, gene = target_id, intgroup = "Condicion", returnData = TRUE)
  d$gene <- gene_nm
  return(d)
}

#Elegir dos genes (uno del heatmap vs. uno clásico)
plot_data <- rbind(get_gene_data(dds, "HOXD9", res_tfs),
                   get_gene_data(dds, "SOX2", res_tfs))

#Crear el gráfico
ggplot(plot_data, aes(x = Condicion, y = count, fill = Condicion)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.5) + #Para ver los puntos (normales)
  facet_wrap(~gene, scales = "free") +
  scale_y_log10() + #Escala logarítmica
  theme_minimal() +
  labs(title = "Expresión de TFs Clave", y = "Conteos Normalizados (log10)")

ggsave("output/Boxplot_TFs.png")


#PCA PLOT
PCA_plot <- plotPCA(vsd, intgroup = "Condicion") + 
  theme_minimal() + 
  labs(title = "PCA: Separación Tumor vs Normal")

#Identificar los 2 grupos de PC1
#Extraer las coordenadas del PCA
pca_data <- plotPCA(vsd, intgroup = "Condicion", returnData = TRUE)

ggsave("output/PCA_plot.pdf", plot = PCA_plot, width = 8, height = 6)


#Agrupar los tumores en 2 grupos usando k-means sobre el PC1
set.seed(123)
clusters <- kmeans(pca_data$PC1, centers = 2)
pca_data$subgrupo <- as.factor(clusters$cluster)

#PCA con nuevos subgrupos
ggplot(pca_data, aes(PC1, PC2, color = subgrupo)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Subtipos Identificados en PC1")

#Añadir la información de los subgrupos a los metadatos del objeto dds
#(Solo para las muestras que son Tumor)
dds_tumor <- dds[, dds$Condicion == "Tumor"]
colData(dds_tumor)$subgrupo <- pca_data$subgrupo[pca_data$group == "Tumor"]

#Comparación subgrupo 2 vs subgrupo 1. UP LFC>0 (derecha), DOWN LFC<0 (izquierda)
design(dds_tumor) <- ~ subgrupo
dds_tumor <- DESeq(dds_tumor)
res_subgrupos <- results(dds_tumor, contrast = c("subgrupo", "2", "1"))
res_sub_df <- as.data.frame(res_subgrupos)
res_sub_df$symbol <- tfs_info$hgnc_symbol[match(rownames(res_sub_df), tfs_info$ensembl_gene_id)]

#Ver TFs que más suben ambos subgrupos
#Filtrar los resultados (TFs)
res_tfs_subgrupos <- res_sub_df[rownames(res_sub_df) %in% tfs_info$ensembl_gene_id, ]

#Añadir los nombres a la lista filtrada
res_tfs_subgrupos$symbol <- tfs_info$hgnc_symbol[match(rownames(res_tfs_subgrupos), tfs_info$ensembl_gene_id)]

#TFs predominantes en el subgrupo 2
top_tfs_sub2 <- res_tfs_subgrupos %>%
  filter(padj < 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  head(50)

print(top_tfs_sub2[, c("symbol", "log2FoldChange", "padj")])


#TFs premodminantes en el subgrupo 1
top_tfs_sub1 <- res_tfs_subgrupos %>%
  filter(padj < 0.05) %>%
  arrange(log2FoldChange) %>% #Busca los más negativos
  head(50)

print(top_tfs_sub1[, c("symbol", "log2FoldChange", "padj")])

#Guardar resultados
write.csv(top_tfs_sub1, "output//TFs_Subgrupo1.csv", row.names = TRUE)
write.csv(top_tfs_sub2, "output/TFs_Subgrupo2.csv", row.names = TRUE)


#ANÁLISIS DE ENRIQUECIMIENTO FUNCIONAL
#General: GBM vs. Normal
#FTs sobreexpresados
tfs_generales_up <- res_tfs_ordered %>%
  filter(padj < 0.05 & log2FoldChange > 1.5) %>%
  pull(hgnc_symbol) #Extraemos los nombres

#Convertir a Entrez IDs
entrez_general <- bitr(tfs_generales_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

#Enriquecimiento GO (Procesos biológicos generales)
ego_general <- enrichGO(gene = entrez_general, OrgDb = org.Hs.eg.db, ont = "BP", 
                        pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

#Graficar y guardar
dotplot_general <- dotplot(ego_general, showCategory = 15) + 
  ggtitle("Vías alteradas GBM vs. Normal")
ggsave("output/GO_General_Tumor_vs_Normal.pdf", plot = dotplot_general, width = 8, height = 6)

#Subgrupo 1
genes_sub1 <- res_tfs_subgrupos %>%
  filter(padj < 0.05 & log2FoldChange < -1) %>% 
  pull(symbol) #Extraer los nombres

#Convertir a Entrez
entrez_sub1 <- bitr(genes_sub1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

#Enriquecimiento GO
ego_sub1 <- enrichGO(gene = entrez_sub1, OrgDb = org.Hs.eg.db, ont = "BP", 
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

#Graficar
dotplot_sub1 <- dotplot(ego_sub1, showCategory = 5) + 
  ggtitle("Vías características del Subgrupo 1")
ggsave("output/GO_Subgrupo1.pdf", plot = dotplot_sub1, width = 8, height = 3.8)

#Subgrupo 2
genes_sub2 <- res_tfs_subgrupos %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>% 
  pull(symbol) #Extraer los nombres

#Convertir a Entrez
entrez_sub2 <- bitr(genes_sub2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

#Enriquecimiento GO
ego_sub2 <- enrichGO(gene = entrez_sub2, OrgDb = org.Hs.eg.db, ont = "BP", 
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

#Graficar
dotplot_sub2 <- dotplot(ego_sub2, showCategory = 5) + 
  ggtitle("Vías características del Subgrupo 2")
ggsave("output/GO_Subgrupo2.pdf", plot = dotplot_sub2, width = 8, height = 3.5)


#OTRAS VISUALIZACIONES
#Calcular la similitud
ego_general_sim <- pairwise_termsim(ego_general)
ego_sub2_sim <- pairwise_termsim(ego_sub2)

#Gráfico cnetplot (subgrupo 2)
redgenes <- cnetplot((ego_sub2_sim), 
         showCategory = 5, 
         layout = "kk",          
         color_category = "maroon", 
         size_item = 1,      
         node_label = "all",) +   
  ggtitle("Red genes-conceptos: subgrupo 2")
ggsave("output/Red_genes.pdf", plot = redgenes, width = 8, height = 6)

#Gráfico treeplot
treeplot <- treeplot(ego_general_sim, 
         showCategory = 20,
         nCluster = 5,
         fontsize_tiplab = 2.5,
         fontsize_cladelab = 2.5,
         tiplab_offset = 0.4,
         cladelab_offset = 5)          
  ggtitle("Jerarquía de procesos: subgrupo 2")
ggsave("output/Treeplot.pdf", plot = treeplot, width = 8, height = 6)


#ANÁLISIS DE SUPERVIVENCIA
df_clinico <- as.data.frame(colData(data_tumor))

#Preparar la tabla de supervivencia
#Unir los subgrupos del PCA con los datos clínicos
df_clinico$subgrupo <- pca_data$subgrupo[match(rownames(df_clinico), rownames(pca_data))]

survival_df <- df_clinico %>%
  mutate(
    #(1 = Dead, 0 = Alive/Not Reported)
    event = ifelse(vital_status == "Dead" | vital_status == "1", 1, 0),
    #Definir el tiempo
    time = pmax(as.numeric(days_to_death), 
                as.numeric(days_to_last_follow_up), 
                na.rm = TRUE)
  ) %>%
  filter(!is.na(time) & time > 0 & !is.na(subgrupo))

#Crear el objeto y ajustar el modelo
surv_obj <- Surv(time = survival_df$time, event = survival_df$event)
fit <- survfit(surv_obj ~ subgrupo, data = survival_df)

#Graficar
km_plot <- ggsurvplot(
  fit,
  data = survival_df,
  palette = c("#E7B800", "#2E9FDF"), 
  pval = TRUE,                       
  conf.int = TRUE,                   
  risk.table = TRUE,                 
  legend.labs = c("Subgrupo 1", "Subgrupo 2"),
  xlab = "Tiempo (Días)",
  ylab = "Probabilidad de Supervivencia",
  title = "Análisis de Supervivencia: Subgrupos Identificados en GBM",
  ggtheme = theme_minimal()
)

#Guardar
print(km_plot)
ggsave("output/Supervivencia_GBM_Subgrupos.pdf", plot = km_plot$plot, width = 8, height = 7)


#ANÁLISIS DE MUTACIONES

#Descargar los datos de mutaciones
query_maf <- GDCquery(project = "TCGA-GBM", 
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation",
                      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")

GDCdownload(query_maf, directory = "GDCdata")
maf_data <- GDCprepare(query_maf, directory = "GDCdata")

#Filtrar IDH1
muestras_idh1_mutadas <- maf_data$Tumor_Sample_Barcode[maf_data$Hugo_Symbol == "IDH1"]
maf_data_filtrado <- maf_data[!maf_data$Tumor_Sample_Barcode %in% muestras_idh1_mutadas, ]

#Ajustar IDs
pca_data$short_id <- substr(rownames(pca_data), 1, 16)
pca_data_filtrado <- pca_data[!pca_data$short_id %in% substr(muestras_idh1_mutadas, 1, 16), ]

clinical_maf <- data.table(Tumor_Sample_Barcode = maf_data_filtrado$Tumor_Sample_Barcode,
                           Subgrupo = pca_data_filtrado$Subgrupo[match(substr(maf_data_filtrado$Tumor_Sample_Barcode, 1, 16), pca_data_filtrado$short_id)])

maf_obj <- read.maf(maf = maf_data_filtrado, clinicalData = clinical_maf)

#Mutaciones que ocurren juntas
pdf("output/Interacción_mutaciones.pdf", width = 12, height = 8)
somaticInteractions(maf = maf_obj, top = 25, pvalue = c(0.05, 0.1))
dev.off()

# Analizar rutas oncogénicas conocidas
pdf("output/Rutas_oncogénicas.pdf", width = 12, height = 8)
pathways(maf = maf_obj, plotType = "bar")
dev.off()


#BÚSQUEDA DE DIANAS FARMACÉUTICAS

dianas_potenciales <- res_tfs_ordered %>%
  filter(Expresión == "UP" & padj < 0.001 & log2FoldChange > 2) %>%
  pull(hgnc_symbol)

#Consulta base de datos farmaceúticos
enriched <- enrichr(tfs_generales_up, dbs)

#Resultados
drug_results <- enriched[["DSigDB"]]

if (!is.null(drug_results) && nrow(drug_results) > 0) {
  
  drug_results_clean <- drug_results %>%
    filter(Adjusted.P.value < 0.05) %>%
    arrange(Adjusted.P.value) %>%
    dplyr::select(Term, Overlap, Adjusted.P.value, Genes)}
  
#Gráfico
plot_drugs <- ggplot(head(drug_results_clean, 15), 
                      aes(x = reorder(Term, -log10(Adjusted.P.value)), 
                          y = -log10(Adjusted.P.value), fill = Adjusted.P.value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_viridis_c(direction = -1) +
    theme_minimal() +
    labs(title = "Fármacos Candidatos (DSigDB)",
         subtitle = "Basado en TFs sobreexpresados en GBM",
         x = "Fármaco / Set de Datos", y = "-log10 P-adj") +
    theme(legend.position = "none")

ggsave("output/Dianas_Terapeuticas_Enrichr.pdf", plot = plot_drugs, width = 10, height = 7)
write.csv(drug_results_clean, "output/Dianas_Farmacos_Enrichr.csv", row.names = FALSE)


