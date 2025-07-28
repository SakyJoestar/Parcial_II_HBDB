#########################################
# Análisis de Expresión Génica con datos reales (GEO)
# Samuel | julio 2025
#########################################

### 1. Instalación y carga de paquetes ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("GEOquery")

library(DESeq2)
library(GEOquery)
library(ggplot2)
library(pheatmap)

### 2. Descargar los datos desde GEO ----
# Carga la plataforma GSE109379
gse <- getGEO("GSE109379", GSEMatrix = TRUE)
gse <- gse[[1]]  # accedemos al primer objeto de la lista

# Extrae los datos de expresión y metadatos
exp_data <- exprs(gse)  # matriz de expresión
meta_data <- pData(gse)  # información de las muestras

# Verifica las condiciones experimentales
table(meta_data$source_name_ch1)

# Creamos un factor para condiciones (LPS vs control)
condition <- ifelse(grepl("LPS", meta_data$source_name_ch1), "treatment", "control")
colData <- data.frame(condition = factor(condition))
rownames(colData) <- colnames(exp_data)

### 3. Crear objeto DESeq2 ----
dds <- DESeqDataSetFromMatrix(countData = round(exp_data),
                              colData = colData,
                              design = ~ condition)

### 4. Análisis diferencial ----
dds <- DESeq(dds)
res <- results(dds)

### 5. Visualización ----

## Volcano plot
res$significant <- ifelse(res$padj < 0.05, "Significant", "Not Significant")

volcano <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("red", "blue")) +
  theme_minimal() +
  labs(title = "Volcano Plot - GSE109379",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value")

print(volcano)

## Heatmap con genes más significativos
top_genes <- head(order(res$padj, na.last = NA), 20)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
heatmap_data <- assay(vsd)[top_genes, ]

annotation_col <- colData
rownames(annotation_col) <- colnames(heatmap_data)

pheatmap(heatmap_data,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = annotation_col)

#########################################
# Fin del análisis real 🔬
#########################################