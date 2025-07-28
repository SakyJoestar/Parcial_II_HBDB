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

# 2. Descargar el archivo suplementario desde GEO
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164073/suppl/GSE164073_Eye_count_matrix.csv.gz"
dest <- "GSE164073_Eye_count_matrix.csv.gz"
download.file(url, destfile = dest)

# 3. Leer la matriz de conteos
counts_df <- read.csv(gzfile(dest), row.names = 1)

countData <- counts_df

# Revisar datos
dim(counts_df)      # genes x muestras
head(counts_df)[,1:5] # primeras filas y columnas

# 4. Definir condiciones experimentales
# ⚠️ Aquí debes mirar la metadata del estudio en GEO
# Vector con las condiciones
conditions <- c(
  rep("control", 9),   # 9 mock
  rep("treatment", 9)  # 9 CoV2
)

colData <- data.frame(
  condition = factor(conditions)
)

rownames(colData) <- colnames(counts_df)

# 5. Crear objeto DESeq2
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)

# 6. Correr análisis diferencial
dds <- DESeq(dds)
res <- results(dds)

# 7. Volcano plot
res$significant <- ifelse(res$padj < 0.05, "Significant", "Not Significant")
volcano <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("red", "blue")) +
  theme_minimal() +
  labs(title = "Volcano Plot - GSE164073",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value")
print(volcano)

# 8. Heatmap con los 20 genes más significativos
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

