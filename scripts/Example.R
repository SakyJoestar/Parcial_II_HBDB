### 2. Simulación de datos de expresión génica ----
set.seed(123)

# Matriz de cuentas (genes x muestras)
counts <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)

# Condiciones experimentales: 5 control, 5 tratamiento
colData <- data.frame(
  condition = factor(rep(c("control", "treatment"), each = 5))
)

# Nombres para genes y muestras
rownames(counts) <- paste0("gene", 1:100)
colnames(counts) <- paste0("sample", 1:10)

### 3. Crear el objeto DESeqDataSet ----
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = colData,
  design = ~ condition
)

### 4. Análisis diferencial de expresión ----
dds <- DESeq(dds)
res <- results(dds)

### 5. Visualización de resultados ----

## Gráfico de volcanes
res$significant <- ifelse(res$padj < 0.05, "Significant", "Not Significant")

volcano <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("red", "blue")) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value")

print(volcano)

## Mapa de calor
# Top 20 genes con menor p-valor ajustado
top_genes <- head(order(res$padj, na.last = NA), 20)
top_counts <- counts[top_genes, ]

# Normalización con transformación de variancia estabilizada
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
heatmap_data <- assay(vsd)[top_genes, ]

# Asegura que las columnas coincidan
annotation_col <- colData
rownames(annotation_col) <- colnames(heatmap_data)

# Graficar mapa de calor
pheatmap(heatmap_data,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = colData)

### 6. Interpretación (comentarios):
# - Volcano Plot: los genes en rojo están diferencialmente expresados con p-valor ajustado < 0.05
# - Mapa de calor: muestra los niveles de expresión normalizados de los 20 genes más significativos

#########################################
# Fin del script 🎉
#########################################