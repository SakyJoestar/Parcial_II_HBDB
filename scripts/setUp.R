# Instalar BiocManager si no está instalado
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Usar BiocManager para instalar DESeq2 (paquete de Bioconductor)
BiocManager::install("DESeq2")

# Instalar paquetes adicionales desde CRAN
install.packages("ggplot2")
install.packages("pheatmap")

# Cargar los paquetes (si ya están instalados)
library(DESeq2)
library(ggplot2)
library(pheatmap)

