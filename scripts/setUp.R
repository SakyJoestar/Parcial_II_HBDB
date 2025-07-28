# Instalar BiocManager si no está instalado
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

install.packages(c('xml2', 'curl', 'rentrez', 'SummarizedExperiment', 'rvest', 'httr2'))

# Usar BiocManager para instalar DESeq2 (paquete de Bioconductor)
BiocManager::install("DESeq2")
BiocManager::install("GEOquery")

# Instalar paquetes adicionales desde CRAN
install.packages("ggplot2")
install.packages("pheatmap")


# Cargar los paquetes (si ya están instalados)
library(DESeq2)
library(GEOquery)
library(ggplot2)
library(pheatmap)

#probar
sessionInfo()

