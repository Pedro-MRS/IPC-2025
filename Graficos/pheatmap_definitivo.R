########### Código para generar el Pheatmap de los datos  ################

# Carga de librerías
library(DESeq2)
library(tidyverse)
library(pheatmap)

# Convertimos a data.frame y añadimos nombres de genes
res_limpios_df <- as.data.frame(res_limpios)

#Secciones de genes expresados diferencialmente
genes_sig <- res_limpios_df %>%
  filter(pvalue < 0.1, abs(log2FoldChange) > 1) %>%
  arrange(pvalue)

# Calculamos los factores de tamaño (diferencias de tamaño relativo entre librerías) 
dds <- estimateSizeFactors(dds)

# Extraemos los conteos normalizados del objeto DESeq2
mat_norm <- counts(dds, normalized = TRUE)
mat_norm_df <- as.data.frame(mat_norm)

# Revisar
head(mat_norm_df)

# Subconjunto con genes significativos
mat_heatmap <- mat_norm_df[rownames(mat_norm_df) %in% rownames(genes_sig), ]

# DEFINIR GRUPOS (CLAVE)
# -------------------------------
grupo <- data.frame(
  Condicion = ifelse(
    colnames(mat_heatmap) %in% c("HomerSimpson", "AbrahamSimpson"),
    "Obeso1",
    "Normopeso"
  )
)

rownames(grupo) <- colnames(mat_heatmap)

annotation_colors <- list(
  Condicion = c(
    Obeso1 = "black",
    Normopeso = "chartreuse"
  )
)


# -------------------------------
# Heatmap FINAL
# -------------------------------
pheatmap(
  mat_heatmap,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = grupo,
  annotation_colors = annotation_colors,
  fontsize_row = 8,
  color = colorRampPalette(c("lightblue", "blue", "darkblue"))(100),
  main = "Heatmap de genes diferencialmente expresados\n(Obeso1 vs Normopeso)",
  border_color = NA
)
