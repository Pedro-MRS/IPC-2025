########### Código para generar el Pheatmap de los datos  ################

# Carga de librerías
library(DESeq2)
library(tidyverse)
library(pheatmap)

# Convertimos a data.frame y añadimos nombres de genes
res_limpios_df <- as.data.frame(res_limpios)
res_limpios_df$gene <- rownames(res_limpios_df)

#Secciones de genes expresados diferencialmente
genes_sig <- res_limpios_df %>%
  filter(pvalue < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(pvalue)

genes_sig <- genes_sig[1:50, ]   # top 50 genes
dds <- estimateSizeFactors(dds)

# Extraemos los conteos normalizados del objeto DESeq2
mat_norm <- counts(dds, normalized = TRUE)
mat_norm_df <- as.data.frame(mat_norm)

# Crear columna "Gene" con los nombres de los genes
mat_norm_df$Gene <- rownames(mat_norm_df)

# Opcional: mover la columna "Gene" al principio
mat_norm_df <- mat_norm_df[, c(ncol(mat_norm_df), 1:(ncol(mat_norm_df)-1))]

# Revisar
head(mat_norm_df)

# Subconjunto con genes significativos
mat_heatmap <- mat_norm_df[genes_sig$gene, ]
mat_heatmap <- mat_norm_df[mat_norm_df$Gene %in% genes_sig$gene, -1]  # excluye columna Gene para heatmap

# Mantener los nombres de los genes como rownames
rownames(mat_heatmap) <- mat_norm_df$Gene[mat_norm_df$Gene %in% genes_sig$gene]

# Transformación logarítmica
mat_log <- log2(mat_heatmap + 1)

# Escalado por gen (fila)
mat_scaled <- t(scale(t(mat_log)))

# DEFINIR GRUPOS (CLAVE)
# -------------------------------
grupo <- data.frame(
  Condicion = ifelse(
    colnames(mat_scaled) %in% c("HomerSimpson", "AbrahamSimpson"),
    "Obeso1",
    "Normopeso"
  )
)

rownames(grupo) <- colnames(mat_scaled)

# Ordenar columnas por grupo (Obeso1 primero)
# Convertir a factor con orden deseado
grupo$Condicion <- factor(
  grupo$Condicion,
  levels = c("Obeso1", "Normopeso")
)
orden_cols <- order(grupo$Condicion)

annotation_colors <- list(
  Condicion = c(
    Obeso1 = "black",     # rojo
    Normopeso = "chartreuse"  # azul
  )
)

mat_scaled_ord <- mat_scaled[, orden_cols]
grupo_ord <- grupo[orden_cols, , drop = FALSE]

# -------------------------------
# Heatmap FINAL
# -------------------------------
pheatmap(
  mat_scaled_ord,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = grupo_ord,
  annotation_colors = annotation_colors,
  fontsize_row = 8,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Heatmap de genes diferencialmente expresados\n(Obeso1 vs Normopeso)",
  border_color = NA
)

<img width="963" height="782" alt="pheatmap" src="https://github.com/user-attachments/assets/1119ea13-3e70-492a-8c98-11a15e1e6fda" />



