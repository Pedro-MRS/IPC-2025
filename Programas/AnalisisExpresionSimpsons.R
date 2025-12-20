# ===========================================================
#       
#        Análisis de expresión diferencial de genes 
#      relacionados con la obesidad mediante RNA-seq
#
# ===========================================================


# 1. Cargamos las librerías.

library(BiocManager)
library(tximport)
library(DESeq2)
library(dplyr)

# 2. Cargamos la matriz de datos indicando el separador (;).

design <- read.csv("Design.csv", sep = ",")
samples <- design$Sample[which(design$Condition == "Normopeso" | design$Condition == "Sobrepeso/Obeso1")]
dir <- paste0(getwd())
files <- c()
for (i in 1:length(samples)) {
  files[i] <- paste0(dir, "/Resultados/Salmon/", samples[i] , "_quant/quant.sf")
}

tx2gene <- read.csv("Transcrito_a_Gen.tsv", sep = "\t")

data.tx <- tximport(files, type = "salmon", tx2gene = tx2gene)
DataSimpsons <- round(data.tx$counts, digits = 0)
colnames(DataSimpsons) <- samples

# 3. Preparamos la tabla de metadatos que usará DESeq con las condiciones que queremos que tome.

Simpson_names <- colnames(DataSimpsons)

vector_cond <- design$Condition[which(design$Condition == "Normopeso" | design$Condition == "Sobrepeso/Obeso1")]

metadatos <- data.frame(row.names = Simpson_names,
                        condition = factor(vector_cond))

metadatos$condition <- as.factor(metadatos$condition)

levels(metadatos$condition) <- c("Normopeso", "Obeso1")

# 4. Analizamos los genes con DESeq2.

# Creamos un objeto DESeqDataSet con el conteo y los metadatos, especificando el genotipo como diseño experimental
dds <- DESeqDataSetFromMatrix(countData = DataSimpsons, colData = metadatos, design=~condition)

# Ejecutamos DESeq2
dds_Simpsons <- DESeq(dds)

# Extraemos la información de los resultados específicos comparando las condiciones del diseño
resultado_simpsons <- results(dds_Simpsons, contrast=c("condition", "Obeso1", "Normopeso"), alpha=1e-5) 

resultado_simpsons

# Lo guardamos en un data frame para poder ver mejor los resultados y establecemos los genes como filas.

DESeq_Simpsons <- as.data.frame(resultado_simpsons)

# También se puede crear otra tabla que elimina los valores N/A y deja unicamente aquellos relevantes para trabajar con ellos.
# 
# Tomamos como referencia para eliminar los N/A el p-valor ajustado (padj) ya que es el único que corrige estadísticamente el riesgo de obtener 
# falsos positivos (errores de Tipo I) al probar la expresión de cientos de genes simultáneamente. 
#
# Entonces, con el código siguiente hacemos que mantenga, usando la columna padj como referencia, solo aquellas filas con un valor numérico válido.

res_limpios <- resultado_simpsons[!is.na(resultado_simpsons$padj), ]

DESeq_Simpsons_limpios <- as.data.frame(res_limpios)

# Comprobamos qué transcritos están diferencialmente expresados (p-valor menor a 0.05)
# Comprobamos qué transcritos están diferencialmente expresados (p-valor menor a 0.05)
rownames(res_limpios[which(res_limpios$padj < 0.05), ])
# [1] "BDNF"  "CADM2" "LEP"   "LEPR"  "MC4R"  "NTRK2" "PCSK1" "POMC"  "SH2B1"

## Boxplots de los genes diferencialmente expresados

degs <- rownames(res_limpios[which(res_limpios$padj < 0.05), ])

counts.norm <- counts(dds_Simpsons, normalize = TRUE)
counts.norm <- counts.norm[degs, ]
counts.norm.t <- as.data.frame(t(counts.norm))
counts.norm.t$Condicion <- as.character(metadatos$condition)

library(reshape2)
counts.norm.melt <- melt(counts.norm.t, value.name = "Expresión")
counts.norm.melt$Genes <- counts.norm.melt$variable
counts.norm.melt$variable <- NULL

ggplot(counts.norm.melt, aes(x = Genes, y = Expresión, colour = Condicion)) + 
  geom_boxplot() +
  theme_classic() +
  scale_colour_manual(values = c(Normopeso = "chartreuse", 
                                 Obeso1 = "black"))

