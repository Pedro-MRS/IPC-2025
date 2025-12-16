# ===========================================================
#       
#        Análisis de expresión diferencial de genes 
#      relacionados con la obesidad mediante RNA-seq
#
# ===========================================================


# 1. Cargamos las librerías.

library(BiocManager)
library(DESeq2)
library(dplyr)

# 2. Cargamos la matriz de datos indicando el separador (;).

DataSimpsons <- read.table("SimpsonsExpression.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE) 

# Ponemos los ID de los genes como nombres de fila y, posteriormente, eliminamos la fila Gene.
rownames(DataSimpsons) = DataSimpsons$Gene
head(DataSimpsons)

genes = DataSimpsons[, c("Gene")]
DataSimpsons = DataSimpsons[, -1]
head(DataSimpsons)

DataSimpsons <- round(DataSimpsons, digits = 0)

# 3. Preparamos la tabla de metadatos que usará DESeq con las condiciones que queremos que tome.

Simpson_names <- colnames(DataSimpsons)

vector_cond <- c("Obeso1", "Normopeso", "Obeso1", "Normopeso", "Normopeso")

metadatos <- data.frame(row.names = Simpson_names,
                        condition = factor(vector_cond, 
                                           levels = c("Normopeso", "Obeso1")))


# 4. Analizamos los genes con DESeq2.

# Creamos un objeto DESeqDataSet con el conteo y los metadatos, especificando el genotipo como diseño experimental
dds <- DESeqDataSetFromMatrix(countData = DataSimpsons, colData = metadatos, design=~condition)

# Ejecutamos DESeq2
dds_Simpsons <- DESeq(dds)

# Extraemos la información de los resultados específicos comparando las condiciones del diseño
resultado_simpsons = results(dds_Simpsons, contrast=c("condition", "Obeso1", "Normopeso"), alpha=1e-5) 

resultado_simpsons

# Lo guardamos en un data frame para poder ver mejor los resultados y establecemos los genes como filas.

write.csv(as.data.frame(resultado_simpsons), file = "DESeq2_simpson.csv")
DESeq_Simpsons <- read.csv("DESeq2_resultados_simpson.csv",
                           header = TRUE, 
                           row.names = 1)

# También se puede crear otra tabla que elimina los valores N/A y deja unicamente aquellos relevantes para trabajar con ellos.
# 
# Tomamos como referencia para eliminar los N/A el p-valor ajustado (padj) ya que es el único que corrige estadísticamente el riesgo de obtener 
# falsos positivos (errores de Tipo I) al probar la expresión de cientos de genes simultáneamente. 
#
# Entonces, con el código siguiente hacemos que mantenga, usando la columna padj como referencia, solo aquellas filas con un valor numérico válido.

res_limpios <- resultado_simpsons[!is.na(resultado_simpsons$padj), ]

write.csv(as.data.frame(res_limpios), file = "DESeq2_simpsons_limpios.csv")
DESeq_Simpsons_limpios <- read.csv("DESeq2_resultados_limpios.csv",
                            header = TRUE, 
                            row.names = 1)



