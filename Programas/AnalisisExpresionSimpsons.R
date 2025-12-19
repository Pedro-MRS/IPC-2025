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

genes <- DataSimpsons[, c("Gene")]
DataSimpsons <- DataSimpsons[, -1]
head(DataSimpsons)

# Redondeamos las estimas de cuantificación de salmon para que DESeq2 pueda interpretarlas como conteos
DataSimpsons <- round(DataSimpsons, digits = 0)

# 3. Preparamos la tabla de metadatos que usará DESeq con las condiciones que queremos que tome.

Simpson_names <- colnames(DataSimpsons)

vector_cond <- c("Obeso1", "Obeso1", "Normopeso", "Normopeso", "Normopeso")

metadatos <- data.frame(row.names = Simpson_names,
                        condition = factor(vector_cond, 
                                           levels = c("Normopeso", "Obeso1")))


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
rownames(res_limpios[which(res_limpios$pvalue < 0.05), ])
# > [1] "NM_001375961.1"   
# NM_001375961.1: cell adhesion molecule 2 CADM2: reducing Cadm2 expression can reverse several traits associated with the metabolic syndrome including obesity, insulin resistance, and impaired glucose homeostasis (DOI: 10.1016/j.molmet.2017.11.010)

# Solo sale un transcrito diferencialmente expresado, así que ampliamos el corte del p-valor a 0.1 para observar transcritos con diferencias marginales de expresión
rownames(res_limpios[which(res_limpios$pvalue < 0.1), ])
# > [1] "NM_001375961.1" "NM_033181.4"    "NM_001003680.3" "XM_005252001.4" "NM_005068.3" "NM_001363705.2"  

# NM_001375961.1: cell adhesion molecule 2 CADM2: reducing Cadm2 expression can reverse several traits associated with the metabolic syndrome including obesity, insulin resistance, and impaired glucose homeostasis (DOI: 10.1016/j.molmet.2017.11.010)
# NM_033181.4: Cannabinoid receptor 1 CNR1: May contribute to the development of diet-induced obesity and several obesity-associated features, such as dyslipidemia and liver steatosis, regulating peripheral lipogenesis, energy expenditure and feeding behavior (DOI: 10.1210/jc.2006-2523)
# NM_001003680.3: Leptin receptor LEPR: Linked to Leptin receptor deficiency (LEPRD), a rare disease characterized by normal levels of serum leptin, hyperphagia and severe obesity from an early age (DOI: 10.1038/32911)
# XM_005252001.4: BDNF/NT-3 growth factors receptor NTRK2: Linked to Obesity, hyperphagia, and developmental delay (OBHD), a disorder characterized by early-onset obesity, hyperphagia, and severe developmental delay in motor function, speech, and language (DOI: 10.1038/nn1336)
# NM_005068.3: Single-minded homolog 1 SIM1: Rare variants in single-minded 1 (SIM1) are associated with severe obesity (DOI: 10.1172/JCI68016)
# NM_001363705.2: E3 ubiquitin-protein ligase UBR2: the risk of a clinical diagnosis of hypertension was increased significantly in carriers of UBR2 and UBR3 [...], but the effect observed in UBR3 Protein Truncating Variants carriers was nearly double that of UBR2(doi: 10.1038/s41588-025-02364-2)






