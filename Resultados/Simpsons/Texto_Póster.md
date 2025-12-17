# ACTIVIDAD 2 - SECUENCIACIÓN Y ÓMICAS
## AUTORES
Paula Damián Moral - Paula Sofía Norniella Jamart - Pedro García Saez - Pedro Manuel Rodríguez Santos - Pilar Navarro Sánchez - Rubén Martín Blázquez

## INTRODUCCIÓN

## METODOLOGÍA
La metodología consistió en el análisis de expresión diferencial de genes relacionados con la obesidad de datos de RNA-seq simulados, incluyendo un control de calidad inicial de las lecturas con FastQC y MultiQC. 
La cuantificación y normalización de la expresión génica se realizó mediante Salmon, generando una matriz de conteos por gen. 
La expresión diferencial, que se realizo entre los grupos Obeso 1 y Normopeso, se evaluó utilizando DESeq2. 
Los resultados se visualizaron mediante gráficos de tipo volcano y mapas de calor (generados con ggplot2, EnhancedVolcano y pheatmap) empleando los entornos R y Linux. 
Los resultados se interpretaron mediante bases de datos de función génica como GeneCards y PubMed y se empleo un analisis de enriquecimiento funcional (opcional).

## RESULTADOS
**PCSK1 (NM_001177875.2 - NM_000439.5)**: Función de codificación para proteasa que interviene en la activación proteolítica de precursores hormonales. Ente los sustratos sobre los que actúa esta proteasa se encuentran la proopiomelanocortina (POMC) y la insulina. Mutaciones en este gen están relacionadas con una mayor susceptibilidad a la obesidad. También está relacionada con este gen la enfermedad Locus de Rasgo Cuantitativo del Índice de Masa Corporal 12 (BMI QTL 12).

**KSR2 (NM_173598.6)**: Favorece la actividad de la proteína serina/treonina quinasa y se asocia con la regulación positiva de la termogénesis, debido a lo cual puede tener un papel en el gasto energético. Patologías relacionadas con este gen incluyen el Locus de Rasgo Cuantitativo del Índice de Masa Corporal 11 (BMI QTL 11), relacionándolo con una posible predisposición genética al peso corporal.

**LEP (NM_000230.3)**: Encargado de la codificación de la proteína Leptina, secretada por los adipocitos blancos y fundamental en la regulación del equilibrio energético y el control del peso corporal. La leptina activa vías de señalización que inhiben la alimentación y promueven el gasto energético. Mutaciones en este gen y regiones reguladoras se asocian a obesidad grave y obesidad mórbida con hipogonadismo. Algunas enfermedades asociadas incluyen el ya mencionado Locus de Rasgo Cuantitativo del Índice de Masa Corporal 11.

## DISCUSIÓN

## CONCLUSIONES

## BIBLIOGRAFÍA
