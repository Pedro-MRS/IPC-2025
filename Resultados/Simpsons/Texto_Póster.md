# ACTIVIDAD 2 - SECUENCIACIÓN Y ÓMICAS
## AUTORES
Paula Damián Moral - Paula Sofía Norniella Jamart - Pedro García Saez - Pedro Manuel Rodríguez Santos - Pilar Navarro Sánchez - Rubén Martín Blázquez

## INTRODUCCIÓN
La obesidad es un importante problema de salud pública con un marcado componente genético y familiar. El análisis de la expresión génica mediante RNAseq permite identificar genes diferencialmente expresados asociados a esta condición. En este estudio se analizan las diferencias de expresión génica entre individuos con obesidad y normopeso de la familia Simpson.

## METODOLOGÍA
La metodología se basó en el análisis de expresión génica diferencial a partir de datos simulados de RNA-seq relacionados con la obesidad. Se realizó un control de calidad inicial de las lecturas mediante FastQC y MultiQC, seguido de la cuantificación y normalización de la expresión génica con Salmon. El análisis de expresión diferencial entre los grupos Obeso 1 y Normopeso se llevó a cabo utilizando DESeq2. Los resultados se visualizaron mediante gráficos volcano y mapas de calor empleando ggplot2, EnhancedVolcano y pheatmap en entornos R y Linux. Finalmente, la interpretación biológica se realizó mediante bases de datos como GeneCards y PubMed, incluyendo de forma opcional un análisis de enriquecimiento funcional. 

## RESULTADOS
**PCSK1 (NM_001177875.2 - NM_000439.5)**: Función de codificación para proteasa que interviene en la activación proteolítica de precursores hormonales. Ente los sustratos sobre los que actúa esta proteasa se encuentran la proopiomelanocortina (POMC) y la insulina. Mutaciones en este gen están relacionadas con una mayor susceptibilidad a la obesidad. También está relacionada con este gen la enfermedad Locus de Rasgo Cuantitativo del Índice de Masa Corporal 12 (BMI QTL 12).

**KSR2 (NM_173598.6)**: Favorece la actividad de la proteína serina/treonina quinasa y se asocia con la regulación positiva de la termogénesis, debido a lo cual puede tener un papel en el gasto energético. Patologías relacionadas con este gen incluyen el Locus de Rasgo Cuantitativo del Índice de Masa Corporal 11 (BMI QTL 11), relacionándolo con una posible predisposición genética al peso corporal.

**LEP (NM_000230.3)**: Encargado de la codificación de la proteína Leptina, secretada por los adipocitos blancos y fundamental en la regulación del equilibrio energético y el control del peso corporal. La leptina activa vías de señalización que inhiben la alimentación y promueven el gasto energético. Mutaciones en este gen y regiones reguladoras se asocian a obesidad grave y obesidad mórbida con hipogonadismo. Algunas enfermedades asociadas incluyen el ya mencionado Locus de Rasgo Cuantitativo del Índice de Masa Corporal 11.

## DISCUSIÓN

## CONCLUSIONES
***Abraham Simpson***: sobreexpresión gen PCSK1 en isoforma NM_001177875.2 y normalidad en isoforma NM_000439.5, normalidad en gen KSR2, normalidad en gen LEP.
Es un perfil mayoritariamente saludable desde el punto de vista genético-metabólico, con una particularidad en la expresión de una variante de PCSK1 que, por sí sola, no define una patología clara en ausencia de síntomas clínicos. Procesamiento hormonal acelerado o respuesta compensatoria de otros desequilibrios metabólicos.

***Homer Simpson***: normalidad gen PCSK1 en isoforma NM_001177875.2 y sobreexpresión en isoforma NM_000439.5, sobreexpresion en gen KSR2 pero menos, sobreexpresión en gen LEP pero menos.
Aumento de la expresión de la isoforma principal del gen PCSK1 con aumento en la producción de insulina con un proceso metabólico a la máxima capacidad. Combinado con la sobreexpresión de KSR2 indica que es un metabolismo que aumenta la quema de grasas. La sobreexpresión de LEP el cerebro deja de responder a la hormona y suele derivar en resistencia a la leptina y además, para compensar la falta de respuesta cerebral, las células grasas producen aún más leptina, agravando la resistencia, por lo que el cuerpo actúa como si estuviera en estado de inanición, aumentando el apetito y reduciendo el metabolismo basal. Alto riesgo de enfermedades cardiovasculares.

***Bart Simpson***: sobreexpresión gen PCSK1 en isoforma NM_001177875.2 y normalidad en isoforma NM_000439.5, normalidad en gen KSR2, normalidad en gen LEP.
Es un perfil mayoritariamente saludable desde el punto de vista genético-metabólico, con una particularidad en la expresión de una variante de PCSK1 que, por sí sola, no define una patología clara en ausencia de síntomas clínicos. Procesamiento hormonal acelerado o respuesta compensatoria de otros desequilibrios metabólicos.

***Lisa Simpson***: normalidad gen PCSK1 en isoforma NM_001177875.2 y sobreexpresión en isoforma NM_000439.5 pero menos, sobreexpresion en gen KSR2, sobreexpresión en gen LEP pero menos.
Indica maquinaria metabólica optimizada. Es el patrón opuesto a los trastornos metabólicos lentos, caracterizado por una respuesta hormonal ágil y un consumo energético celular elevado. El gen PCSK1 mantiene la estabilidad del procesamiento hormonal, la sobreexpresión de KSR2 facilita que el cuerpo queme ácidos grasos y utilice la glucosa en los músculos de forma muy eficaz, y una sobreexpresión de LEP moderada provoca una reducción del apetito y favorece mantener el peso.

***Maggie Simpson***: normalidad gen PCSK1 en isoforma NM_001177875.2 y sobreexpresión en isoforma NM_000439.5 pero menos, sobreexpresion en gen KSR2 pero menos, sobreexpresión en gen LEP pero menos.
Como Lisa pero con niveles un poco mas altos en la sobreexpresión de KSR2 y LEP.

## BIBLIOGRAFÍA
