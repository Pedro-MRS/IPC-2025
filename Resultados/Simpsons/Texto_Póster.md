# ACTIVIDAD 2 - SECUENCIACIÓN Y ÓMICAS
## AUTORES
Paula Damián Moral - Paula Sofía Norniella Jamart - Pedro García Saez - Pedro Manuel Rodríguez Santos - Pilar Navarro Sánchez - Rubén Martín Blázquez

## INTRODUCCIÓN
La obesidad es un importante problema de salud pública con un marcado componente genético y familiar. El análisis de la expresión génica mediante RNAseq permite identificar genes diferencialmente expresados asociados a esta condición. En este estudio se analizan las diferencias de expresión génica entre individuos con obesidad y normopeso de la familia Simpson.

Nueva propuesta: La obesidad es un importante problema de salud pública con un marcado componente genético y familiar. El análisis de la expresión génica mediante RNA-seq permite identificar genes diferencialmente expresados asociados a esta condición. En este estudio se tiene como objetivo analizar la expresión diferencial de genes relacionados con la obesidad entre individuos obesos y normopeso de la familia Simpson.

## METODOLOGÍA
La metodología se basó en el análisis de expresión génica diferencial a partir de datos simulados de RNA-seq relacionados con la obesidad. Se realizó un control de calidad inicial de las lecturas mediante FastQC y MultiQC, seguido de la cuantificación y normalización de la expresión génica con Salmon. El análisis de expresión diferencial entre los grupos Obeso 1 y Normopeso se llevó a cabo utilizando DESeq2. Los resultados se visualizaron mediante gráficos volcano y mapas de calor empleando ggplot2, EnhancedVolcano y pheatmap en entornos R y Linux. Finalmente, la interpretación biológica se realizó mediante bases de datos como GeneCards y PubMed, incluyendo de forma opcional un análisis de enriquecimiento funcional. 

Propuesta: En este trabajo se llevó a cabo un análisis de expresión génica diferencial a partir de datos simulados de RNA-seq centrados en genes asociados a la obesidad. Las lecturas crudas se sometieron a un control de calidad mediante FastQC y MultiQC. Posteriormente, la cuantificación de la expresión génica se realizó con Salmon y el análisis diferencial entre los grupos Obeso1 y Normopeso se llevó a cabo utilizando DESeq2. Los resultados se visualizaron mediante gráficos volcano y mapas de calor generados en R empleando las librerías ggplot2, EnhancedVolcano y pheatmap. Finalmente, la interpretación biológica de los genes diferencialmente expresados se apoyó en bases de datos especializadas, como GeneCards y PubMed, incluyendo de forma opcional análisis de enriquecimiento funcional.

## RESULTADOS
Tras realizar el análisis y representaciones gráficas pertinentes, se obtienen los siguientes genes principales:

**PCSK1 (NM_001177875.2 - NM_000439.5)**: Función de codificación para proteasa que interviene en la activación proteolítica de precursores hormonales. Entre los sustratos sobre los que actúa esta proteasa se encuentran la proopiomelanocortina (POMC) y la insulina (1).

**KSR2 (NM_173598.6)**: Favorece la actividad de la proteína serina/treonina quinasa y se asocia con la regulación positiva de la termogénesis, debido a lo cual puede tener un papel en el gasto energético. (2).

**LEP (NM_000230.3)**: Encargado de la codificación de la proteína Leptina, secretada por los adipocitos blancos y fundamental en la regulación del equilibrio energético y el control del peso corporal (3). 

## DISCUSIÓN

Los genes obtenidos tras el análisis de los diferentes perfiles tienen una importante relevancia en la regulación del peso corporal y el metabolismo. Estos procesos dependes de una compleja interacción entre procesamientos hormonales (PCSK1), señalización celular (KSR2) y regulación del apetito (LEP) (1, 2, 3). 
En el caso de PCSK1, activa importantes precursores hormonales, pero su deficiencia severa puede provocar obesidad y trastornos endocrinos. Mutaciones en este gen están relacionadas con una mayor susceptibilidad a la obesidad. También está relacionada con este gen la enfermedad Locus de Rasgo Cuantitativo del Índice de Masa Corporal 12 (BMI QTL 12)  (1, 4). 
Por otro lado, el gen KSR2 coordina la señalización entre AMPK y Raf/MEK/ERK, regulando la termogénesis y el gasto energético. Patologías relacionadas con este gen incluyen el Locus de Rasgo Cuantitativo del Índice de Masa Corporal 11 (BMI QTL 11), relacionándolo con una posible predisposición genética al peso corporal (2, 5). 
Y por último el gen LEP, encargado de la producción de leptina. La leptina activa vías de señalización que inhiben la alimentación y promueven el gasto energético. Niveles elevados de leptina se han asociado a la resistencia de dicha hormona, provocando una desregulación en el ciclo de saciedad y fomentando la hiperfagia y la obesidad. Mutaciones en este gen y regiones reguladoras se asocian a obesidad grave y obesidad mórbida con hipogonadismo. Algunas enfermedades asociadas incluyen el ya mencionado Locus de Rasgo Cuantitativo del Índice de Masa Corporal 11. Sin embargo (3, 6).

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


Alternativa a conclusiones:

**Homer Simpson**: Perfil con alto riesgo metabólico y resistencia a la leptina que deriva en un aumento del apetito y disminución del gasto energético.

**Abraham y Bart Simpson**: Perfil saludable y equilibrado, rápido procesamiento de hormonas debido a expresión de variante de PCSK1.

**Lisa Simpson**: Perfil metabólico muy optimizado, quema de grasa eficaz y regulación del apetito adecuada gracias a una adecuada respuesta hormonal.

**Maggie Simpson**: Perfil metabólico muy optimizado con respuesta hormonal más ágil debido a una mayor expresión en KSR2 y LEP.

En el caso de los pacientes que presentaban normopeso, desde el punto de vista genético se puede explicar debido a la optimización de los ciclos metabólicos mediante KSR2 y LEP (Lisa y Maggie) y también por la expresión de una de las variantes específicas de PCSK1 (Bart). Por el contrario, en los perfiles de obesidad, como es el caso de Homer, pueden deberse a una sobreexposición y posterior resistencia a la leptina. En el caso de Abraham según nuestros hallazgos, que presenta un perfil de expresión igual que Bart (normopeso), no podemos asegurar que su obesidad sea debida exclusivamente a causas genéticas.

Otra conclusion:

Los resultados sugieren una desregulación profunda de la vía de la leptina y el control del apetito. La sobre-expresión del receptor de leptina (LEPR) junto con la sub-expresión de la leptina misma (LEP) y otros genes clave como MC4R en el grupo Obeso1, apunta a una alteración en los mecanismos de señalización metabólica que diferencian claramente ambos estados fenotípicos.

## BIBLIOGRAFÍA

1.  GeneCards. PCSK1 gene [Internet]. [Citado el 19 de diciembre de 2025]. Disponible en: https://www.genecards.org/cgi-bin/carddisp.pl?gene=PCSK1&keywords=NM_001177875.2
2.  GeneCards. KSR2 gene [Internet]. [Citado el 19 de diciembre de 2025]. Disponible en: https://www.genecards.org/cgi-bin/carddisp.pl?gene=KSR2&keywords=NM_173598.6
3.  GeneCards. LEP gene [Internet]. [Citado el 19 de diciembre de 2025]. Disponible en: https://www.genecards.org/cgi-bin/carddisp.pl?gene=LEP&keywords=NM_000230.3
4.  Ramos-Molina B, Martin MG, Lindberg I. PCSK1 Variants and Human Obesity. Prog Mol Biol Transl Sci. 2016;140:47-74. doi: 10.1016/bs.pmbts.2015.12.001.
5.  Guo L, Costanzo-Garvey DL, Smith DR, Neilsen BK, MacDonald RG, Lewis RE. Kinase Suppressor of Ras 2 (KSR2) expression in the brain regulates energy balance and glucose homeostasis. Mol Metab. 2016 Dec 18;6(2):194-205. doi: 10.1016/j.molmet.2016.12.004.
6.  Zhang Y, Scarpace PJ. The role of leptin in leptin resistance and obesity. Physiol Behav. 2006 Jun 30;88(3):249-56. doi: 10.1016/j.physbeh.2006.05.038.

