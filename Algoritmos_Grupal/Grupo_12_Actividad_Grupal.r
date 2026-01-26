##############################################################################
##############################################################################
###                                                                        ###
### ACTIVIDAD GRUPAL 3. ANALISIS CONJUNTO DATOS ORIGEN BIOLOGICO MEDIANTE  ###
### MEDIANTE TECNICAS DE MACHINE LEARNING SUPERVISADAS Y NO SUPERVISADAS   ###
###                                                                        ###
##############################################################################
##############################################################################


# Se establece el directorio de trabajo.

setwd("C:/Users/laura/Documents/Unir/Algoritmos/ACTIVIDAD GRUPAL")

# Carga de librerias

library(tidyverse)   # Libreria analisis de datos
library(dplyr)       # Libreria para hacer la representación gráfica
library(stats)       # Libreria para el PCA
library(plotly)      # Libreria para MDS en 3D
library(Rtsne)       # Libreria para el t-SNE
library(uwot)        # Libreria para UMAP
library(RDRToolbox)  # Libreria para Isomap
library (factoextra) # Libreria para clustering k-means
library(ggplot2)     # Libreria para hacer la representación gráfica

# Generar dataframe

genes <- readLines("column_names.txt")  # Leer los nombres de los genes (columnas)
expresion <- read.csv("gene_expression.csv", sep = ";", header = FALSE, check.names = FALSE) # Leer la matriz de expresión
colnames(expresion) <- genes            # Asignar los nombres de los genes a las columnas
clases <- read.csv("classes.csv", sep = ";", header = FALSE)   # Leer las clases e IDs de las muestras
colnames(clases) <- c("ID", "clase")
df <- cbind(clases, expresion)     #  Combinar todo en un solo DataFrame
rownames(df) <- df$ID # Establecer el ID de la muestra como nombre de fila

write_csv(df, "expresion_genica.csv")     # Se guarda el dataframe como .csv

# Se evalua la existencia de NAS
anyNA(df)
na_counts <- colSums((is.na(df)))
na_counts

# Se comprueba si hay valores 0
any(df == 0)
# Se han encontrado valores 0, que pueden ser por errores en la toma de datos o porque son genes
# que no se han expresado, y, debido a esta ultima circunstancia, se cuenta con ellos
# para realizar el análisis

zero_counts <- colSums(df == 0)
zero_counts
# Graficamos los estadísticos

zero_df <- data.frame(
  Variable = names(zero_counts),
  Zeros = as.numeric(zero_counts)
)
ggplot(zero_df, aes(x = Variable, y = Zeros, fill = Variable)) +
  geom_bar(stat = "identity") +
  labs(title = "Cantidad de ceros por columna",
       x = "Variable",
       y = "Número de ceros") +
  theme_minimal() +
  theme(legend.position = "none")  # Oculta la leyenda

# Se calcula la varianza de cada gen(excluyendo las dos primeras columnas)
varianzas <- apply(df[, -c(1, 2)], 2, var)

# Identificacion de cuantos genes tienen varianza cero
var_genes <- names(varianzas[varianzas == 0])
cat("Se han encontrado", length(var_genes), "genes con varianza cero.")
print(head(var_genes, 10))
# Se eliminan los genes con varianza 0 y se genera un nuevo dataframe sin esos genes
# Esto es necesario para poder realizar metodos de reduccion de dimensionalidad, ya que los genes se van a ordenar en función de su varianza 
# Mantenemos las columnas 1 y 2 (ID y Clase) y solo los genes con varianza > 0
# Creamos el dataframe final manteniendo solo genes con varianza > 0
# Mantenemos columnas 1 y 2 (ID y clase) y filtramos el resto
df_final <- df[, c(colnames(df)[1:2], names(varianzas[varianzas > 0]))]

# Se comprueba la eliminacion de los genes con varianza 0

cat("Columnas antes:", ncol(df), "\n")
cat("Columnas ahora:", ncol(df_final), "\n")

# Se realiza un diagrama de cajas para para las variables y vemos los estadísticos y outliers

boxplot(df[, 3:13], main = "Boxplot de los 10 primeros genes")

#################################################################
#################################################################
##            MÉTODOS DE APRENDIZAJE NO SUPERVISADO            ##
#################################################################
#################################################################


#################################################################
##          0. MÉTODO REDUCCIÓN DIMENSIONALIDAD PCA            ##
#################################################################

### Funcion prcomp()
#   df_final: conjunto de datos
#   center: si queremos que las variables esten centradas en cero
#   scale: si queremos que las variables tengan varianza 1

# Calculo de componentes principales con a funcion prcomp
# Se ejecuta PCA excluyendo las columnas 1 y 2 que corresponden a ID y clase
pca.results <- prcomp(df_final[, -c(1, 2)], center = TRUE, scale. = TRUE)

# Ver resumen de PCA
summary(pca.results)

# Resultado de las componentes principales
pca.df <- data.frame(pca.results$x)

# Varianza de los genes en los componenetes principales (cuadrado de la desviacion tipica)
varianzas <- pca.results$sdev^2

# Total de la varianza de los datos
total.varianza <- sum(varianzas)

# Varianza explicada por cada componente principal
varianza.explicada <- varianzas/total.varianza

# Se calcula la varianza acumulada
varianza.acumulada <- cumsum(varianza.explicada)

# Se toma el numero de componentes principales que explican el 90% de la varianza
n.pc <- min(which(varianza.acumulada > 0.90))

# Etiquetas de los ejes del gráfico
x_label <- paste0(paste('PC1', round(varianza.explicada[1] * 100, 2)), '%')
y_label <- paste0(paste('PC2', round(varianza.explicada[2] * 100, 2)), '%')

# Representación gráfica de las primeras dos componentes principales respecto a los datos

ggplot(pca.df, aes(x=PC1, y=PC2, color=df_final$clase)) +
  geom_point(size=3) +
  scale_color_manual(values=c('red', 'blue', 'green', 'orange', 'purple')) +
  labs(title='PCA - Clases tumores hepáticos', x=x_label, y=y_label, color='Grupo') +
  theme_classic() +
  theme(panel.grid.major = element_line(color="gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title = element_text(hjust = 0.5))



#################################################################
##          1. MÉTODO REDUCCIÓN DIMENSIONALIDAD MDS            ##
#################################################################

# Funcion cmdscale()
#   d: matriz de distancias (usaremos la funcion dist)
#   k: numero que indica el tamaño final de los datos (max num de variables)
#   eig: si calculamos autovalores de las variables. Nos sirve para el calculo 
#        de la varianza explicada, es decir, para coger las columnas de mayor
#        variabilidad
#   x.ret: para devolver la matriz de distancias que calcule el algoritmo

#   points: dataframe de tamaño k que representa las nuevas coordenadas
#   eig: vector con los autovalores para elegir el numero de dimensiones


# Se utiliza la funcion dist para calcular la matriz de distancias euclideas
# matriz NxN de distancias entre todos los puntos

distances <- dist(df_final, method = 'euclidean')

# Utilizamos la funcon cmdscale para realizar el MSD

mds.results <- cmdscale(distances, eig=TRUE, k=2, x.ret=TRUE)

# Calculamos la varianza explicada

varianza.explicada <- mds.results$eig/sum(mds.results$eig) * 100

# Sacamos en un dataframe los puntos del mds

mds.df <- data.frame(mds.results$points)


# Grafico

ggplot(mds.df, aes(x=X1, y=X2, color=df_final$clase)) +
  geom_point(size=3) + 
  scale_color_manual(values=c("red", "blue", "green", "orange", "purple")) +
  labs(title="MDS - Clases tumores hepáticos", x="Dimension 1 (X1)", y="Dimension 2 (X2)",color = "Grupo") +
  theme_classic() + 
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))

###############################
# Grafico varianza n-dim MDS  #
###############################

# Lectura de datos
datos_raw <- read.csv("expresion_genica.csv")

# Guardado en un dataframe
data <- sapply(datos_raw[, 3:ncol(datos_raw)], as.numeric)
distancia <- dist(data)

# Aplicar MDS a 3 dimensiones

mds_result <- cmdscale(distancia, k = 3)

# Convertir el resultado a un data frame

mds_df <- as.data.frame(mds_result)
colnames(mds_df) <- c("Dim1", "Dim2", "Dim3")

# Añadir etiquetas de clase de tumores
mds_df$clase_tumores <- df_final$clase



#################################################################
##          2. MÉTODO REDUCCIÓN DIMENSIONALIDAD t-SNE          ##
#################################################################

# Funcion Rtsne()
#   X: datos sobre los que reduciremos la dimensionalidad
#   dims: tamaño final del conjunto de datos (mejor <=3) por eficiencia
#   num_threads: hilos a utilizar (procesadores). No hace falta usarlo de momento
# 
#   Variable Y con matriz del t-SNE

tsne <- Rtsne(X=df_final, dims=2)
tsne_result <- data.frame(tsne$Y)

# Graficamos
tsne_graf <- ggplot(tsne_result, aes(x = X1, y = X2, color = df_final$clase)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "t-SNE - Clases tumores hepáticos", x = "Dim 1", y = "Dim 2", color = "Grupo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))

# Mostrar el gráfico
tsne_graf



#################################################################
##          3. MÉTODO REDUCCIÓN DIMENSIONALIDAD UMAP           ##
#################################################################

# Funcion umap()
#     x: reduccion de la dimensionalidad
#     n_neighbours: entero que indica el numero de vecinos cercanos
#     n_componentes: entero que determina el tamaño del espacio de salida
#     metric: define la distancia entre puntos
#     min_dist: distancia minima permitida entre puntos
#     local_conectivity: cuantos puntos están conectados "fuertemente"
#     scale: tipos de escalado
#     ret_model: para sacar los resultados del embedding
#     verbose: tiempo hasta que se complete el calculo

#     Y resultado

umap.results <- umap(df_final, n_neighbors=0.2 * nrow(data),
                     n_components = 2, min_dist = 0.1, local_connectivity=1, 
                     ret_model = TRUE, verbose = TRUE)

umap.df <- data.frame(umap.results$embedding)

m_dist <- dist(umap.df)

# Graficamos

ggplot(umap.df, aes(x = X1, y = X2, color = df_final$clase)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "UMAP - Clases tumores hepáticos", x = "X1", y = "X2", color = "Grupo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))



#################################################################
##          4. MÉTODO REDUCCIÓN DIMENSIONALIDAD ISOMAP         ##
#################################################################

# Funcion Isomap()
#   data -> datos (matriz) sobre los que haremos reduccion de dimensionalidad
#   dim -> dimensiones de las columnas del espacio reducido
#   k -> numero de vecinos cercanos a cada punto. A mayor k mayor computacion
#   potResiduals -> devuelve la varianza explicada por las diferentes dimensiones

#   Si se ha elegido una unica dimension devuelve una matriz
#   Si se ha elegido un vector de dimensiones devolvera una matriz por cada elemento del vector

# Lectura de datos
datos_raw <- read.csv("expresion_genica.csv")

# Guardado en un dataframe
data <- sapply(datos_raw[, 3:ncol(datos_raw)], as.numeric)

# Calculamos isomap de 1 a 10 dimensiones y con 5 vecinos (y 30?)
isomap.results = Isomap(data=data, dims=1:10, k=15, plotResiduals=TRUE)

# Dataframe con los puntos que queremos dibujar en el plano 2D
#     (elegiriamos otro si queremos otra dimension)
isomap.df <- data.frame(isomap.results$dim2) 

# Graficamos
ggplot(isomap.df, aes(x = X1, y = X2, color = df_final$clase)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Isomap - Clases tumores hepáticos", x = 'Dim 1', y = 'Dim 2', color = "Grupo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))



#################################################################
##          5. MÉTODO CLUSTERING NO JERÁRQUICO K-MEANS         ##
#################################################################

# Seleccionar desde la columna 3 hasta la última columna que corresponde a los genes
df_kmeans <- df_final %>% 
  dplyr::select(3:last_col())
df_final_scale <- scale(df_kmeans)  # Normalización z-score

#### Clustering no jerarquico: K-means

# funcion kmeans
#   x: matriz numerica o dataframe
#   centers: numero de clusteres que se desean formar
#   iter.max: numero de iteraciones para el algoritmo
#   nstart: veces que se ejecuta el algoritmo con diferentes centroides (10 - 25)
#   algorithm: varias opciones (Hartingan-wong por defecto)


# k=3
kmeans.result <- kmeans(df_final_scale, centers = 3, iter.max = 100, nstart = 25)

# Visualizacion
fviz_cluster(kmeans.result, df_final_scale, 
             xlab = '', ylab = '',
             geom = "point") + # <--- Esta línea quita las etiquetas
  ggtitle("Cluster plot, centers = 5", subtitle = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(b = -10)))

#################################################################
##          5. MÉTODO CLUSTERING JERARQUICO                    ##
#################################################################

# Realizamos el análisis cluster para agrupar las muestras en racimos (clusters): 
# primero calculando la matriz de distancias euclídeas con dist(method = "euclidean"),
# después se genera el dendrograma con el método de Ward's D con hclust(method = "ward.D")
hclust_model <- df_final %>%
                dist(method = "euclidean") %>% 
                hclust(method = "ward.D")

# Generamos paleta de colores para las clases
colors <- rainbow(5)

# Representamos el dendrograma gráficamente
fviz_dend(hclust_model,
          k = 5,
          palette = colors,
          main = "Dist. Euclídea, Ward's D",
          xlab = "Muestras",
          ylab = "Distancia",
          color_labesls_by_k = TRUE) +
  theme_classic()


#################################################################
#################################################################
##              MÉTODOS DE APRENDIZAJE SUPERVISADO             ##
#################################################################
#################################################################

