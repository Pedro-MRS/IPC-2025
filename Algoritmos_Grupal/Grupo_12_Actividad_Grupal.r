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
library(ggplot2)     # Libreria para hacer la representación gráfica
library(stats)       # Libreria para el PCA
library(plotly)      # Libreria para MDS en 3D
library(Rtsne)       # Libreria para el t-SNE
library(uwot)        # Libreria para UMAP
library(RDRToolbox)  # Libreria para Isomap
library(factoextra)  # Libreria para clustering k-means
library(caret)       # Librería para usar funciones de partición de datos de entrenamiento y test, y acceder a los modelos de ML
library(kernlab)     # Librería que es dependencia del SVM de caret

# Generar dataframe

genes <- readLines("column_names.txt")  # Leer los nombres de los genes (columnas)
expresion <- read.csv("gene_expression.csv", sep = ";", header = FALSE, check.names = FALSE) # Leer la matriz de expresión
colnames(expresion) <- genes            # Asignar los nombres de los genes a las columnas
clases <- read.csv("classes.csv", sep = ";", header = FALSE)   # Leer las clases e IDs de las muestras
colnames(clases) <- c("ID", "clase")
df <- cbind(clases, expresion)     #  Combinar todo en un solo DataFrame
rownames(df) <- df$ID # Establecer el ID de la muestra como nombre de fila

# Se evalua la existencia de NAS
anyNA(df)
na_counts <- colSums((is.na(df)))
na_counts

# Los datos no presentan valores perdidos, por lo que no es necesario imputar datos

# Se comprueba si hay muestras con alto contenido en valores 0
zero_counts <- colSums(df == 0)
zero_counts

# Los valores de cero indican la no detección de expresión para un gen en una muestra,
# por lo que se cuenta con ellos para realizar el análisis

# Graficamos la presencia de valores cero
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

# Aunque hay muestras con alto valor de valores cero, quizás sean genes que se 
# expresen en ciertos individuos, y contribuyen a la variabilidad biológica de los datos,
# por lo que mantenemos todas las muestras

# Comprobamos qué genes tienen varianza de cero (es decir, tienen un valor constante a lo largo de todas las muestras)
# Se calcula la varianza de cada gen (excluyendo las dos primeras columnas)
varianzas <- apply(df[, -c(1, 2)], 2, var)

# Identificacion de cuantos genes tienen varianza cero
var_genes <- names(varianzas[varianzas == 0])
cat("Se han encontrado", length(var_genes), "genes con varianza cero.")
print(head(var_genes, 10))

# Se eliminan los genes con varianza 0 y se genera un nuevo dataframe sin esos genes
# Esto es necesario para poder realizar metodos de reducción de dimensionalidad, 
# ya que los genes sin variación no nos aportan información.
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
##          1. MÉTODO REDUCCIÓN DIMENSIONALIDAD t-SNE          ##
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
##          2. MÉTODO REDUCCIÓN DIMENSIONALIDAD UMAP           ##
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

umap.results <- umap(df_final, n_neighbors=0.2 * nrow(df_final),
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
##          3. MÉTODO CLUSTERING NO JERÁRQUICO K-MEANS         ##
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


# k=5
kmeans.result <- kmeans(df_final_scale, centers = 5, iter.max = 100, nstart = 25)

# Visualizacion
fviz_cluster(kmeans.result, df_final_scale, 
             xlab = '', ylab = '',
             geom = "point") + # <--- Esta línea quita las etiquetas
  ggtitle("Cluster plot, centers = 5", subtitle = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(b = -10)))

#################################################################
##          4. MÉTODO CLUSTERING JERARQUICO                    ##
#################################################################

# Realizamos el análisis cluster para agrupar las muestras en racimos (clusters): 
# primero calculando la matriz de distancias euclídeas con dist(method = "euclidean"),
# después se genera el dendrograma con el método de Ward's D con hclust(method = "ward.D")
hclust_model <- df_final_scale %>%
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

# Preparamos los datos para ser analizados con métodos de aprendizaje automatizado supervisado

# Establecemos una semilla para poder replicar los resultados
set.seed(1234)

# Cargamos la librería 'caret' para usar funciones de partición de datos de entrenamiento y test, y acceder a los modelos de ML
library(caret)

# Cargamos la librería 'pROC' para calcular la curva ROC
library(pROC)

# Ceneramos una matriz de expresión con una columna con la clase de tumor
df_model <- as.data.frame(df_final[, -c(1,2)]) %>%
  mutate(Clase = df_final$clase)

# Convertimos la columna clase en un factor
df_model$Clase <- as.factor(df_model$Clase)

# Particionamos los datos en sets de Training (70 % de los individuos) y Testing (30 % de individuos)
idx <- createDataPartition(df_model$Clase, p = 0.7, list = FALSE)

train <- df_model[idx, ]
test  <- df_model[-idx, ]

# Usaremos los datos train y test para entrenar modelos de aprendizaje supervisado

#################################################################
##                     1. RANDOM FOREST                        ##
#################################################################

# Generamos un modelo de Random Forest para predecir la variable clase
rf_fit <- train(Clase ~ ., 
                data = train,
                method = "rf",
                trControl = trainControl(method = "cv", number = 10),
                preProcess = c("center", "scale"),
                tuneLength = 10,
                importance = TRUE
                )

rf_fit

plot(rf_fit) # La precisión para el RF fue óptima para el valor de mtry = 57, con 0.8445

# Calculamos las predicciones para los datos test
rf_predictions <- predict(rf_fit,
                          newdata = test
                          )

table(rf_predictions) # clases predichas por el RF
table(test$Clase)     # clases reales para datos test

# Comparamos nuestras predicciones con la variable Diagnosis mediante una 
# matriz de confusión:
confusionMatrix(
  rf_predictions,
  test$Clase
)

# Calculamos la probabilidad predicha para cada individuo de pertenecer a una clase de tumor
rf_probabilities <- predict(
  rf_fit,
  newdata = test,
  type = "prob"
)

rf_probabilities # clases predichas, pero en formato de probabilidad

# Generamos la curva ROC, esta vez usamos multiclass.roc() para poder comparar entre más de dos clases 
rf_roc <- multiclass.roc(
  response  = test$Clase,
  predictor = as.numeric(test$Clase)
)

# Representamos gráficamente las curvas ROC cuando comparamos con AGH
rocRF_AGH <- roc(test$Clase == "AGH", rf_probabilities[, "AGH"])
plot(rocRF_AGH, 
     col="red", 
     lwd=2, 
     legacy.axes=TRUE,
     main = paste("ROC del modelo RF"),
     )

rocRF_CFB <- roc(test$Clase == "CFB", rf_probabilities[, "CFB"])
plot.roc(
  rocRF_CFB,
  col  = "blue",
  add = TRUE,
  lwd = 2
  )

rocRF_CGC <- roc(test$Clase == "CGC", rf_probabilities[, "CGC"])
plot.roc(
  rocRF_CGC,
  col  = "green",
  add = TRUE,
  lwd = 2
  )

rocRF_CHC <- roc(test$Clase == "CHC", rf_probabilities[, "CHC"])
plot.roc(
  rocRF_CHC,
  col  = "orange",
  add = TRUE,
  lwd = 2
  )

rocRF_HPB <- roc(test$Clase == "HPB", rf_probabilities[, "HPB"])
plot.roc(
  rocRF_HPB,
  col  = "purple",
  add = TRUE,
  lwd = 2
)

# El ajuste de las curvas ROC sale perfecto para todas las comparaciones. 
# Descartamos el sobreajuste, dado que la precisión es alta y es distinta del azar cuando 
# la comparamos con la tasa de no información. 

# Calculamos el AUC
rf_auc <- auc(rf_roc)
cat("AUC Random Forest (raw):", rf_auc, "\n")

# El área debajo de la curva (AUC) sale también perfecto, como se esperaba tas ver las curvas ROC.

#################################################################
##                     2. SVM MODEL LINEAL                     ##
#################################################################

# Generamos un modelo de SVM Lineal para predecir la variable clase

svmModelLineal <- train(Clase ~.,
                        data = train,
                        method = "svmLinear",
                        trControl = trainControl(method = "cv", number = 10),
                        preProcess = c("center", "scale"),
                        tuneGrid = expand.grid(C = seq(0, 2, length = 20)),
                        prob.model = TRUE) 
svmModelLineal

plot(svmModelLineal) # La precisión óptima se consiguió con un valor de C = 0.11, con 0.8633

# Calculo de la precisión del modelo en el conjunto de prueba utilizando el modelo entrenado
# Se usa el modelo entrenado para predecir las etiquetas del conjunto de prueba

svm_predictions <- predict(svmModelLineal, newdata = test )
svm_predictions

# Evaluacion de la precisión del modelo utilizando la matriz de confusión

confusionMatrix(svm_predictions, test$Clase)

# SVM lineal. Calculo de las probabilidades del modelo de que sea la clase correcta

svm_probabilities <- predict(svmModelLineal, newdata = test, type = "prob")
svm_probabilities

### Curva ROC de SVM Lineal y AUC
# Clase AGH
rocSVM_AGH <- roc(test$Clase == "AGH", svm_probabilities[, "AGH"])
plot(rocSVM_AGH, col="red", lwd=2, legacy.axes=TRUE,
     main = paste("Curva ROC: Clase AGH\nAUC =", round(auc(rocSVM_AGH), 4)),
     xlab = "Especificidad",
     ylab = "Sensibilidad")

# Clase CFB
rocSVM_CFB <- roc(test$Clase == "CFB", svm_probabilities[, "CFB"])
plot(rocSVM_CFB, col="blue", lwd=2, legacy.axes=TRUE,
     main = paste("Curva ROC: Clase CFB\nAUC =", round(auc(rocSVM_CFB), 4)),
     xlab = "Especificidad",
     ylab = "Sensibilidad")

# Clase CGC
rocSVM_CGC <- roc(test$Clase == "CGC", svm_probabilities[, "CGC"])
plot(rocSVM_CGC, col="green", lwd=2, legacy.axes=TRUE,
     main = paste("Curva ROC: Clase CGC\nAUC =", round(auc(rocSVM_CGC), 4)),
     xlab = "Especificidad",
     ylab = "Sensibilidad")

# Clase CHC
rocSVM_CHC <- roc(test$Clase == "CHC", svm_probabilities[, "CHC"])
plot(rocSVM_CHC, col="orange", lwd=2, legacy.axes=TRUE,
     main = paste("Curva ROC: Clase CHC\nAUC =", round(auc(rocSVM_CHC), 4)),
     xlab = "Especificidad",
     ylab = "Sensibilidad")

# Clase HPB
rocSVM_HPB <- roc(test$Clase == "HPB", svm_probabilities[, "HPB"])
plot(rocSVM_HPB, col="purple", lwd=2, legacy.axes=TRUE,
     main = paste("Curva ROC: Clase HPB\nAUC =", round(auc(rocSVM_HPB), 4)),
     xlab = "Especificidad",
     ylab = "Sensibilidad")

###################################################################
##                     3. k-NEAREST NEIGHBORS                    ##
###################################################################

model_knn <- train(
  Clase ~ .,
  data=train,
  method="knn",
  trControl=trainControl(method = "cv", number=10),
  preProcess=c("center", "scale"),
  tuneLength=30
)

model_knn

# Se aplica el modelo entrenado con los datos nuevos para obtener las predicciones para evaluar cómo se comporta el modelo
predicciones <- predict(model_knn, newdata=test)

# Con las predicciones previas se calcula la matriz de confusión, que compara las predicciones con valores reales
confusionMatrix(predicciones, test$Clase)

# Se representa la precisión del modelo en base al número de vecinos
plot(model_knn) # La precisión óptima se consiguió con un valor de k = 11, con 0.9964

# Se representa la curva ROC para evaluar el modelo planteado
predicciones_prob <- predict(model_knn, newdata=test, type="prob")

roc_AGH <- plot(roc(test$Clase == "AGH", predicciones_prob[,1]),
                col="red",
                lwd=2,
                main = "ROC del modelo kNN"
                )
roc_CFB <- plot(roc(test$Clase == "CFB", predicciones_prob[,2]),
                col="blue",
                add=TRUE
                )
roc_CGC <- plot(roc(test$Clase == "CGC", predicciones_prob[,3]),
                col="green",
                add=TRUE
)
roc_CHC <- plot(roc(test$Clase == "CHC", predicciones_prob[,4]),
                col="orange",
                add=TRUE
)
roc_HPB <- plot(roc(test$Clase == "HPB", predicciones_prob[,5]),
                col="purple",
                add=TRUE
)

legend("bottomright",
       legend=c("AGH", "CFB", "CGC", "CHC", "HPB"),
       col=c("red", "blue", "green", "orange", "purple"),
       lwd=2)
