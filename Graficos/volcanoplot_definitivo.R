########### Código para generar el VolcanoPlot de los datos  ################

# Cargamos la libreria
library(EnhancedVolcano)

# Generamos data frame con el nombre de los genes
res_limpios_df <- as.data.frame(res_limpios)
res_limpios_df$gene <- rownames(res_limpios_df)

# Se genera el gráfico
EnhancedVolcano(res_limpios_df,
                lab = res_limpios_df$gene,
                title = "Expresión diferencial de genes",
                subtitle = "Volcano Plot",
                titleLabSize = 20,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'Log2foldchange (Obeso1 vs Normopeso)'),
                ylab = bquote(~-Log[10]~ 'p-value'),
                ylim = c(0, 2),     # ponemos 2 como límite del eje Y
                pCutoff = 0.05,     # perfilamos umbral de significancia <0,05
                labSize = 4.0,      # tamaño nombre genes
                axisLabSize = 10,
                colAlpha = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.5)



<img width="963" height="782" alt="volcanoplot" src="https://github.com/user-attachments/assets/c15e975c-7fa4-4c94-8f0f-1f9efe6d0ac3" />

