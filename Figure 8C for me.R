library(ggplot2)


# Leer el archivo Excel (ajusta la ruta al archivo)
data <- read.csv("Data_ABI5-ABA.tsv", sep = "\t", header = TRUE)

# Verifica la estructura de tus datos para asegurarte de que estén correctamente cargados
head(data)

ggplot(data, aes(x = reorder(Symbol, log2FoldChange), y = log2FoldChange, color = log2FoldChange > 0)) +
  geom_segment(aes(xend = Symbol, y = 0, yend = log2FoldChange), size = 1.2) +  # Líneas más gruesas
  geom_point(size = 4, shape = 21, fill = "white") +  # Puntos más grandes y huecos
  coord_flip() +
  scale_color_manual(values = c("steelblue", "firebrick")) +  # Otros colores
  theme_minimal() +
  labs(x = "Genes", y = "log2FoldChange", title = "Lollipop Plot de Expresión Génica") +
  theme(axis.text.y = element_text(size = 10, face = "bold"),  # Ajustar tamaño de los nombres de genes
        plot.title = element_text(hjust = 0.5))  # Centrar el título



# Crear el Lollipop Plot ajustado
ggplot(data, aes(x = reorder(Symbol, log2FoldChange), y = log2FoldChange)) +
  
# Segmento que conecta el eje con el punto
geom_segment(aes(xend = Symbol, color = log2FoldChange > 0, y = 0, yend = log2FoldChange), size = 2) +
  
# Los puntos, con jitter para evitar superposición, color según FoldChange y tamaño según adjPvalue
geom_point(aes(color = log2FoldChange > 0, size = -log10(padj)), 
           shape = 21, stroke = 1, fill = "white") +
  
  # Voltear el gráfico para tener genes en el eje Y
  coord_flip() +
  
  # Colores para valores positivos y negativos de FoldChange (sin mostrar la leyenda de color)
  scale_color_manual(values = c("lightblue", "lightcoral"), guide = "none") +
  
  # Ajustar el tamaño de los círculos (escala de tamaño)
  scale_size_continuous(range = c(2, 8)) +  # Ajustar el rango de tamaños de los puntos
  
  # Etiquetas y título del gráfico
  labs(x = "Genes", y = "L2FC", title = "Possible ABI5 targets- ABA responding genes",
       size = "-log10(padj)") +
  
  # Ajustar el estilo de los nombres de los genes (cursiva y tamaño más pequeño)
  theme_minimal() +
  
  theme(
    
    # Ajustar el tamaño de las leyendas de los ejes
    axis.title.x = element_text(size =20, face = "bold"),  # Tamaño de la leyenda del eje X
    axis.title.y = element_text(size = 20, face = "bold", vjust = 0.5, hjust = 0.5),  # Tamaño de la leyenda del eje Y
    
    # Ajustar el tamaño de los números de la escala
    axis.text.x = element_text(size = 18),  # Tamaño de los números del eje X
    axis.text.y = element_text(size = 18, face = "italic"),  # Tamaño de los nombres de los genes
    
    # Ajustar el tamaño de la leyenda de los puntos
    legend.title = element_text(size = 18, face = "bold"),  # Tamaño del título de la leyenda
    legend.text = element_text(size = 18),  # Tamaño del texto de la leyenda
    
    # Eliminar la cuadrícula
    panel.grid.minor = element_blank(),
    
    # Añadir borde negro alrededor de la figura
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Borde negro

    
    # Centrar el título
    plot.title = element_text(hjust = 0.5, size= 18, face= "bold"),
    
    # Ajustar los márgenes: (top, right, bottom, left)
    plot.margin = margin(8, 8,20, 80))
    
  
