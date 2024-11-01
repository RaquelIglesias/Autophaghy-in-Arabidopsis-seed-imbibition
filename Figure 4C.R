

library(ggplot2)


data <- read.csv("GenesSel_GO.tsv", sep = "\t", header = TRUE)
head(data)# Asegúrate de que 'data' contenga la columna de categorías llamada 'categories'
# Asegúrate de que 'data' contenga la columna de categorías llamada 'categories'



data$categories <- factor(data$categories, 
                          levels = c("Lipid storage", 
                                     "Seed maturation", 
                                     "Selected", 
                                     "Response to ER-stress"))

# Crear el Lollipop Plot ajustado
ggplot(data, aes(x = reorder(symbol, log2FoldChange), y = log2FoldChange, fill = categories)) + 
  # Segmento que conecta el eje con el punto
  geom_segment(aes(xend = symbol, y = 0, yend = log2FoldChange, color = log2FoldChange > 0), size = 2) +
  
  # Los puntos, con tamaño según adjPvalue y color según log2FoldChange
  geom_point(aes(size = -log10(padj), color = log2FoldChange > 0), shape = 21, stroke = 1, fill = "white") + 
  
  # Voltear el gráfico para tener genes en el eje Y
  coord_flip() + 
  
  # Colores para valores positivos y negativos de FoldChange
  scale_color_manual(values = c("lightblue", "lightcoral"), guide = "none") + 
  
  # Ajustar el tamaño de los círculos (escala de tamaño)
  scale_size_continuous(range = c(2, 6)) + 
  
  # Ajustar el estilo de los nombres de los genes (cursiva y tamaño más pequeño)
  theme_minimal() + 
  
  theme(
    # Ajustar el tamaño de las leyendas del eje Y
    axis.title.y = element_blank(),  # Eliminar el título del eje Y
    
    # Ajustar el tamaño de los números de la escala
    axis.text.x = element_text(size = 18),  # Tamaño de los números del eje X
    axis.text.y = element_text(size = 7.5, face = "bold.italic"),  # Tamaño de los nombres de los genes más reducido
    
    # Ajustar el tamaño de la leyenda de las categorías
    legend.title = element_text(size = 20, face = "bold"),  # Tamaño del título de la leyenda
    legend.text = element_text(size = 18),  # Tamaño del texto de la leyenda
    
    # Eliminar la cuadrícula
    panel.grid.minor = element_blank(),
    
    # Añadir borde negro alrededor de la figura
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    
    # Centrar el título
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    
    # Ajustar los márgenes: más espacio alrededor del gráfico
    plot.margin = margin(20, 20, 20, 100)
  ) + 
  
  # Ajustar la posición de los genes en el eje Y según su categoría
  facet_grid(categories ~ ., scales = "free_y", space = "free_y") +  # Crea facetas por categoría
  
  # Ajustar el espacio adicional en el eje X (para el log2FoldChange)
  scale_y_continuous(limits = c(-3, 7), expand = expansion(mult = c(0.1, 0.1))) +  # Aumentar espacio extra arriba y abajo
  
  # Ajustar el tamaño de la leyenda de -log10(padj)
  guides(size = guide_legend(title = "-log10(padj)", title.theme = element_text(size = 16), label.theme = element_text(size = 12))) +
  
  # Eliminar título del eje X
  labs(x = NULL, y = NULL)  # Sin título para eje X e Y

