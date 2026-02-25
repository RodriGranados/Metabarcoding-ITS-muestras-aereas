###############################################################
# Cargar datos
###############################################################
library(readxl)
library(dplyr)

conteo_generos <- read_excel("Tabla conteo de generos y esporas.xlsx")





# Mostrar primeras filas para revisar estructura
head(conteo_generos)
View(conteo_generos)

## Eliminar las esporas si se requiere
conteo_generos_filtrado <- conteo_generos %>%
  filter(!grepl("Esporas", Genero))

conteo_generos <- conteo_generos_filtrado

###############################################################
# Librerías necesarias
###############################################################

library(ggplot2)
library(dplyr)


###############################################################
# Limpieza y preparación de datos
###############################################################

# Asegurar que Semana sea un factor (categoría), no número
conteo_generos$Semana <- factor(conteo_generos$Semana)

# Convertir Total.con.correccion a número (si viene con coma decimal)
conteo_generos$`Total con correccion` <- as.numeric(
  gsub(",", ".", conteo_generos$`Total con correccion`)
)

###############################################################
# 1) Gráfico horizontal: Total por género según semana
###############################################################

ggplot(conteo_generos, aes(x = reorder(Genero, Total),  # Ordenar por total
                           y = Total,
                           fill = Semana)) +
  geom_col() +                  # Barras
  coord_flip() +                # Hacer el gráfico horizontal
  labs(title = "Conteo de géneros y esporas",
       x = "Géneros y esporas",
       y = "Total",
       fill = "Semana") +
  theme_bw() +
  theme(text = element_text(size = 14))


###############################################################
# 2) Gráfico por temporada (facetas)
###############################################################

ggplot(conteo_generos, aes(x = reorder(Genero, Total),
                           y = Total,
                           fill = Semana)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Temporada, scales = "free_y") +  # Un panel por temporada
  labs(title = "Conteos según temporada",
       x = "Géneros y esporas",
       y = "Total") +
  theme_bw() +
  theme(text = element_text(size = 14))


###############################################################
# 3) Gráfico usando el valor corregido
###############################################################

ggplot(conteo_generos, aes(x = reorder(Genero, `Total con correccion`),
                           y = `Total con correccion`,
                           fill = Semana)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Temporada, scales = "free_y") +
  labs(title = "Conteo corregido según temporada",
       x = "Géneros y esporas",
       y = "Total (corregido)") +
  theme_bw()


###############################################################
# 4) Heatmap usando Z-score por semana dentro de cada temporada
###############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# Convertir semana a numérica para eje x del heatmap
conteo_generos$Semana <- as.numeric(conteo_generos$Semana)

# Calcular Z-score dentro de cada temporada para cada género
conteo_z <- conteo_generos %>%
  group_by(Temporada, Genero) %>%
  mutate(Zscore = scale(`Total con correccion`)) %>%
  ungroup()

# Heatmap
ggplot(conteo_z, aes(x = Semana, y = Genero, fill = Zscore)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, name = "Z-score"
  ) +
  facet_wrap(~ Temporada, scales = "free_y") +
  labs(title = "Heatmap por temporada usando Z-score",
       x = "Semana", y = "Género y esporas") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 0),
        panel.grid = element_blank())


###############################################################
# 5) Heatmap comparando temporadas (sin semanas)
###############################################################

# Sumar total corregido por temporada y género
conteo_temp <- conteo_generos %>%
  group_by(Temporada, Genero) %>%
  summarise(Total_corr = sum(`Total con correccion`), .groups = "drop")

# Calcular Z-score por género (comparando estaciones)
conteo_temp_z <- conteo_temp %>%
  group_by(Genero) %>%
  mutate(Zscore = scale(Total_corr)) %>%
  ungroup()

# Heatmap único
ggplot(conteo_temp_z, aes(x = Temporada, y = Genero, fill = Zscore)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "Z-score") +
  labs(title = "Abundancia relativa por temporada (Z-score)",
       x = "Temporada", y = "Género y esporas") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())


###############################################################
# 6) Barras por temporada
###############################################################

comp_temporada <- conteo_generos %>%
  group_by(Temporada, Genero) %>%
  summarise(Total = sum(`Total con correccion`), .groups = "drop")

ggplot(comp_temporada, aes(x = Temporada, y = Total, fill = Genero)) +
  geom_bar(stat = "identity") +
  labs(title = "Composición de géneros y esporas por temporada",
       x = "Temporada", y = "Total corregido") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


###############################################################
# 7) Gráfico de áreas apiladas por semana
###############################################################

ggplot(conteo_generos, aes(x = Semana, y = `Total con correccion`, fill = Genero)) +
  geom_area(alpha = 0.8) +
  facet_wrap(~ Temporada, scales = "free_x") +
  labs(title = "Dinámica semanal de géneros y esporas (área apilada)",
       x = "Semana", y = "Total corregido") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 0))


###############################################################
# 8) Boxplots por género comparando temporadas
###############################################################

ggplot(conteo_generos, aes(x = Temporada, y = `Total con correccion`, fill = Temporada)) +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_wrap(~ Genero, scales = "free_y") +
  labs(title = "Variación semanal por temporada (boxplots por género y esporas)",
       x = "Temporada", y = "Total corregido") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


###############################################################
# 9) Proporciones por temporada
###############################################################

prop_data <- conteo_generos %>%
  group_by(Temporada, Genero) %>%
  summarise(Total = sum(`Total con correccion`), .groups = "drop") %>%
  group_by(Temporada) %>%
  mutate(Prop = Total / sum(Total))

ggplot(prop_data, aes(x = Temporada, y = Prop, fill = Genero)) +
  geom_bar(stat = "identity") +
  labs(title = "Proporción de géneros y esporas por temporada",
       y = "Proporción") +
  theme_minimal()


###############################################################
# 10) Sankey plot (flujo Temporada → Género)
###############################################################

library(ggalluvial)

sankey_data <- conteo_generos %>%
  group_by(Temporada, Genero) %>%
  summarise(Total = sum(`Total con correccion`), .groups = "drop")

ggplot(sankey_data,
       aes(axis1 = Temporada, axis2 = Genero, y = Total)) +
  geom_alluvium(aes(fill = Genero), alpha = 0.7) +
  geom_stratum() +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  labs(title = "Sankey - Flujo entre Temporada y Género")


###############################################################
# 11) Chord Diagram
###############################################################

library(circlize)
library(tibble)

mat <- sankey_data %>%
  tidyr::pivot_wider(names_from = Genero, values_from = Total,
                     values_fill = 0) %>%
  column_to_rownames("Temporada") %>%
  as.matrix()

chordDiagram(mat)


###############################################################
# 12) Barras finales por temporada
###############################################################

ggplot(conteo_generos,
       aes(x = Genero, y = `Total con correccion`, fill = Genero)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Temporada, scales = "free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Abundancia por género y esporas en las 4 estaciones")


