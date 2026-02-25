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


library(dplyr)
library(stringr)

conteo_generos <- conteo_generos %>%
  mutate(
    Genero = case_when(
      Genero == "Dictyoartrinium sp" ~ "Dictyoarthrinium sp.",
      Genero == "Lasiodiplidia sp" ~ "Lasiodiplodia sp.",
      Genero == "Tetraploide sp" ~ "Tetraploa sp.",
      Genero == "Fusarium spp" ~ "Fusarium sp.",
      Genero == "Hifa de Aspergillus" ~ "Aspergillus sp.",
      Genero == "Curvularia spp" ~ "Curvularia sp.",
      TRUE ~ Genero
    )
  )

conteo_generos <- conteo_generos %>%
  mutate(
    Genero = str_replace(Genero, "\\bsp$", "sp.")
  )

library(dplyr)

conteo_generos <- conteo_generos %>%
  rename(
    Hongos = Genero
  )

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

ggplot(conteo_generos, aes(x = reorder(Hongos, Total),  # Ordenar por total
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

ggplot(conteo_generos, aes(x = reorder(Hongos, Total),
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

ggplot(conteo_generos, aes(x = reorder(Hongos, `Total con correccion`),
                           y = `Total con correccion`,
                           fill = Semana)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Temporada, scales = "free_y") +
  labs(title = "Conteo corregido según temporada",
       x = "Hongos",
       y = "Total (corregido)") +
  theme_bw() +
  theme(
    axis.text.y = element_text(face = "italic")
  )


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
  group_by(Temporada, Hongos) %>%
  mutate(Zscore = scale(`Total con correccion`)) %>%
  ungroup()

# Heatmap
ggplot(conteo_z, aes(x = Semana, y = Hongos, fill = Zscore)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, name = "Z-score"
  ) +
  facet_wrap(~ Temporada, scales = "free_y") +
  labs(title = "Heatmap por temporada usando Z-score",
       x = "Semana", y = "Hongos") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 0),
        panel.grid = element_blank(), axis.text.y = element_text(face = "italic"))


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
  group_by(Temporada, Hongos) %>%
  summarise(Total = sum(`Total con correccion`), .groups = "drop")

ggplot(comp_temporada, aes(x = Temporada, y = Total, fill = Hongos)) +
  geom_bar(stat = "identity") +
  labs(title = "Composición de géneros y esporas por temporada",
       x = "Temporada", y = "Total corregido") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


###############################################################
# 7) Gráfico de áreas apiladas por semana
###############################################################

ggplot(conteo_generos, aes(x = Semana, y = `Total con correccion`, fill = Hongos)) +
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
  facet_wrap(~ Hongos, scales = "free_y") +
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



plot_generos_temporada_all <- function(data,
                                       group = "Temporada",
                                       tax = "Hongos",
                                       value = "Total con correccion",
                                       threshold = 0.03) {
  
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  
  # Orden de las temporadas
  orden_grupo <- unique(data[[group]])
  
  # Calcular proporciones
  df <- data %>%
    group_by(!!sym(group), !!sym(tax)) %>%
    summarise(Total = sum(.data[[value]]), .groups = "drop") %>%
    group_by(!!sym(group)) %>%
    mutate(Abundance = Total / sum(Total)) %>%
    ungroup()
  
  # Reordenar géneros por abundancia global
  taxa_order <- df %>%
    group_by(!!sym(tax)) %>%
    summarise(total = sum(Abundance)) %>%
    arrange(desc(total)) %>%
    pull(!!sym(tax))
  
  df[[tax]] <- factor(df[[tax]], levels = taxa_order)
  
  # Reordenar temporadas
  df[[group]] <- factor(df[[group]], levels = orden_grupo)
  
  # Etiquetas (%)
  df <- df %>%
    mutate(label = ifelse(Abundance >= threshold,
                          paste0(round(Abundance * 100, 1), "%"),
                          NA))
  
  # Paleta (15 colores → interpolación suave)
  n_taxa <- length(levels(df[[tax]]))
  colores <- if (n_taxa <= 12) {
    brewer.pal(n_taxa, "Paired")
  } else {
    colorRampPalette(brewer.pal(12, "Paired"))(n_taxa)
  }
  
  # Gráfico
  ggplot(df, aes(x = !!sym(group), y = Abundance, fill = !!sym(tax))) +
    geom_bar(stat = "identity", color = "black") +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              size = 3,
              color = "black",
              na.rm = TRUE) +
    scale_fill_manual(values = colores) +
    labs(
      title = "Composición de hongos por conteo morfologíco",
      x = NULL,
      y = "Abundancia relativa"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = element_text(hjust = 0.5),
      legend.position = "right", legend.text = element_text(face = "italic")
    ) +
    coord_cartesian(ylim = c(0, 1))
}

plot_generos_temporada_all(conteo_generos, threshold = 0.03)


###############################################################
# 10) Sankey plot (flujo Temporada → Género)
###############################################################

library(ggalluvial)

sankey_data <- conteo_generos %>%
  group_by(Temporada, Hongos) %>%
  summarise(Total = sum(`Total con correccion`), .groups = "drop")

ggplot(sankey_data,
       aes(axis1 = Temporada, axis2 = Hongos, y = Total)) +
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
  tidyr::pivot_wider(names_from = Hongos, values_from = Total,
                     values_fill = 0) %>%
  column_to_rownames("Temporada") %>%
  as.matrix()

chordDiagram(mat)



library(circlize)

circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = 4,          # más espacio entre sectores
  track.margin = c(0.01, 0.01)
)

chordDiagram(
  mat,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.15)
)

# Etiquetas manuales
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector_name = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.text(
      x = mean(xlim),
      y = ylim[1] + 0.2,
      labels = sector_name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.6,   # 🔹 aquí controlas el tamaño
      font = ifelse(grepl("sp|spp", sector_name), 3, 1)
    )
  },
  bg.border = NA
)



###############################################################
# 12) Barras finales por temporada
###############################################################

ggplot(conteo_generos,
       aes(x = Hongos, y = `Total con correccion`, fill = Hongos)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Temporada, scales = "free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Abundancia por género y esporas en las 4 estaciones")






library(dplyr)

conteo_hongo_semana <- conteo_generos %>%
  group_by(Temporada, Semana, Hongos) %>%
  summarise(
    Total_hongo = sum(`Total`, na.rm = TRUE),
    .groups = "drop"
  )


conteo_total_semana <- conteo_generos %>%
  group_by(Temporada, Semana) %>%
  summarise(
    Total_hongos = sum(`Total`, na.rm = TRUE),
    .groups = "drop"
  )


conteo_total_temporada <- conteo_generos %>%
  group_by(Temporada) %>%
  summarise(
    Total_hongos = sum(`Total`, na.rm = TRUE),
    .groups = "drop"
  )

