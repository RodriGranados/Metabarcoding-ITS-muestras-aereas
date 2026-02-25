
###############################################################
# 7️⃣ Cargar archivo de salida de FUNGuild
#    (el archivo .guilds.txt generado por FUNGuild en el paso anterior)
###############################################################
funguild_table <- read.delim("OTU_abundance_taxonomy_funguild.guilds.txt",
                             header = TRUE, sep = "\t")

head(funguild_table)


###############################################################
# 8️⃣ Limpiar espacios extra en la columna "Trophic.Mode"
###############################################################
library(stringr)
funguild_table$Trophic.Mode <- str_trim(funguild_table$Trophic.Mode,
                                        side = "both")




library(dplyr)

funguild_summary <- funguild_table %>%
  mutate(
    funguild_status = ifelse(Guild == "-" | is.na(Guild),
                             "No identificado",
                             "Identificado")
  )


conteo_otus <- funguild_summary %>%
  count(funguild_status)

conteo_otus

otus_no_identificados <- funguild_summary %>%
  filter(funguild_status == "No identificado") %>%
  select(OTU_ID, taxonomy, Guild)

otus_no_identificados



library(dplyr)
library(stringr)

no_id_breakdown <- funguild_summary %>%
  filter(funguild_status == "No identificado") %>%
  mutate(
    no_id_tipo = case_when(
      is.na(taxonomy) ~ "NA",
      str_detect(taxonomy, "Incertae") ~ "Incertae sedis",
      TRUE ~ "Otros no identificados"
    )
  ) %>%
  count(no_id_tipo)

no_id_breakdown


otros_no_identificados <- funguild_summary %>%
  filter(
    funguild_status == "No identificado",
    !is.na(taxonomy),
    !str_detect(taxonomy, "Incertae")
  ) %>%
  select(OTU_ID, taxonomy, Guild, Trophic.Mode)

otros_no_identificados


otros_no_identificados %>%
  distinct(taxonomy) %>%
  nrow()


lista_taxonomy_otros <- otros_no_identificados %>%
  distinct(taxonomy) %>%
  arrange(taxonomy)

lista_taxonomy_otros

###############################################################
# 9️⃣ Convertir abundancias a formato largo (long format)
#    - Esto permite hacer sumas por Guild, gráficos, etc.
###############################################################
library(tidyverse)

fungi_long <- funguild_table %>%
  pivot_longer(cols = PG1:PG7_8,
               names_to = "Sample",
               values_to = "Abundance")

###############################################################
# 🔟 Identificar los 10 Guilds más abundantes (Top 10)
#    - Se suman abundancias por Guild a través de todas las muestras.
###############################################################
top10_guilds <- fungi_long %>%
  group_by(Guild) %>%
  summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total)) %>%
  slice(1:10) %>%
  pull(Guild)

###############################################################
# 1️⃣1️⃣ Crear categoría "Other" para todos los Guilds que NO están
#      en el Top 10.
###############################################################
fungi_long_top <- fungi_long %>%
  mutate(Guild2 = ifelse(Guild %in% top10_guilds, Guild, "Other"))

# Ordenar niveles: Top10 primero, luego Other
fungi_long_top$Guild2 <- factor(
  fungi_long_top$Guild2,
  levels = c(top10_guilds, "Other")
)

###############################################################
# 1️⃣2️⃣ Gráfico de barras de abundancia por Guild (Top10 + Other)
###############################################################
ggplot(fungi_long_top, aes(x = Sample, y = Abundance, fill = Guild2)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Abundancia por Guild (Top 10 + Other)",
       y = "Abundancia",
       x = "",
       fill = "Guild") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###############################################################
# 1️⃣3️⃣ AGRUPAR POR Trophic.Mode (Saprotroph, Pathotroph, Symbiotroph...)
#      - Suma abundancias por muestra.
###############################################################
trophic <- fungi_long %>%
  group_by(Sample, `Trophic.Mode`) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Gráfico por modo trófico
ggplot(trophic, aes(x = Sample, y = Abundance, fill = `Trophic.Mode`)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Trophic Modes por muestra", y = "Abundancia")








###############################################################
# 1️⃣4️⃣ Calcular abundancia total por Guild (en toda la tabla)
#      - rowSums: abundancia total por OTU
#      - luego sumar por Guild
###############################################################
guild_total <- funguild_table %>%
  mutate(Total = rowSums(across(PG1:PG7_8))) %>%
  group_by(Guild) %>%
  summarise(Abundance = sum(Total), .groups = "drop")


###############################################################
# 1️⃣5️⃣ Selección del Top10 Guilds totales
###############################################################
top10 <- guild_total %>%
  arrange(desc(Abundance)) %>%
  slice(1:10) %>%
  pull(Guild)


###############################################################
# 1️⃣6️⃣ Reasignar Guilds fuera del Top10 → "Other"
###############################################################
guild_total_top <- guild_total %>%
  mutate(Guild2 = ifelse(Guild %in% top10, Guild, "Other"))

# Sumar por categoría agrupada
guild_total_top_sum <- guild_total_top %>%
  group_by(Guild2) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Ordenar niveles para graficar
guild_total_top_sum$Guild2 <- factor(
  guild_total_top_sum$Guild2,
  levels = c(top10, "Other")
)


###############################################################
# 1️⃣7️⃣ Gráfico circular (pie chart) de abundancia total por Guild
###############################################################
ggplot(guild_total_top_sum, aes(x = "", y = Abundance, fill = Guild2)) +
  geom_col() +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Distribución total de Guilds (Top 10 + Other)")
###############################################################
# 🔍 A) TOP taxones por Trophic.Mode
###############################################################

top_taxa_trophic <- funguild_table %>%
  mutate(Total = rowSums(across(PG1:PG7_8))) %>%   # Abundancia total por taxón
  group_by(Trophic.Mode, taxonomy) %>%             # Agrupar por modo trófico
  summarise(Abundance = sum(Total), .groups = "drop") %>%
  arrange(Trophic.Mode, desc(Abundance))

# Mostrar los 10 más abundantes por Trophic.Mode
top_taxa_trophic %>% 
  group_by(Trophic.Mode) %>%
  slice_max(order_by = Abundance, n = 10)


library(ggplot2)

top_taxa_trophic %>%
  group_by(Trophic.Mode) %>%
  slice_max(order_by = Abundance, n = 10) %>%
  ggplot(aes(x = reorder(taxonomy, Abundance), y = Abundance, fill = Trophic.Mode)) +
  geom_col() +
  facet_wrap(~Trophic.Mode, scales="free_y") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 taxones por Trophic.Mode",
       x = "Taxón",
       y = "Abundancia total")





###############################################################
# 🔍 B) TOP taxones por Guild
###############################################################

top_taxa_guild <- funguild_table %>%
  mutate(Total = rowSums(across(PG1:PG7_8))) %>%   # Abundancia total
  group_by(Guild, taxonomy) %>%                    # Agrupar por guild
  summarise(Abundance = sum(Total), .groups = "drop") %>%
  arrange(Guild, desc(Abundance))

# Mostrar los 10 más abundantes por Guild
top_taxa_guild %>%
  group_by(Guild) %>%
  slice_max(order_by = Abundance, n = 10)

top_taxa_guild %>%
  group_by(Guild) %>%
  slice_max(order_by = Abundance, n = 10) %>%
  ggplot(aes(x = reorder(taxonomy, Abundance), y = Abundance, fill = Guild)) +
  geom_col() +
  facet_wrap(~Guild, scales="free_y") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 taxones por Guild",
       x = "Taxón",
       y = "Abundancia total")


###############################################################
# TOP 10 GUILDS MÁS ABUNDANTES
###############################################################

top10_guilds <- funguild_table %>%
  mutate(Total = rowSums(across(PG1:PG7_8))) %>%   # Sumar abundancia total por OTU
  group_by(Guild) %>%                              # Agrupar por Guild
  summarise(Abundance = sum(Total), .groups = "drop") %>% 
  arrange(desc(Abundance)) %>%                     # Ordenar descendente
  slice(1:10)                                      # Tomar Top 10

top10_guilds

funguild_top10 <- funguild_table %>%
  mutate(Total = rowSums(across(PG1:PG7_8))) %>%
  filter(Guild %in% top10_guilds$Guild)


###############################################################
# TOP ESPECIES DENTRO DE CADA GUILD DEL TOP 10
###############################################################

top_species_per_guild <- funguild_top10 %>%
  group_by(Guild, taxonomy) %>% 
  summarise(Abundance = sum(Total), .groups = "drop") %>%
  arrange(Guild, desc(Abundance))

top10_species_per_guild <- top_species_per_guild %>%
  group_by(Guild) %>%
  slice_max(order_by = Abundance, n = 10)


ggplot(top10_species_per_guild,
       aes(x = reorder(taxonomy, Abundance), y = Abundance, fill = Guild)) +
  geom_col() +
  facet_wrap(~ Guild, scales = "free_y") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 especies dentro de los Guilds más abundantes",
       x = "Especie",
       y = "Abundancia")


### NUEVO INTENTO CON FUNGUILD ###3

library(tidyverse)

abund_tbl <- feature_table %>%
  rownames_to_column(var = "OTU_ID")


tax_tbl <- tax_split_2 %>%
  select(Feature.ID, Taxon_7) %>%
  rename(
    OTU_ID = Feature.ID,
    taxonomy = Taxon_7
  )

final_tbl <- abund_tbl %>%
  left_join(tax_tbl, by = "OTU_ID") %>%
  select(
    OTU_ID,
    PG1, PG2, PG3, PG4, PG5_6, PG7_8,
    taxonomy
  )


write.table(
  final_tbl,
  file = "OTU_abundance_taxonomy_funguild.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)







library(dplyr)
library(stringr)

tax_split_2 <- tax_split_1 %>%
  mutate(
    Taxon_1 = ifelse(is.na(Taxon_1), NA,
                     paste0("k__", str_replace_all(Taxon_1, " ", "_"))),
    
    Taxon_2 = ifelse(is.na(Taxon_2), NA,
                     paste0("p__", str_replace_all(Taxon_2, " ", "_"))),
    
    Taxon_3 = ifelse(is.na(Taxon_3), NA,
                     paste0("c__", str_replace_all(Taxon_3, " ", "_"))),
    
    Taxon_4 = ifelse(is.na(Taxon_4), NA,
                     paste0("o__", str_replace_all(Taxon_4, " ", "_"))),
    
    Taxon_5 = ifelse(is.na(Taxon_5), NA,
                     paste0("f__", str_replace_all(Taxon_5, " ", "_"))),
    
    Taxon_6 = ifelse(is.na(Taxon_6), NA,
                     paste0("g__", str_replace_all(Taxon_6, " ", "_"))),
    
    Taxon_7 = ifelse(is.na(Taxon_7), NA,
                     paste0("s__", str_replace_all(Taxon_7, " ", "_")))
  )




### GRAFICOS DE COMPOSICION MODO TROFICO Y GUILD ###

metadata <- tibble(
  Muestra   = c("PG1","PG2","PG3","PG4","PG5_6","PG7_8"),
  Mes       = c("Junio","Julio","Septiembre","Octubre","Febrero","Abril"),
  Temporada = c("Invierno","Invierno","Primavera","Primavera","Verano","Otoño")
)


plot_funguild_level <- function(
    funguild_table,
    metadata,
    rank = "Guild",
    x_var = c("Mes", "Temporada"),
    top = 10,
    threshold = 0.03
) {
  
  x_var <- match.arg(x_var)
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(RColorBrewer)
  
  ## 1️⃣ Formato largo + metadata
  df_long <- funguild_table %>%
    pivot_longer(
      cols = PG1:PG7_8,
      names_to = "Muestra",
      values_to = "Abundance"
    ) %>%
    left_join(metadata, by = "Muestra") %>%
    mutate(
      rank_val = case_when(
        is.na(.data[[rank]]) ~ "No asignados",
        .data[[rank]] == "-" ~ "No asignados",
        TRUE ~ as.character(.data[[rank]])
      )
    )
  
  ## 2️⃣ Abundancia relativa por muestra
  df_rel <- df_long %>%
    group_by(.data[[x_var]], rank_val) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    group_by(.data[[x_var]]) %>%
    mutate(Abundance = Abundance / sum(Abundance)) %>%
    ungroup()
  
  ## 3️⃣ Top N (promedio entre muestras)
  top_levels <- df_rel %>%
    group_by(rank_val) %>%
    summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = top) %>%
    pull(rank_val)
  
  ## 4️⃣ Agrupar el resto como "Other"
  df_rel <- df_rel %>%
    mutate(
      rank2 = ifelse(rank_val %in% top_levels, rank_val, "Other")
    ) %>%
    group_by(.data[[x_var]], rank2) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  ## 5️⃣ Reordenar niveles (Other al final)
  orden_rank <- df_rel %>%
    group_by(rank2) %>%
    summarise(total = sum(Abundance)) %>%
    arrange(desc(total)) %>%
    pull(rank2)
  
  df_rel$rank2 <- factor(
    df_rel$rank2,
    levels = c(setdiff(orden_rank, "Other"), "Other")
  )
  
  ## 6️⃣ Etiquetas (> threshold)
  df_rel <- df_rel %>%
    mutate(
      label = ifelse(
        Abundance >= threshold,
        paste0(round(Abundance * 100, 1), "%"),
        NA
      )
    )
  
  ## 7️⃣ Colores
  n_levels <- length(levels(df_rel$rank2))
  n_sin_other <- n_levels - 1
  
  paleta <- if (n_sin_other <= 12) {
    brewer.pal(n_sin_other, "Paired")
  } else {
    colorRampPalette(brewer.pal(12, "Paired"))(n_sin_other)
  }
  
  colores <- c(paleta, "grey70")
  
  ## 8️⃣ Gráfico
  ggplot(df_rel, aes(x = .data[[x_var]], y = Abundance, fill = rank2)) +
    geom_bar(stat = "identity", color = "black") +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 3,
      na.rm = TRUE
    ) +
    scale_fill_manual(values = colores) +
    labs(
      title = paste("Composición funcional por", rank, "según", x_var),
      x = x_var,
      y = "Abundancia relativa",
      fill = rank
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5
      ),
      legend.position = "right"
    ) +
    coord_cartesian(ylim = c(0, 1))
}



plot_funguild_level(
  funguild_table,
  metadata,
  rank = "Guild",
  x_var = "Mes",
  top = 15,
  threshold = 0.03
)



plot_funguild_level(
  funguild_table,
  metadata,
  rank = "Trophic.Mode",
  x_var = "Mes",
  top = 10,
  threshold = 0.05
)


plot_funguild_level(
  funguild_table,
  metadata,
  rank = "Growth.Morphology",
  x_var = "Mes",
  top = 10,
  threshold = 0.05
)

### ESPECIES QUE GOBIERNAN CADA MODO TROFICO Y GUILD

library(dplyr)
library(tidyr)

fungi_long <- funguild_table %>%
  pivot_longer(
    cols = PG1:PG7_8,
    names_to = "Muestra",
    values_to = "Abundance"
  ) %>%
  mutate(
    taxonomy = ifelse(is.na(taxonomy), "No asignados", taxonomy),
    Guild = ifelse(is.na(Guild) | Guild == "-", "No asignados", Guild),
    Trophic.Mode = ifelse(is.na(Trophic.Mode) | Trophic.Mode == "-", "No asignados", Trophic.Mode)
  )

top_species_trophic <- fungi_long %>%
  group_by(Trophic.Mode, taxonomy) %>%
  summarise(Total_abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Trophic.Mode) %>%
  arrange(desc(Total_abundance)) %>%
  slice_head(n = 5)

top_species_trophic


top_species_guild <- fungi_long %>%
  group_by(Guild, taxonomy) %>%
  summarise(Total_abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Guild) %>%
  arrange(desc(Total_abundance)) %>%
  slice_head(n = 5)

top_species_guild





library(dplyr)
library(tidyr)
library(ggplot2)

# Datos largos + limpieza
fungi_long <- funguild_table %>%
  pivot_longer(
    cols = PG1:PG7_8,
    names_to = "Muestra",
    values_to = "Abundance"
  ) %>%
  mutate(
    taxonomy = ifelse(is.na(taxonomy), "No asignados", taxonomy),
    Guild = ifelse(is.na(Guild) | Guild == "-", "No asignados", Guild),
    Trophic.Mode = ifelse(is.na(Trophic.Mode) | Trophic.Mode == "-", "No asignados", Trophic.Mode)
  )

# Top 5 especies por Trophic.Mode (abundancia relativa)
top_trophic <- fungi_long %>%
  group_by(Trophic.Mode, taxonomy) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  group_by(Trophic.Mode) %>%
  mutate(Rel = Total / sum(Total)) %>%
  arrange(desc(Rel)) %>%
  slice_head(n = 5)

# Gráfico
ggplot(top_trophic,
       aes(x = reorder(taxonomy, Rel), y = Rel, fill = taxonomy)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~ Trophic.Mode, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Especies dominantes por modo trófico",
    x = NULL,
    y = "Abundancia relativa"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )



## intento de grafico Nro. 3 - EL MEJOR

limpiar_especie <- function(x) {
  x <- gsub("^s__", "", x)                # quitar s__
  x <- gsub("_sp\\.", " sp.", x)           # _sp. → sp.
  x <- gsub("_", " ", x)                   # _ → espacio
  x <- gsub("Incertae sedis", "Incertae_sedis", x)
  x
}


top_trophic <- fungi_long %>%
  group_by(Trophic.Mode, taxonomy) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  group_by(Trophic.Mode) %>%
  mutate(Rel = Total / sum(Total)) %>%
  arrange(Trophic.Mode, desc(Rel)) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    especie = limpiar_especie(taxonomy),
    x_label = especie
  )

top_trophic$x_label <- factor(
  top_trophic$x_label,
  levels = top_trophic$x_label
)

separadores <- top_trophic %>%
  count(Trophic.Mode) %>%
  mutate(pos = cumsum(n)) %>%
  slice(-n())


ggplot(top_trophic,
       aes(x = x_label, y = Rel, fill = Trophic.Mode)) +
  
  geom_col(color = "black", linewidth = 0.2) +
  
  geom_vline(
    xintercept = separadores$pos + 0.5,
    linetype = "dashed",
    color = "grey40"
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  labs(
    title = "Especies dominantes por modo trófico",
    x = NULL,
    y = "Abundancia relativa",
    fill = "Modo trófico"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      face = "italic"        # 👈 CURSIVA
    ),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )


#####

top15_guilds <- fungi_long %>%
  group_by(Guild) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = 15) %>%
  pull(Guild)

top_guild <- fungi_long %>%
  filter(Guild %in% top15_guilds) %>%          # 👈 FILTRO CLAVE
  group_by(Guild, taxonomy) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  group_by(Guild) %>%
  mutate(Rel = Total / sum(Total)) %>%
  arrange(Guild, desc(Rel)) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    especie = limpiar_especie(taxonomy),
    x_label = especie
  )


top_guild$x_label <- factor(
  top_guild$x_label,
  levels = top_guild$x_label
)

separadores_guild <- top_guild %>%
  count(Guild) %>%
  mutate(pos = cumsum(n)) %>%
  slice(-n())


ggplot(top_guild,
       aes(x = x_label, y = Rel, fill = Guild)) +
  
  geom_col(color = "black", linewidth = 0.2) +
  
  geom_vline(
    xintercept = separadores_guild$pos + 0.5,
    linetype = "dashed",
    color = "grey40"
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  labs(
    title = "Especies dominantes en los 15 guilds más abundantes",
    x = NULL,
    y = "Abundancia relativa",
    fill = "Guild"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      face = "italic"
    ),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

######

top10_morph <- fungi_long %>%
  group_by(Growth.Morphology) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = 10) %>%
  pull(Growth.Morphology)

top_morph <- fungi_long %>%
  filter(Growth.Morphology %in% top10_morph) %>%   # 👈 filtro clave
  group_by(Growth.Morphology, taxonomy) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  group_by(Growth.Morphology) %>%
  mutate(Rel = Total / sum(Total)) %>%
  arrange(Growth.Morphology, desc(Rel)) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    especie = limpiar_especie(taxonomy),
    x_label = especie
  )

top_morph$x_label <- factor(
  top_morph$x_label,
  levels = top_morph$x_label
)

separadores_morph <- top_morph %>%
  count(Growth.Morphology) %>%
  mutate(pos = cumsum(n)) %>%
  slice(-n())


ggplot(top_morph,
       aes(x = x_label, y = Rel, fill = Growth.Morphology)) +
  
  geom_col(color = "black", linewidth = 0.2) +
  
  geom_vline(
    xintercept = separadores_morph$pos + 0.5,
    linetype = "dashed",
    color = "grey40"
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  labs(
    title = "Especies dominantes en las morfologías de crecimiento más abundantes",
    x = NULL,
    y = "Abundancia relativa",
    fill = "Morfología de crecimiento"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      face = "italic"
    ),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

#####

library(scales)

heat_df <- fungi_long %>%
  group_by(Guild, taxonomy) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  group_by(Guild) %>%
  mutate(Rel = Total / sum(Total)) %>%
  ungroup()

ggplot(heat_df,
       aes(x = Guild, y = taxonomy, fill = Rel)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Distribución funcional de especies fúngicas",
    x = "Guild",
    y = "Especie",
    fill = "Abundancia relativa"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )






library(tidyverse)
library(pheatmap)

fungi_long <- funguild_table %>%
  pivot_longer(
    cols = PG1:PG7_8,
    names_to = "Muestra",
    values_to = "Abundancia"
  ) %>%
  mutate(
    Especie = ifelse(is.na(taxonomy), "No asignados", taxonomy),
    trophic.mode = ifelse(
      is.na(Trophic.Mode) | Trophic.Mode == "-",
      "No asignados",
      Trophic.Mode
    )
  )

abund_tm_especie <- fungi_long %>%
  group_by(trophic.mode, Especie) %>%
  summarise(Abundancia = sum(Abundancia), .groups = "drop")

top_15_tm <- abund_tm_especie %>%
  group_by(trophic.mode) %>%
  summarise(Total_TM = sum(Abundancia)) %>%
  arrange(desc(Total_TM)) %>%
  slice_head(n = 15)

abund_tm_especie_filt <- abund_tm_especie %>%
  filter(trophic.mode %in% top_15_tm$trophic.mode)


top_5_especies_tm <- abund_tm_especie_filt %>%
  group_by(trophic.mode) %>%
  arrange(desc(Abundancia)) %>%
  slice_head(n = 5) %>%
  ungroup()

heatmap_data <- top_5_especies_tm %>%
  pivot_wider(
    names_from = Especie,
    values_from = Abundancia,
    values_fill = 0
  ) %>%
  column_to_rownames("trophic.mode")

heatmap_data_log <- log10(heatmap_data + 1)

pheatmap(
  heatmap_data_log,
  scale = "row",
  clustering_method = "complete",
  fontsize_row = 10,
  fontsize_col = 8,
  angle_col = 45,
  main = "Top 15 trophic.mode y Top 5 especies más abundantes"
)






abund_tm_especie <- data %>%
  group_by(trophic.mode, Especie) %>%
  summarise(Abundancia = sum(Abundancia), .groups = "drop")

top_tm <- abund_tm_especie %>%
  group_by(trophic.mode) %>%
  summarise(Total = sum(Abundancia)) %>%
  slice_max(Total, n = 15) %>%
  pull(trophic.mode)

heatmap_df <- abund_tm_especie %>%
  filter(trophic.mode %in% top_tm) %>%
  group_by(trophic.mode) %>%
  slice_max(Abundancia, n = 5) %>%
  ungroup()

heatmap_df <- heatmap_df %>%
  group_by(trophic.mode) %>%
  mutate(Abund_z = as.numeric(scale(Abundancia))) %>%
  ungroup()

library(ggplot2)

ggplot(heatmap_df,
       aes(x = Especie,
           y = trophic.mode,
           fill = Abund_z)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Z-score"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title = element_blank()
  )




library(dplyr)

tm_temporada <- funguild_table %>%
  group_by(trophic.mode, Temporada) %>%
  summarise(Abundancia = sum(Abundancia), .groups = "drop")


