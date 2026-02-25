### CARGA DE LIBRERIAS
library(tidyverse)
library(UpSetR)

### ORGANIZACION DE LA METADATA
metadata <- tibble(
  Muestra   = c("PG1","PG2","PG3","PG4","PG5_6","PG7_8"),
  Mes       = c("Junio","Julio","Septiembre","Octubre","Febrero","Abril"),
  Temporada = c("Invierno","Invierno","Primavera","Primavera","Verano","Otoño")
)


### FUNCION GENERAL PARA OBTENER PA (PRESENCIA/AUSENCIA) DESDE
### EL ARCHIVO DE ABUNDANCIA RELATIVA
get_pa_from_rel_abund <- function(file, rank, metadata, group_var) {
  
  # Leer archivo
  df <- read.delim(file, header = TRUE, sep = "\t", check.names = FALSE)
  
  colnames(df)[1] <- rank
  
  # Convertir abundancias en presencia/ausencia
  pa <- df %>%
    mutate(across(-all_of(rank), ~ ifelse(.x > 0, 1, 0)))
  
  # Reestructurar y unirse a la metadata
  pa2 <- pa %>%
    pivot_longer(-all_of(rank), names_to = "Muestra", values_to = "Presencia") %>%
    left_join(metadata, by = "Muestra") %>%
    group_by(!!sym(rank), !!sym(group_var)) %>%
    summarise(Presencia = as.numeric(any(Presencia > 0)), .groups = "drop") %>%
    pivot_wider(names_from = !!sym(group_var), values_from = Presencia, values_fill = 0)
  
  # Asignar nombres de fila
  rownames(pa2) <- pa2[[rank]]
  pa2 <- pa2 %>% select(-all_of(rank))
  
  return(as.data.frame(pa2))
}


### GENERAR TABLAS DE PA POR MES
phylum_pa_mes  <- get_pa_from_rel_abund("Resultados/Abundancias separadas por nivel/abundancia_phylum.tsv", "Phylum", metadata, "Mes")
clase_pa_mes  <- get_pa_from_rel_abund("Resultados/Abundancias separadas por nivel/abundancia_clase.tsv", "Clase", metadata, "Mes")
familia_pa_mes  <- get_pa_from_rel_abund("Resultados/Abundancias separadas por nivel/abundancia_familia.tsv", "Familia", metadata, "Mes")
especie_pa_mes   <- get_pa_from_rel_abund("Resultados/Abundancias separadas por nivel/abundancia_especie.tsv",  "Especie",  metadata, "Mes")
genero_pa_mes   <- get_pa_from_rel_abund("Resultados/Abundancias separadas por nivel/abundancia_genero.tsv",  "Genero",  metadata, "Mes")
orden_pa_mes   <- get_pa_from_rel_abund("Resultados/Abundancias separadas por nivel/abundancia_orden.tsv",  "Orden",  metadata, "Mes")


### GENERAR TABLAS DE PA POR TEMPORADA
phylum_pa_temp  <- get_pa_from_rel_abund("Resultados/Abundancias separadas por nivel/abundancia_phylum.tsv", "Phylum", metadata, "Temporada")
clase_pa_temp  <- get_pa_from_rel_abund("Resultados/Abundancias separadas por nivel/abundancia_clase.tsv", "Clase", metadata, "Temporada")
familia_pa_temp <- get_pa_from_rel_abund("Resultados/Abundancias separadas por nivel/abundancia_familia.tsv", "Familia", metadata, "Temporada")
genero_pa_temp  <- get_pa_from_rel_abund("Resultados/Abundancias separadas por nivel/abundancia_genero.tsv",  "Genero",  metadata, "Temporada")
especie_pa_temp <- get_pa_from_rel_abund("Resultados/Abundancias separadas por nivel/abundancia_especie.tsv", "Especie", metadata, "Temporada")
orden_pa_temp <- get_pa_from_rel_abund("Resultados/Abundancias separadas por nivel/abundancia_orden.tsv", "Orden", metadata, "Temporada")


### GENERACION DE GRAFICOS UPSET POR MES

upset(phylum_pa_mes,
      sets = colnames(phylum_pa_mes),
      main.bar.color = "#1b9e77",
      sets.bar.color = "#7570b3",
      mainbar.y.label = "Phylum compartidos",
      sets.x.label = "Phylum por mes",
      order.by = "freq"
)

upset(especie_pa_mes,
      sets = colnames(especie_pa_mes),
      main.bar.color = "#1b9e77",
      sets.bar.color = "#7570b3",
      mainbar.y.label = "Especies compartidos",
      sets.x.label = "Especies por mes",
      order.by = "freq"
)


upset(genero_pa_mes,
      sets = colnames(genero_pa_mes),
      main.bar.color = "#1b9e77",
      sets.bar.color = "#7570b3",
      mainbar.y.label = "Genero compartidos",
      sets.x.label = "Genero por mes",
      order.by = "freq"
)

upset(orden_pa_mes,
      sets = colnames(orden_pa_mes),
      main.bar.color = "#1b9e77",
      sets.bar.color = "#7570b3",
      mainbar.y.label = "Orden compartidos",
      sets.x.label = "Orden por mes",
      order.by = "freq"
)


upset(clase_pa_mes,
      sets = colnames(clase_pa_mes),
      main.bar.color = "#1b9e77",
      sets.bar.color = "#7570b3",
      mainbar.y.label = "Clase compartidos",
      sets.x.label = "Clase por mes",
      order.by = "freq"
)


upset(familia_pa_mes,
      sets = colnames(familia_pa_mes),
      main.bar.color = "#1b9e77",
      sets.bar.color = "#7570b3",
      mainbar.y.label = "Familia compartidos",
      sets.x.label = "Familia por mes",
      order.by = "freq"
)


##### PARA LOS DIAGRAMAS DE VENN EN TEMPORADAS #######


### FUNCION PARA GENERAR LISTAS PARA VENN
get_sets_for_venn <- function(pa_table) {
  lapply(colnames(pa_table), function(col) {
    rownames(pa_table)[pa_table[[col]] == 1]
  }) %>% setNames(colnames(pa_table))
}

### ORGANIZAR VENN POR TEMPORADA
familia_sets_temp  <- get_sets_for_venn(familia_pa_temp)
genero_sets_temp   <- get_sets_for_venn(genero_pa_temp)
especie_sets_temp <- get_sets_for_venn(especie_pa_temp)
clase_sets_temp <- get_sets_for_venn(clase_pa_temp)
orden_sets_temp <- get_sets_for_venn(orden_pa_temp)
phylum_sets_temp <- get_sets_for_venn(phylum_pa_temp)

### GRAFICOS DE VENN
library(ggvenn)

ggvenn(clase_sets_temp,
       fill_color = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"),
       stroke_size = 0.7) +
  ggtitle("Clases compartidas entre temporadas") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggvenn(familia_sets_temp,
       fill_color = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"),
       stroke_size = 0.7) +
  ggtitle("Familias compartidas entre temporadas") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggvenn(genero_sets_temp,
       fill_color = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"),
       stroke_size = 0.7) +
  ggtitle("Géneros compartidos entre temporadas") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggvenn(especie_sets_temp,
       fill_color = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"),
       stroke_size = 0.7) +
  ggtitle("Especies compartidas entre temporadas") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggvenn(orden_sets_temp,
       fill_color = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"),
       stroke_size = 0.7) +
  ggtitle("Orden compartidos entre temporadas") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggvenn(phylum_sets_temp,
       fill_color = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"),
       stroke_size = 0.7) +
  ggtitle("Phylum compartidas entre temporadas") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))





#### ANALISIS DE ASVs COMPARTIDOS Y EXCLUSIVOS ###

asv_pa_mes <- get_pa_from_rel_abund(
  file = "feature-table.tsv",
  rank = "ASV",
  metadata = metadata,
  group_var = "Mes"
)

asv_pa_temp <- get_pa_from_rel_abund(
  file = "feature-table.tsv",
  rank = "ASV",
  metadata = metadata,
  group_var = "Temporada"
)


asv_sets_temp <- get_sets_for_venn(asv_pa_temp)
library(ggvenn)

ggvenn(asv_sets_temp,
       fill_color = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"),
       stroke_size = 0.7) +
  ggtitle("ASVs compartidos entre temporadas") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))



asv_exclusivos_mes <- asv_pa_mes %>%
  rownames_to_column("ASV") %>%
  mutate(Total = rowSums(across(-ASV))) %>%
  filter(Total == 1)

asv_exclusivos_mes


asv_exclusivos_mes_long <- asv_exclusivos_mes %>%
  pivot_longer(-c(ASV, Total), names_to = "Mes", values_to = "Presencia") %>%
  filter(Presencia == 1) %>%
  select(ASV, Mes)

asv_exclusivos_mes_long


asv_compartidos_todos_mes <- asv_pa_mes %>%
  rownames_to_column("ASV") %>%
  filter(rowSums(across(-ASV)) == ncol(asv_pa_mes))

asv_compartidos_todos_mes


asv_compartidos_parcial_mes <- asv_pa_mes %>%
  rownames_to_column("ASV") %>%
  mutate(Total = rowSums(across(-ASV))) %>%
  filter(Total >= 2 & Total < ncol(asv_pa_mes))

asv_compartidos_parcial_mes


upset(asv_pa_mes,
      sets = colnames(asv_pa_mes),
      main.bar.color = "#1b9e77",
      sets.bar.color = "#7570b3",
      mainbar.y.label = "ASVs compartidos",
      sets.x.label = "ASVs por mes",
      order.by = "freq"
)




library(tidyverse)

asv_long <- as.data.frame(feature_table) %>%
  rownames_to_column("ASV") %>%
  pivot_longer(
    cols = -ASV,
    names_to = "Muestra",
    values_to = "Abundance"
  ) %>%
  filter(Abundance > 0)


asv_long <- asv_long %>%
  left_join(metadata %>% select(Muestra, Mes, Temporada),
            by = "Muestra")



asvs_por_mes <- asv_long %>%
  distinct(ASV, Mes) %>%
  group_by(Mes) %>%
  summarise(
    ASVs_totales = n()
  )
asvs_por_temporada <- asv_long %>%
  distinct(ASV, Temporada) %>%
  group_by(Temporada) %>%
  summarise(
    ASVs_totales = n()
  )


asvs_exclusivos_mes <- asv_long %>%
  distinct(ASV, Mes) %>%
  group_by(ASV) %>%
  filter(n() == 1) %>%     # solo aparece en un mes
  ungroup() %>%
  group_by(Mes) %>%
  summarise(
    ASVs_exclusivos = n()
  )
asvs_exclusivos_temporada <- asv_long %>%
  distinct(ASV, Temporada) %>%
  group_by(ASV) %>%
  filter(n() == 1) %>%     # solo aparece en una temporada
  ungroup() %>%
  group_by(Temporada) %>%
  summarise(
    ASVs_exclusivos = n()
  )


tabla_mes <- asvs_por_mes %>%
  left_join(asvs_exclusivos_mes, by = "Mes") %>%
  replace_na(list(ASVs_exclusivos = 0))


tabla_temporada <- asvs_por_temporada %>%
  left_join(asvs_exclusivos_temporada, by = "Temporada") %>%
  replace_na(list(ASVs_exclusivos = 0))
