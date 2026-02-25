#######################################################################
# SCRIPT PARA GENERAR TABLAS DE PRESENCIA/ABSENCIA POR MES, TEMPORADA,
# Y COMBINACIONES ÚNICAS DE ESPECIES
#
# Incluye:
# 1) Tabla PA por mes
# 2) Listas de especies por mes (Excel)
# 3) Tabla PA por temporada
# 4) Listas de especies por temporada (Excel)
# 5) Identificación de combinaciones únicas de meses donde aparece cada especie
#######################################################################

library(tidyverse)
library(writexl)
library(phyloseq)

#######################################################################
# 1️⃣ Crear tabla PA (Presencia/Ausencia) por especie y mes
#######################################################################
# - tax_glom: colapsa el objeto phyloseq a nivel de "Especie".
# - psmelt: convierte el objeto phyloseq a formato largo (tabla típica).
# - group_by: agrupa por Mes y Especie sumando abundancias.
# - mutate(Presencia): convierte abundancia > 0 en presencia (=1).
# - pivot_wider: convierte la tabla a formato de matriz con meses como columnas.
#######################################################################

species_pa_mes <- tax_glom(ps, taxrank = "Especie") %>%
  psmelt() %>%
  group_by(Meses, Especie) %>%
  summarise(Abundancia = sum(Abundance), .groups = "drop") %>%
  mutate(Presencia = ifelse(Abundancia > 0, 1, 0)) %>%
  pivot_wider(
    names_from = Meses,
    values_from = Presencia,
    values_fill = 0
  )

#######################################################################
# 2️⃣ Combinar especies duplicadas
#######################################################################
# - Algunas especies pueden aparecer en múltiples filas por cómo psmelt
#   maneja datos taxonómicos.
# - Se agrupan filas duplicadas por Especie.
# - across(any()): si al menos una fila tenía 1 en un mes, la especie se mantiene presente.
#######################################################################

species_pa_mes <- species_pa_mes %>%
  group_by(Especie) %>%
  summarise(across(where(is.numeric), ~ as.numeric(any(.x > 0))), .groups = "drop")

# Eliminar columna Abundancia si quedó por accidente
species_pa_mes <- species_pa_mes %>% select(-Abundancia)

# Verificar estructura final
head(species_pa_mes)

#######################################################################
# 3️⃣ Convertir tabla PA por mes a formato largo y generar listas por mes
#######################################################################
# - Se pasa a formato largo para identificar especies presentes en cada mes.
# - Se filtra solo las especies presentes (Presencia == 1).
# - group_by(Mes) permite crear una lista por cada mes.
# - Luego se rellena con NA para igualar longitudes y poder exportar a Excel.
#######################################################################

species_long <- species_pa_mes %>%
  pivot_longer(
    cols = -Especie,
    names_to = "Mes",
    values_to = "Presencia"
  ) %>%
  filter(Presencia == 1) %>%   # solo especies presentes
  select(-Presencia)

species_list <- species_long %>%
  group_by(Mes) %>%
  summarise(Especies = list(Especie)) %>%
  deframe()

# Igualar largo de listas para exportación
max_len <- max(sapply(species_list, length))
species_equal <- lapply(species_list, function(x) c(x, rep(NA, max_len - length(x))))
species_df <- as.data.frame(species_equal)

# Exportar a Excel
write_xlsx(species_df, "RESULTADOS/Taxones compartidos/especies_por_mes.xlsx")

#######################################################################
# 4️⃣ Generar tabla PA por Temporada
#######################################################################
# - Se definen las temporadas basadas en los meses donde hay muestreo.
# - Verano: Febrero
# - Otoño: Abril
# - Invierno: Junio y Julio
# - Primavera: Septiembre y Octubre
# - Si aparece en cualquier mes de una temporada, se asigna 1.
#######################################################################

species_pa_temp <- species_pa_mes %>%
  mutate(
    Verano = ifelse(Febrero == 1, 1, 0),
    Otoño = ifelse(Abril == 1, 1, 0),
    Invierno = ifelse(Junio == 1 | Julio == 1, 1, 0),
    Primavera = ifelse(Septiembre == 1 | Octubre == 1, 1, 0)
  ) %>%
  select(Especie, Verano, Otoño, Invierno, Primavera)


#######################################################################
# 5️⃣ Crear listas por temporada y exportar a Excel
#######################################################################
# - Cada temporada se convierte en un vector de especies presentes.
# - Se rellenan las columnas con NA para alinearlas.
# - write_xlsx genera el archivo final.
#######################################################################

verano <- species_pa_temp %>% filter(Verano == 1) %>% pull(Especie)
otoño <- species_pa_temp %>% filter(Otoño == 1) %>% pull(Especie)
invierno <- species_pa_temp %>% filter(Invierno == 1) %>% pull(Especie)
primavera <- species_pa_temp %>% filter(Primavera == 1) %>% pull(Especie)

max_len <- max(length(verano), length(otoño), length(invierno), length(primavera))
species_by_temp <- data.frame(
  Verano = c(verano, rep(NA, max_len - length(verano))),
  Otoño = c(otoño, rep(NA, max_len - length(otoño))),
  Invierno = c(invierno, rep(NA, max_len - length(invierno))),
  Primavera = c(primavera, rep(NA, max_len - length(primavera)))
)


write_xlsx(species_by_temp, "RESULTADOS/Taxones compartidos/especies_por_temporada.xlsx")

#######################################################################
# 6️⃣ Generar combinaciones únicas de meses por especie
#######################################################################
# - Para cada especie se identifica exactamente en qué meses aparece.
# - Esto permite encontrar patrones de ocurrencia únicos.
# - Se agrupan las especies que comparten la misma combinación.
# - El resultado se ordena por número de especies por combinación.
# - Se exporta en formato Excel.
#######################################################################

library(openxlsx)

species_comb <- species_pa_mes %>%
  mutate(
    Combinacion = apply(select(., -Especie), 1, function(x) {
      meses_presentes <- names(x)[which(x == 1)]
      if (length(meses_presentes) == 0) {
        return(NA)
      } else {
        return(paste(sort(meses_presentes), collapse = "_"))
      }
    })
  )

species_unique_combos <- species_comb %>%
  filter(!is.na(Combinacion)) %>%
  group_by(Combinacion) %>%
  summarise(
    N_especies = n(),
    Especies = paste(Especie, collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(N_especies))


write.xlsx(
  species_unique_combos,
  file = "RESULTADOS/Taxones compartidos/Especies_combinaciones_unicas_meses.xlsx",
  asTable = TRUE
)




#######################################################################
# COMBINACIONES ÚNICAS DE TEMPORADAS POR ESPECIE
# Todas las especies en columnas (sin límite)
#######################################################################

library(tidyverse)
library(openxlsx)

#######################################################################
# 1️⃣ Crear combinación de temporadas por especie
#######################################################################

species_comb_temp <- species_pa_temp %>%
  mutate(
    Combinacion = apply(select(., -Especie), 1, function(x) {
      temporadas_presentes <- names(x)[x == 1]
      if (length(temporadas_presentes) == 0) {
        NA
      } else {
        paste(sort(temporadas_presentes), collapse = "_")
      }
    })
  ) %>%
  filter(!is.na(Combinacion))

#######################################################################
# 2️⃣ Agrupar especies por combinación
#######################################################################

species_by_comb <- species_comb_temp %>%
  group_by(Combinacion) %>%
  summarise(
    N_especies = n(),
    Especies = list(sort(Especie)),
    .groups = "drop"
  ) %>%
  arrange(desc(N_especies))

#######################################################################
# 3️⃣ Expandir lista de especies a columnas (dinámico)
#######################################################################

max_len <- max(lengths(species_by_comb$Especies))

species_comb_expanded <- species_by_comb %>%
  mutate(
    Especies = map(
      Especies,
      ~ c(.x, rep(NA, max_len - length(.x)))
    )
  ) %>%
  unnest_wider(Especies, names_sep = "_")

#######################################################################
# 4️⃣ Renombrar columnas de especies
#######################################################################

colnames(species_comb_expanded) <- c(
  "Combinacion",
  "N_especies",
  paste0("Especie_", seq_len(max_len))
)

#######################################################################
# 5️⃣ Exportar a Excel
#######################################################################

write.xlsx(
  species_comb_expanded,
  file = "RESULTADOS/Taxones compartidos/Especies_combinaciones_unicas_temporadas.xlsx",
  asTable = TRUE
)
