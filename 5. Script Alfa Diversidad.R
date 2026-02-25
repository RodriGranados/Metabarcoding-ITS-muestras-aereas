### ============================================================
###     ANÁLISIS DE DIVERSIDAD ALFA
### ============================================================

library(vegan)
library(tidyverse)

### FILTRO
library(dplyr)


# === 5. Rarefacción ===

colSums(feature_table)
min_depth <- min(colSums(feature_table))
cat("🔹 Profundidad mínima para rarefacción:", min_depth, "\n")


# Aplicar rarefacción a igual profundidad
otu_rarefied <- rrarefy(t(feature_table), sample = min_depth)
otu_rarefied <- t(otu_rarefied)

colSums(otu_rarefied)
# === 6. Calcular métricas de diversidad alfa ===

## Observed (riqueza observada)
observed <- specnumber(t(otu_rarefied))

## Chao1 y ACE
estimadores <- estimateR(t(otu_rarefied))
chao1 <- estimadores["S.chao1", ]
ace <- estimadores["S.ACE", ]

## Shannon, Simpson, Inverse Simpson, Fisher
shannon <- diversity(t(otu_rarefied), index = "shannon")
simpson <- diversity(t(otu_rarefied), index = "simpson")
invsimpson <- 1 / simpson
fisher <- fisher.alpha(t(otu_rarefied))



num_individuos <- colSums(otu_rarefied)
num_individuos_original <- colSums(feature_table)

# === 7. Combinar métricas en una tabla final ===
alpha_metrics <- data.frame(
  Sample = names(observed),
  Observed = observed,
  Chao1 = chao1,
  ACE = ace,
  Shannon = shannon,
  Simpson = simpson,
  InvSimpson = invsimpson,
  Fisher = fisher
)

# === 8. Mostrar y guardar resultados ===
print(alpha_metrics)

library(openxlsx)
write.xlsx(alpha_metrics, "Indices_alfa_diversidad_molecular.xlsx")

# === 7. Combinar métricas en una tabla final ===
alpha_metrics <- data.frame(
  Sample = names(observed),
  Individuos = num_individuos,   # <-- NUEVA COLUMNA
  Observed = observed,
  Chao1 = chao1,
  ACE = ace,
  Shannon = shannon,
  Simpson = simpson,
  InvSimpson = invsimpson,
  Fisher = fisher
)




alpha_metrics <- data.frame(
  Sample = names(observed),
  Individuos_Original = num_individuos_original[names(observed)],  # <-- NUEVA COLUMNA
  Individuos_Rarefactados = num_individuos,                       # <-- YA EXISTENTE
  Observed = observed,
  Chao1 = chao1,
  ACE = ace,
  Shannon = shannon,
  Simpson = simpson,
  InvSimpson = invsimpson,
  Fisher = fisher
)

print(alpha_metrics)

## EVALUACION DE NORMALIDAD
indices <- c("Shannon", "Simpson", "InvSimpson",
             "Chao1", "ACE", "Observed", "Fisher")

normalidad <- lapply(indices, function(i) {
  shapiro.test(alpha_metrics[[i]])
})

names(normalidad) <- indices
normalidad


### EVALUACION DE HOMOGENEIDAD
library(car)

homogeneidad <- lapply(indices, function(i) {
  leveneTest(alpha_metrics[[i]] ~ alpha_metrics$Sample)
})

names(homogeneidad) <- indices
homogeneidad


## Guardar metricas de alfa diversidad
write.table(alpha_metrics, "Resultados/Metricas Alfa Diversidad1.csv",
            sep = "\t", quote = FALSE, row.names = FALSE)






# --- 9. Asignar temporadas a las muestras ---
alpha_metrics <- alpha_metrics %>%
  mutate(Season = case_when(
    Sample %in% c("PG1", "PG2") ~ "Invierno",
    Sample %in% c("PG3", "PG4") ~ "Primavera",
    Sample == "PG5_6" ~ "Verano",
    Sample == "PG7_8" ~ "Otoño",
    TRUE ~ "Desconocido"
  ))

alpha_metrics$Season <- factor(alpha_metrics$Season,
                               levels = c("Invierno", "Primavera", "Verano", "Otoño"))


# --- 10. Preparar datos para gráficos ---
alpha_long <- alpha_metrics %>%
  pivot_longer(cols = c(Observed, Chao1, ACE, Shannon, Simpson, InvSimpson, Fisher),
               names_to = "Index", values_to = "Value")

# --- 11. Crear gráfico boxplot facetado ---
p <- ggplot(alpha_long, aes(x = Season, y = Value, fill = Season)) +
  geom_boxplot(alpha = 0.75, width = 0.7, color = "black") +
  facet_wrap(~ Index, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(
    "Invierno" = "#5DADE2",   # azul
    "Primavera" = "#58D68D",  # verde
    "Verano" = "#F4D03F",     # amarillo
    "Otoño" = "#CA6F1E"       # marrón
  )) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  ) +
  labs(title = "Índices de Diversidad Alfa por Temporada (ITS)",
       x = "Temporada", y = "Valor del índice")

print(p)


# asegurarse que Season sea factor
alpha_metrics$Season <- factor(alpha_metrics$Season)

# lista de índices
indices <- c("Shannon", "Simpson", "InvSimpson",
             "Observed", "Chao1", "ACE", "Fisher")

# Kruskal–Wallis para todos los índices
kw_results <- lapply(indices, function(idx) {
  kruskal.test(as.formula(paste(idx, "~ Season")), data = alpha_metrics)
})

names(kw_results) <- indices
kw_results




welch_results <- lapply(indices, function(idx) {
  oneway.test(
    as.formula(paste(idx, "~ Season")),
    data = alpha_metrics,
    var.equal = FALSE
  )
})

names(welch_results) <- indices
welch_results

library(rstatix)

games_results <- lapply(indices, function(idx) {
  games_howell_test(
    alpha_metrics,
    formula = as.formula(paste(idx, "~ Season"))
  )
})

names(games_results) <- indices
games_results


### ALFA DIVERSIDAD CON MICROECO ####

library(microeco)
metadata_df      <- as.data.frame(metadata)
rownames(metadata_df) <- metadata_df$Muestra
metadata_df$Muestra <- NULL

feature_table_df <- as.data.frame(feature_table)

tax_split_df     <- as.data.frame(tax_split)
rownames(tax_split_df) <- tax_split_df$Feature.ID
tax_split_df$Feature.ID <- NULL

mt <- microtable$new(
  sample_table = metadata_df,
  otu_table    = feature_table_df,
  tax_table    = tax_split_df
)

mt

mt$cal_abund()
mt$cal_alphadiv()
mt$cal_betadiv()


colSums(mt$otu_table)

tn <- trans_norm$new(dataset = mt)
mt_rarefy <- tn$norm(method="rarefy", sample.size=87821)



mt_rarefy$cal_alphadiv()
mt_rarefy$alpha_diversity

alpha_df <- mt_rarefy$alpha_diversity

library(dplyr)

alpha_df <- alpha_df %>%
  rownames_to_column(var = "Muestra") %>%
  left_join(
    metadata_df %>% rownames_to_column("Muestra"),
    by = "Muestra"
  )



t1 <- trans_alpha$new(dataset = mt_rarefy, group = "Temporada")

t1$cal_diff(method = "KW")
t1$res_diff

t1$cal_diff(method = "wilcox")
t1$res_diff


t1$cal_diff(method = "anova")
t1$res_diff


t1 <- trans_alpha$new(dataset = mt_rarefy, group = "Temporada")

t1$cal_diff(method="anova")

t1$res_diff


t1$plot_alpha(measure = "Simpson")

t1$plot_alpha(measure = "Fisher")

t1$plot_alpha(measure = "Chao1")

### BETA DIVERSIDAD CON MICROECO ####

mt_rarefy$cal_betadiv()
mt_rarefy$beta_diversity$bray %>% View()

t2 <- trans_beta$new(dataset = mt_rarefy, group = "Temporada", measure = "bray")

t2$cal_ordination(method="PCoA")

t2$plot_ordination()

t2$plot_ordination(plot_color = "Temporada", plot_type = c("point","ellipse"))

p <- t2$plot_ordination(
  plot_color = "Mes",
  plot_shape = "Temporada",
  plot_type = c("point", "ellipse"),
  point_size = 4
)
p +
  scale_shape_manual(values = c(
    "Verano"     = 17,  # triangulo lleno
    "Otoño"      = 5,  # rombo vacio
    "Invierno"   = 16,  # circulo lleno
    "Primavera"  = 15   # cuadrado lleno
  )) +
  labs(shape = "Temporada")




t2$cal_group_distance(within_group = F)
t2$cal_group_distance_diff(method="anova")
t2$res_group_distance_diff

t2$plot_group_distance()




t2$plot_clustering()
t2$plot_clustering(replace_name = c("Temporada"), group = "Temporada")



t2$cal_ordination(
  method = "NMDS")

p1 <- t2$plot_ordination(
  plot_color = "Mes",
  plot_shape = "Temporada",
  plot_type = c("point", "ellipse"),
  point_size = 4
)
p1 +
  scale_shape_manual(values = c(
    "Verano"     = 17,  # triangulo lleno
    "Otoño"      = 5,  # rombo vacio
    "Invierno"   = 16,  # circulo lleno
    "Primavera"  = 15   # cuadrado lleno
  )) +
  labs(shape = "Temporada")




t2$cal_betadiv_stat(
  method = "permanova",
  group = "Temporada"
)

