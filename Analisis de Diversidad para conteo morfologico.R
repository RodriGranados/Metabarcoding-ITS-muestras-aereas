library(dplyr)
library(tidyr)
library(vegan)
library(readxl)
library(fossil)
library(tibble)

# Cargar archivo
##df <- read_excel("conteo_generos.xlsx")##


df_tot <- conteo_generos %>%
  select(Temporada, Hongos, Total)

# Agrupar porque cada temporada tiene varias semanas
df_summarized_tot <- df_tot %>%
  group_by(Temporada, Hongos) %>%
  summarise(Total = sum(Total), .groups = "drop")

# Matriz de abundancia (filas = temporada, columnas = géneros)
mat_abund <- df_summarized_tot %>%
  pivot_wider(names_from = Hongos,
              values_from = Total,
              values_fill = 0) %>%
  column_to_rownames("Temporada")





mat_abund_int <- apply(mat_abund, 2, as.integer)
rownames(mat_abund_int) <- rownames(mat_abund)



min_n <- min(rowSums(mat_abund_int))
min_n



set.seed(123)  # reproducibilidad

mat_raref <- rrarefy(mat_abund_int, sample = min_n)

observed <- specnumber(mat_raref)
shannon  <- diversity(mat_raref, index = "shannon")
simpson  <- diversity(mat_raref, index = "simpson")
invsimp  <- diversity(mat_raref, index = "invsimpson")
chao1    <- estimateR(mat_raref)["S.chao1", ]
ace      <- apply(mat_raref, 1, fossil::ACE)
fisher   <- fisher.alpha(mat_raref)

alpha_raref <- data.frame(
  Sample = rownames(mat_raref),
  Observed = observed,
  Chao1 = chao1,
  ACE = ace,
  Shannon = shannon,
  Simpson = simpson,
  InvSimpson = invsimp,
  Fisher = fisher
)

alpha_raref

### ANALISIS ALFA DIVERSIDAD USANDO SEMANAS COMO REPLICA

library(dplyr)
library(tidyr)
library(vegan)
library(tibble)
library(fossil)

df_tot <- conteo_generos %>%
  select(Temporada, Semana, Hongos, Total)

df_week <- df_tot %>%
  group_by(Temporada, Semana, Hongos) %>%
  summarise(Total = sum(Total), .groups = "drop")

mat_abund_week <- df_week %>%
  unite("Sample", Temporada, Semana, sep = "_") %>%
  pivot_wider(
    names_from = Hongos,
    values_from = Total,
    values_fill = 0
  ) %>%
  column_to_rownames("Sample")

mat_abund_int <- apply(mat_abund_week, 2, as.integer)
rownames(mat_abund_int) <- rownames(mat_abund_week)

rowSums(mat_abund_int)
        
min_n <- min(rowSums(mat_abund_int))
min_n

set.seed(123)
mat_raref <- rrarefy(mat_abund_int, sample = min_n)


alpha_week <- data.frame(
  Sample = rownames(mat_raref),
  Observed = specnumber(mat_raref),
  Chao1 = estimateR(mat_raref)["S.chao1", ],
  ACE = apply(mat_raref, 1, fossil::ACE),
  Shannon = diversity(mat_raref, index = "shannon"),
  Simpson = diversity(mat_raref, index = "simpson"),
  InvSimpson = diversity(mat_raref, index = "invsimpson"),
  Fisher = fisher.alpha(mat_raref)
)

alpha_week


library(openxlsx)
write.xlsx(alpha_week, "Indices_alfa_diversidad_morfologico.xlsx")


# recuperar Temporada y Semana
alpha_week <- alpha_week %>%
  separate(Sample, into = c("Temporada", "Semana"), sep = "_")

library(ggplot2)

ggplot(alpha_week, aes(x = Temporada, y = Shannon)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_minimal()



by(alpha_week$Fisher, alpha_week$Temporada, shapiro.test)

car::leveneTest(Fisher ~ Temporada, data = alpha_week)


anova_chao1 <- aov(Chao1 ~ Temporada, data = alpha_week)
summary(anova_chao1)


anova_shannon <- aov(Shannon ~ Temporada, data = alpha_week)
summary(anova_shannon)


anova_simpson <- aov(Simpson ~ Temporada, data = alpha_week)
summary(anova_simpson)





kruskal_observed <- kruskal.test(Observed ~ Temporada, data = alpha_week)
kruskal_observed

kruskal_ace <- kruskal.test(ACE ~ Temporada, data = alpha_week)
kruskal_ace

kruskal_invsimpson <- kruskal.test(InvSimpson ~ Temporada, data = alpha_week)
kruskal_invsimpson

kruskal_fisher <- kruskal.test(Fisher ~ Temporada, data = alpha_week)
kruskal_fisher


### TABLA SIN NORMALIZAR CON RAREFACCION ###

# Observed richness
observed <- specnumber(mat_abund_int)

# Shannon
shannon <- diversity(mat_abund_int, index = "shannon")

# Simpson (1 - D)
simpson <- diversity(mat_abund_int, index = "simpson")

# Inverse Simpson
invsimpson <- diversity(mat_abund_int, index = "invsimpson")

# Chao1
chao1 <- estimateR(mat_abund_int)["S.chao1", ]

# ACE
ace <- apply(mat_abund_int, 1, fossil::ACE)

# Fisher alpha
fisher <- fisher.alpha(mat_abund_int)

# Tabla final
alpha <- data.frame(
  Sample = rownames(mat_abund_int),
  Observed = observed,
  Chao1 = chao1,
  ACE = ace,
  Shannon = shannon,
  Simpson = simpson,
  InvSimpson = invsimpson,
  Fisher = fisher
)

alpha




### ============================================================
###     ANÁLISIS DE DIVERSIDAD BETA (PCoA – Bray–Curtis)
###     Conteo morfológico
### ============================================================

library(vegan)
library(ggplot2)
library(ggrepel)

# Calcular matriz de disimilitud Bray–Curtis a partir del conteo morfológico
bray_dist_morf <- vegdist(mat_raref, method = "bray")

# Análisis de Coordenadas Principales (PCoA)
# Se extraen los dos primeros ejes para visualización
pcoa_res <- cmdscale(bray_dist_morf, eig = TRUE, k = 2)

# Convertir las coordenadas en un data.frame
pcoa_points <- as.data.frame(pcoa_res$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")

# Asociar cada punto con su respectiva temporada
pcoa_points$Temporada <- rownames(pcoa_points)

# Gráfico PCoA
ggplot(pcoa_points, aes(PCoA1, PCoA2, color = Temporada)) +
  geom_point(size = 5) +
  scale_color_manual(values = c(
    "Invierno"   = "#5DADE2",
    "Primavera"  = "#58D68D",
    "Verano"     = "#F4D03F",
    "Otoño"      = "#CA6F1E"
  )) +
  theme_minimal(base_size = 13) +
  labs(
    title = "PCoA Bray–Curtis (conteo morfológico)",
    x = paste0("pCo1 (", round(pcoa_res$eig[1] / sum(pcoa_res$eig) * 100, 1), "%)"),
    y = paste0("pCo2 (", round(pcoa_res$eig[2] / sum(pcoa_res$eig) * 100, 1), "%)")
  ) +
  theme(plot.title = element_text(face = "bold"))


# PERMANOVA para evaluar diferencias entre temporadas
library(vegan)

# --- 1. Crear tabla de metadata ---
# Asegúrate de que 'mat_raref' tiene filas = muestras
meta <- data.frame(
  Sample = rownames(mat_raref),
  Temporada = c("Invierno", "Invierno", "Primavera", "Primavera", 
                "Verano", "Verano", "Otoño", "Otoño") # Ajusta según tus muestras
)

# --- 2. Calcular matriz de disimilitud Bray-Curtis ---
bray_dist_morf <- vegdist(mat_raref, method = "bray")

# --- 3. PERMANOVA ---
permanova_morf <- adonis2(bray_dist_morf ~ Temporada, data = meta, permutations = 999)
print(permanova_morf)






meta <- data.frame(
  Temporada = factor(c("Invierno", "Primavera", "Verano", "Otoño"))
)

rownames(meta) <- rownames(mat_raref)

all(rownames(mat_raref) == rownames(meta))


permanova <- adonis2(
  bray_dist_morf ~ Temporada,
  data = meta,
  permutations = 999
)

permanova



#####
library(dplyr)
library(tidyr)

df_week <- conteo_generos %>%
  mutate(SampleID = paste(Temporada, Semana, sep = "_")) %>%
  select(SampleID, Temporada, Hongos, Total)


mat_week <- df_week %>%
  pivot_wider(
    names_from = Hongos,
    values_from = Total,
    values_fill = 0
  )

mat_week_clean <- mat_week %>%
  column_to_rownames("SampleID") %>%   # SampleID → rownames
  select(-Temporada) %>%               # quitar metadata
  mutate(across(everything(), as.numeric))



library(vegan)

min_n_week <- min(rowSums(mat_week_clean))
min_n_week
set.seed(123)

mat_week_raref <- rrarefy(mat_week_clean, sample = min_n_week)

bray_week <- vegdist(mat_week_raref, method = "bray")

pcoa_week <- cmdscale(bray_week, eig = TRUE, k = 2)

pcoa_week_df <- as.data.frame(pcoa_week$points)
colnames(pcoa_week_df) <- c("PCoA1", "PCoA2")
pcoa_week_df$SampleID <- rownames(pcoa_week_df)



pcoa_week_df <- pcoa_week_df %>%
  mutate(
    Temporada = sub("_.*", "", SampleID),
    Semana = sub(".*_", "", SampleID)
  )


ggplot(pcoa_week_df, aes(PCoA1, PCoA2, color = Temporada)) +
  geom_point(size = 4, alpha = 0.9) +
  stat_ellipse(
    aes(group = Temporada),
    type = "t",
    linetype = 2,
    level = 0.95,
    linewidth = 1
  ) +
  scale_color_manual(values = c(
    "Invierno"   = "#5DADE2",
    "Primavera"  = "#58D68D",
    "Verano"     = "#F4D03F",
    "Otoño"      = "#CA6F1E"
  )) +
  theme_minimal(base_size = 13) +
  labs(
    title = "PCoA Bray–Curtis por semanas (conteo morfológico)",
    x = paste0("PCoA1 (", round(pcoa_week$eig[1] / sum(pcoa_week$eig) * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(pcoa_week$eig[2] / sum(pcoa_week$eig) * 100, 1), "%)")
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  )



meta_week <- data.frame(
  Temporada = pcoa_week_df$Temporada,
  row.names = pcoa_week_df$SampleID
)

adonis2(
  bray_week ~ Temporada,
  data = meta_week,
  permutations = 999
)




# Análisis NMDS Bray-Curtis

set.seed(123)

nmds_res <- metaMDS(mat_raref, distance = "bray", k = 2)

nmds_points <- as.data.frame(nmds_res$points)
colnames(nmds_points) <- c("NMDS1", "NMDS2")
nmds_points$Temporada <- rownames(nmds_points)

ggplot(nmds_points, aes(NMDS1, NMDS2, color = Temporada)) +
  geom_point(size = 5) +
  scale_color_manual(values = c(
    "Invierno"   = "#5DADE2",
    "Primavera"  = "#58D68D",
    "Verano"     = "#F4D03F",
    "Otoño"      = "#CA6F1E"
  )) +
  theme_minimal(base_size = 13) +
  labs(
    title = paste0("NMDS Bray–Curtis (conteo morfológico, stress = ",
                   round(nmds_res$stress, 3), ")")
  ) +
  theme(plot.title = element_text(face = "bold"))




