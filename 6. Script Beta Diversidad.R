### ============================================================
###     ANÁLISIS DE DIVERSIDAD BETA
### ============================================================

library(vegan)
library(tidyverse)
library(ggrepel)

# --- 1. Usar la matriz de abundancia rarefactada ---
otu_beta <- t(otu_rarefied)   # muestras en filas

# --- 2. Calcular distancias de disimilitud ---
## Bray-Curtis (más usada para abundancia)
bray_dist <- vegdist(otu_beta, method = "bray")


# --- 3. Análisis de ordenación (visualización) ---

## a) PCoA (Análisis de coordenadas principales)
pcoa_res <- cmdscale(bray_dist, eig = TRUE, k = 2)
pcoa_points <- as.data.frame(pcoa_res$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")
pcoa_points$Sample <- rownames(pcoa_points)

# Añadir temporadas
pcoa_points <- pcoa_points %>%
  left_join(alpha_metrics %>% select(Sample, Season), by = "Sample")

# --- 4. Gráfico PCoA ---
ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = Season)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c(
    "Invierno" = "#5DADE2",
    "Primavera" = "#58D68D",
    "Verano" = "#F4D03F",
    "Otoño" = "#CA6F1E"
  )) +
  theme_minimal(base_size = 13) +
  labs(
    title = "PCoA basada en distancia Bray-Curtis (ITS)",
    x = paste0("Eje 1 (", round(pcoa_res$eig[1] / sum(pcoa_res$eig) * 100, 1), "%)"),
    y = paste0("Eje 2 (", round(pcoa_res$eig[2] / sum(pcoa_res$eig) * 100, 1), "%)")
  )


# --- 5. Análisis NMDS (Non-metric Multidimensional Scaling) ---
set.seed(123)  # Para reproducibilidad
nmds_res <- metaMDS(bray_dist, k = 2, trymax = 100)

# Extraer coordenadas NMDS
nmds_points <- as.data.frame(nmds_res$points)
colnames(nmds_points) <- c("NMDS1", "NMDS2")
nmds_points$Sample <- rownames(nmds_points)

# Añadir temporadas
nmds_points <- nmds_points %>%
  left_join(alpha_metrics %>% select(Sample, Season), by = "Sample")


# --- Define tu propio título aquí ---
custom_title <- "NMDS basada en distancia Bray-Curtis (ITS)"

# --- 6. Gráfico NMDS ---
ggplot(nmds_points, aes(x = NMDS1, y = NMDS2, color = Season)) +
  geom_point(size = 4, alpha = 0.8) +
  # Etiquetas con nombres de muestra
  geom_text_repel(aes(label = Sample), size = 3.5, show.legend = FALSE, max.overlaps = 100) +
  scale_color_manual(values = c(
    "Invierno" = "#5DADE2",
    "Primavera" = "#58D68D",
    "Verano"   = "#F4D03F",
    "Otoño"    = "#CA6F1E"
  )) +
  theme_minimal(base_size = 13) +
  labs(
    title = paste0(custom_title, " (stress = ", round(nmds_res$stress, 3), ")"),
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14)
  )

# --- 7. Análisis estadístico: PERMANOVA ---
## Para probar si las comunidades difieren significativamente entre temporadas
permanova <- adonis2(bray_dist ~ Season, data = pcoa_points, permutations = 999)
print(permanova)

# Convertir el resultado a un data frame
permanova_df <- as.data.frame(permanova)

# Exportar a Excel
library(openxlsx)

write.xlsx(permanova_df, file = "permanova_results.xlsx", row.names = TRUE)






## --- NMDS ---
set.seed(123)  # importante para reproducibilidad

nmds_res <- metaMDS(
  otu_beta,
  distance = "bray",
  k = 2,
  trymax = 100,
  autotransform = FALSE
)

nmds_res$stress


nmds_points <- as.data.frame(nmds_res$points)
nmds_points$Sample <- rownames(nmds_points)


# Añadir metadata (temporada)
nmds_points <- nmds_points %>%
  left_join(alpha_metrics %>% select(Sample, Season), by = "Sample")

ggplot(nmds_points, aes(x = MDS1, y = MDS2, color = Season)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c(
    "Invierno" = "#5DADE2",
    "Primavera" = "#58D68D",
    "Verano" = "#F4D03F",
    "Otoño" = "#CA6F1E"
  )) +
  theme_minimal(base_size = 13) +
  labs(
    title = "NMDS basada en distancia Bray-Curtis (ITS)",
    subtitle = paste("Stress =", round(nmds_res$stress, 3)),
    x = "NMDS1",
    y = "NMDS2"
  )


