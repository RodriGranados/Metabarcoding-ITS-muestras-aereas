library(phyloseq)
library(tidyverse)

# --- 1. Definir función para exportar ---
export_taxa_abundance <- function(ps_1, rank, relative = TRUE, outdir = "Figuras y Tablas/") {
  # Crear carpeta de salida si no existe
  if (!dir.exists(outdir)) dir.create(outdir)
  
  # Agrupar el phyloseq por nivel taxonómico
  ps_glom <- tax_glom(ps_1, taxrank = rank)
  
  # Convertir a abundancias relativas si se solicita
  if (relative) {
    ps_glom <- transform_sample_counts(ps_glom, function(x) x / sum(x))
  }
  
  # Convertir a formato largo (melted)
  df <- psmelt(ps_glom)
  
  # Calcular abundancia promedio por muestra y taxón
  df_summary <- df %>%
    group_by(Sample, !!sym(rank)) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    pivot_wider(names_from = Sample, values_from = Abundance, values_fill = 0)
  
  # Guardar como archivo TSV
  outfile <- file.path(outdir, paste0("abundancia_", tolower(rank), ".tsv"))
  write.table(df_summary, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("✅ Exportado: ", outfile)
  return(df_summary)
}

# --- 2. Obtener los niveles taxonómicos disponibles ---
niveles <- rank_names(ps_1)
print(niveles)

# --- 3. Exportar todos los niveles disponibles ---
tablas_abundancia <- lapply(niveles, function(rank) {
  export_taxa_abundance(ps_1, rank = rank, relative = TRUE)
})

names(tablas_abundancia) <- niveles
