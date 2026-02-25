library(tidyverse)
library(phyloseq)

# Leer tabla de abundancia


feature_table <- read.table("feature-table.tsv",
                            sep = "\t",       # separador tabulador
                            header = TRUE,    # usa la primera fila como nombres de columna
                            row.names = 1,    # usa la primera columna (#OTU ID) como nombres de fila
                            check.names = FALSE)  # mantiene nombres de muestras como están

# Leer tabla de taxonomía
taxonomy <- read.table("taxonomy.tsv", sep="\t", header = TRUE)


# Leer metadatos
metadata <- read.table("metadata.tsv", sep="\t", header=TRUE, row.names=1)

tax_split <- taxonomy %>%
  separate_wider_delim(
    cols = Taxon,
    delim = ";",
    names_sep = "_",    # genera columnas como Taxon_1, Taxon_2, etc.
    too_few = "align_start"
  )

##modificacion del tax split
tax_cols <- paste0("Taxon_", 1:8)

tax_split <- tax_split %>%
  rowwise() %>%
  mutate(
    # detectar dónde aparece por primera vez "Incertae_sedis"
    first_incert = which(grepl("Incertae_sedis", c_across(all_of(tax_cols))))[1],
    
    # reemplazar celdas posteriores por NA
    across(all_of(tax_cols), 
           ~ if (!is.na(first_incert) && 
                 which(tax_cols == cur_column()) >= first_incert) NA_character_ else .x)
  ) %>%
  ungroup() %>%
  select(-first_incert)   # remover columna auxiliar



tax_split <- tax_split %>%
  mutate(
    Taxon_7 = ifelse(
      is.na(Taxon_7) & !is.na(Taxon_6) & grepl("^g__", Taxon_6),
      paste0("s__", sub("^g__", "", Taxon_6), "_sp"),
      Taxon_7
    )
  )


tax_split_1 <- tax_split %>%
  mutate(
    across(
      Taxon_1:Taxon_7,
      ~ ifelse(
        is.na(.),
        NA,
        sub("^[a-z]__", "", .)
      )
    ),
    Taxon_7 = ifelse(
      is.na(Taxon_7),
      NA,
      gsub("_", " ", Taxon_7)
    )
  )

tax_split_1 <- tax_split_1 %>%
  mutate(
    Taxon_7 = ifelse(
      !is.na(Taxon_6) &
        grepl("Incertae_sedis", Taxon_6) &
        !is.na(Taxon_7) &
        grepl("^.+ sp$", Taxon_7),
      paste0(
        sub("_gen_Incertae_sedis", "", Taxon_6),
        "_sp_Incertae_sedis"
      ),
      Taxon_7
    )
  )


tax_split_1 <- tax_split_1 %>%
  mutate(
    Taxon_7 = ifelse(
      !is.na(Taxon_7) &
        grepl("\\bsp$", Taxon_7),
      sub("\\bsp$", "sp.", Taxon_7),
      Taxon_7
    )
  )



# Convertir a matriz y poner los IDs como nombres de fila
tax_mat <- as.matrix(tax_split[, -1])  # eliminamos la columna Feature ID
rownames(tax_mat) <- taxonomy$Feature.ID

# Crear objetos individuales
OTU <- otu_table(as.matrix(feature_table), taxa_are_rows=TRUE)
TAX <- tax_table(tax_mat_1)
META <- sample_data(metadata)

# Combinar en un objeto phyloseq
ps <- phyloseq(OTU, TAX, META)

sample_names(OTU)[1:5]
sample_names(META)[1:5]# Verificar
ps


colnames(tax_table(ps)) <- c(
  "Reino", "Phylum", "Clase", "Orden", "Familia", "Genero", "Especie", "Subespecie", "Confidence"
)


## Funcion para construir los graficos
plot_tax_level <- function(ps, rank="Orden", group="Temporada", top=10, threshold=0.03) {
  # Asegurar que el grupo existe en metadata
  if (!(group %in% colnames(sample_data(ps)))) {
    stop(paste("No se encontró la columna", group, "en los metadatos"))
  }
  
  # Obtener el orden de las muestras según la metadata
  orden_muestras <- unique(sample_data(ps)[[group]])
  
  # Procesar el phyloseq
  df <- ps %>%
    tax_glom(taxrank = rank) %>%
    transform_sample_counts(function(x) x / sum(x)) %>%
    psmelt() %>%
    group_by(!!sym(group), !!sym(rank)) %>%
    summarise(Abundance = mean(Abundance), .groups = "drop")
  
  # Top N taxones
  top_taxa <- df %>%
    group_by(!!sym(rank)) %>%
    summarise(mean_abund = mean(Abundance)) %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = top) %>%
    pull(!!sym(rank))
  
  # Agrupar los demás en "Otros"
  df <- df %>%
    mutate(!!rank := ifelse(!!sym(rank) %in% top_taxa, !!sym(rank), "Otros")) %>%
    group_by(!!sym(group), !!sym(rank)) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  # Reordenar niveles de taxones (Otros al final)
  taxa_order <- df %>%
    group_by(!!sym(rank)) %>%
    summarise(total = sum(Abundance)) %>%
    arrange(desc(total)) %>%
    pull(!!sym(rank))
  df[[rank]] <- factor(df[[rank]], levels = c(setdiff(taxa_order, "Otros"), "Otros"))
  
  # Reordenar muestras según metadata
  df[[group]] <- factor(df[[group]], levels = orden_muestras)
  
  # Añadir etiquetas de porcentaje (>3%)
  df <- df %>%
    mutate(label = ifelse(Abundance >= threshold, paste0(round(Abundance*100, 1), "%"), NA))
  
  # Colores
  colores <- c(RColorBrewer::brewer.pal(min(length(unique(df[[rank]])) - 1, 12), "Paired"), "grey70")
  n_taxa <- length(levels(df[[rank]]))
  n_sin_otros <- n_taxa - 1
  
  paleta_base <- if (n_sin_otros <= 12) {
    RColorBrewer::brewer.pal(n_sin_otros, "Paired")
  } else {
    colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_sin_otros)
  }
  
  colores <- c(paleta_base, "grey70")
  # Graficar
  ggplot(df, aes(x = !!sym(group), y = Abundance, fill = !!sym(rank))) +
    geom_bar(stat = "identity", color = "black") +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              size = 3,
              color = "black",
              na.rm = TRUE) +
    scale_fill_manual(values = colores) +
    labs(title = paste("Composición de hongos aereos por", rank),
         x = NULL, y = "Abundancia relativa") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(
        hjust = 0.5,       # centrar título
        face = "bold",     # negrita
        size = 16          # tamaño un poco mayor
      ),
      axis.text.x = element_text(
        angle = 0,         # horizontal
        hjust = 0.5,       # centrado
        vjust = 1
      ),
      legend.position = "right"
    ) +
    coord_cartesian(ylim = c(0, 1))
}

## Graficar con el siguiente ejemplo
plot_tax_level(ps, rank="Especie", group="Temporada", top=15, threshold=0.03)

