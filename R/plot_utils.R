# plot_utils.R
# Hilfsfunktionen zur Visualisierung von Methylierungsdaten
# Autor: Pauline Hübner – zuletzt bearbeitet: 05.06.2025

library(ggplot2)
library(pheatmap)
library(matrixStats)

#' Hole alle CpGs aus meth_df, die innerhalb einer gegebenen Region liegen
#'
#' @param meth_df DataFrame mit Methylierungsdaten (inkl. chr, start, end, *_frac-Spalten)
#' @param region_row Ein einzelnes DataFrame-Row-Objekt (z. B. region_df[i, ])
#' @return Ein DataFrame mit allen CpGs in dieser Region
get_region_cpgs <- function(meth_df, region_row) {
  meth_df %>%
    filter(
      chr == region_row$chr,
      start >= region_row$start,
      end <= region_row$end
    )
}

#' Erstelle Boxplot der Methylierungsverteilung für eine Region
#'
#' @param matches CpG-Treffer innerhalb der Region
#' @param region Genomische Region (einzelne Zeile aus region_df)
#' @param sample_name_map Umbenennung von Spaltennamen in Klartextnamen
#' @return Boxplot-Objekt (ggplot2) oder NULL, wenn keine Daten vorhanden
plot_region_methylation <- function(matches, region, sample_name_map) {
  frac_cols <- names(matches)[grepl("_frac$", names(matches))]
  
  if (length(frac_cols) == 0 || nrow(matches) == 0) {
    message("⚠️  Keine CpG-Daten in dieser Region – Boxplot wird übersprungen.")
    return(NULL)
  }
  
  # Daten ins lange Format bringen – nötig für ggplot
  frac_df_long <- matches[, frac_cols, drop = FALSE] %>%
    tidyr::pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = "methylation"
    ) %>%
    dplyr::mutate(
      sample = dplyr::recode(sample, !!!sample_name_map)
    )
  
  # Boxplot zeichnen
  p <- ggplot(frac_df_long, aes(x = sample, y = methylation)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(
      title = sprintf("Region %s:%d-%d", region$chr, region$start, region$end),
      x = "Sample",
      y = "Methylierungsgrad"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  return(p)
}

#' Extrahiere und filtere Methylierungsmatrix für Heatmap
#'
#' @param df DataFrame mit *_frac_mean-Spalten und Koordinaten in den ersten 3 Spalten (chr, start, end)
#' @return Gefilterte Methylierungsmatrix mit Koordinaten als Zeilennamen
prepare_methylation_matrix <- function(df) {
  message("Eingabedaten: ", nrow(df), " Zeilen, ", ncol(df), " Spalten")
  
  # Koordinaten-Spalte erzeugen
  coords <- paste0(df[[1]], ":", df[[2]], "-", df[[3]])
  
  # Extrahiere *_frac_mean-Spalten als Matrix
  meth_mat <- as.matrix(df[, grep("_frac_mean$", colnames(df))])
  rownames(meth_mat) <- coords  # Koordinaten als Zeilennamen
  
  message("Extrahierte Methylierungs-Matrix: ", nrow(meth_mat), " Zeilen, ", ncol(meth_mat), " Spalten")
  
  n_before_na <- nrow(meth_mat)
  meth_mat <- meth_mat[complete.cases(meth_mat), ]
  n_after_na <- nrow(meth_mat)
  message("Entfernte Zeilen mit NA: ", n_before_na - n_after_na)
  
  n_before_var <- nrow(meth_mat)
  row_vars <- apply(meth_mat, 1, var)
  meth_mat <- meth_mat[row_vars > 0, ]
  n_after_var <- nrow(meth_mat)
  message("Entfernte Zeilen mit Varianz = 0: ", n_before_var - n_after_var)
  
  message("Final: ", nrow(meth_mat), " Zeilen, ", ncol(meth_mat), " Spalten")
  
  return(meth_mat)
}



#' Benenne Spalten der Methylierungsmatrix um
#'
#' @param meth_mat Matrix mit Spaltennamen (Samples)
#' @param sample_name_map Named Vector zur Umbenennung
#' @return Matrix mit umbenannten Spalten
rename_sample_columns <- function(meth_mat, sample_name_map) {
  colnames(meth_mat) <- sample_name_map[colnames(meth_mat)]
  return(meth_mat)
}

#' Zeichne Heatmap der Methylierungsdaten
#'
#' @param meth_mat Gefilterte Methylierungsmatrix
#' @param filename Optionaler Dateiname zur Speicherung
plot_methylation_heatmap <- function(meth_mat, filename = NULL, plot_title = NULL) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  
  # Plot-Parameter
  plot_args <- list(
    meth_mat,
    scale = "row",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "correlation",
    show_rownames = FALSE,
    fontsize_col = 10,
    main = plot_title
  )
  
  # Speichern je nach Endung
  if (!is.null(filename)) {
    if (grepl("\\.pdf$", filename)) {
      plot_args$filename <- filename
      do.call(pheatmap::pheatmap, plot_args)
    } else if (grepl("\\.png$", filename)) {
      png(filename, width = 1000, height = 800, res = 150)
      do.call(pheatmap::pheatmap, plot_args)
      dev.off()
    } else {
      stop("❌ Nur .pdf oder .png werden als Ausgabeformate unterstützt.")
    }
    cat("Heatmap gespeichert unter:", filename, "\n")
  } else {
    do.call(pheatmap::pheatmap, plot_args)
  }
}


#' Extrahiere CpG-Matrix aus rohen *_frac-Spalten
#'
#' @param df DataFrame mit *_frac-Spalten
#' @return Gefilterte Methylierungsmatrix
prepare_cpg_matrix <- function(df) {
  message("Eingabedaten: ", nrow(df), " Zeilen, ", ncol(df), " Spalten")
  
  meth_mat <- as.matrix(df[, grep("_frac$", colnames(df))])
  message("Extrahierte CpG-Matrix: ", nrow(meth_mat), " Zeilen, ", ncol(meth_mat), " Spalten")
  
  n_before_na <- nrow(meth_mat)
  meth_mat <- meth_mat[complete.cases(meth_mat), ]
  n_after_na <- nrow(meth_mat)
  message("Entfernte Zeilen mit NA: ", n_before_na - n_after_na)
  
  n_before_var <- nrow(meth_mat)
  row_vars <- apply(meth_mat, 1, var)
  meth_mat <- meth_mat[row_vars > 0, ]
  n_after_var <- nrow(meth_mat)
  message("Entfernte Zeilen mit Varianz = 0: ", n_before_var - n_after_var)
  
  message("Final: ", nrow(meth_mat), " Zeilen, ", ncol(meth_mat), " Spalten")
  
  return(meth_mat)
}

# Methylierungs-PCA

# Lade und bereite Methylierungsdaten vor
read_methylation_data <- function(filepath) {
  df <- read.csv(filepath, sep = "\t", header = TRUE)
  df <- dplyr::rename_with(df, tolower)
  return(df)
}

# Entferne unvollständige Zeilen
filter_complete_cases <- function(df) {
  n_before <- nrow(df)
  df_clean <- df[complete.cases(df), ]
  n_after <- nrow(df_clean)
  message("Entfernte Zeilen mit NA: ", n_before - n_after)
  message("Verbleibende Zeilen: ", n_after)
  return(df_clean)
}

# Wähle Top-N variabelste CpGs mit Mindest-Coverage
select_top_variable_cpgs <- function(df, n_top, min_coverage = 10) {
  frac_cols <- grep("_frac$", colnames(df), value = TRUE)
  cov_cols <- gsub("_frac$", "_cov", frac_cols)
  cov_data <- df[, cov_cols]
  
  # Coverage-Filter
  coverage_ok <- cov_data >= min_coverage
  keep_rows <- rowSums(coverage_ok) == length(cov_cols)
  df_filt <- df[keep_rows, ]
  
  message("Coverage-Filter Top-CpGs: ", sum(keep_rows), " / ", nrow(df), " CpGs behalten")
  
  # Methylierungsmatrix extrahieren
  meth_mat <- as.matrix(df_filt[, frac_cols])
  row_vars <- rowVars(meth_mat, na.rm = TRUE)
  top_idx <- order(row_vars, decreasing = TRUE)[1:min(n_top, nrow(meth_mat))]
  
  return(list(
    meth = meth_mat[top_idx, ],
    meta = df_filt[top_idx, c("chr", "start", "end")]
  ))
}

# Wähle alle CpGs mit Mindest-Coverage
select_all_cpgs_with_coverage <- function(df, min_coverage = 10) {
  frac_cols <- grep("_frac$", colnames(df), value = TRUE)
  cov_cols <- gsub("_frac$", "_cov", frac_cols)
  cov_data <- df[, cov_cols]
  
  coverage_ok <- cov_data >= min_coverage
  keep_rows <- rowSums(coverage_ok) == length(cov_cols)
  df_filt <- df[keep_rows, ]
  
  message("Coverage-Filter Alle-CpGs: ", sum(keep_rows), " / ", nrow(df), " CpGs behalten")
  
  meth_mat <- as.matrix(df_filt[, frac_cols])
  
  return(list(
    meth = meth_mat,
    meta = df_filt[, c("chr", "start", "end")]
  ))
}

# Entferne Proben (Spalten) mit Varianz = 0 (inkl. Prüfung auf numerische Spalten)
filter_zero_variance_samples <- function(mat) {
  # Sicherstellen, dass es sich um eine numerische Matrix handelt
  mat <- as.matrix(mat)
  is_numeric <- apply(mat, 2, is.numeric)
  mat <- mat[, is_numeric, drop = FALSE]
  
  # Varianz berechnen (NA-sicher)
  sample_vars <- apply(mat, 2, var, na.rm = TRUE)
  
  # Filtern
  n_before <- length(sample_vars)
  mat_filtered <- mat[, sample_vars > 0, drop = FALSE]
  n_after <- ncol(mat_filtered)
  
  message("Entfernte Proben mit Varianz = 0 oder NA: ", n_before - n_after)
  message("Verbleibende Proben: ", n_after)
  
  return(mat_filtered)
}

#' Führe PCA durch und speichere Plot
#'
#' @param mat Matrix mit Methylierungswerten (Zeilen = CpGs, Spalten = Samples)
#' @param sample_mapping Named Vector zur Umbenennung der Spaltennamen
#' @param filename Dateipfad zur Speicherung des PCA-Plots (PDF)
#' @return Ergebnisobjekt von prcomp()
run_pca <- function(mat, 
                    sample_mapping = NULL, 
                    filename = "results/plots/pca_plot.pdf", 
                    plot_title = "PCA der Methylierungsmatrix") {
  
  # Zielordner anlegen
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  
  # Matrix vorbereiten
  mat <- as.matrix(mat)
  
  # Proben & CpGs mit Varianz = 0 oder NA entfernen
  sample_vars <- apply(mat, 2, var, na.rm = TRUE)
  mat <- mat[, sample_vars > 0, drop = FALSE]
  row_vars <- apply(mat, 1, var, na.rm = TRUE)
  mat <- mat[row_vars > 0, , drop = FALSE]
  
  if (ncol(mat) < 2 || nrow(mat) < 2) {
    stop("❌ Zu wenig valide Daten (n < 2) für PCA.")
  }
  
  # PCA berechnen
  pca_result <- prcomp(t(mat), scale. = TRUE)
  message("PCA erfolgreich berechnet!")
  
  # Daten für Plot vorbereiten
  pc_df <- as.data.frame(pca_result$x[, 1:2])
  pc_df$sample <- rownames(pc_df)
  pc_df$label <- if (!is.null(sample_mapping)) sample_mapping[pc_df$sample] else pc_df$sample
  
  # Farben manuell definieren
  group_colors <- c(
    "naive" = "#e41a1c",               # rot
    "reprimed_formative" = "#fb9a99",  # rosa
    "neuron_d26_reprimed" = "#f781bf", # pink
    "primed" = "#377eb8",              # blau
    "neuron_d26_primed" = "#a6cee3"    # hellblau
  )
  
  # Plot
  p <- ggplot(pc_df, aes(x = PC1, y = PC2, color = label)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(aes(label = label), size = 4, box.padding = 0.4) +
    scale_color_manual(values = group_colors) +
    labs(
      title = plot_title,
      x = "PC1", y = "PC2",
      color = "Zelltyp"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.margin = margin(t = 20, r = 30, b = 40, l = 40)
    )
  
  # Plot speichern
  if (grepl("\\.pdf$", filename)) {
    ggsave(filename, plot = p, width = 7.5, height = 6.5, units = "in", dpi = 300,
           limitsize = FALSE, useDingbats = FALSE)
  } else {
    ggsave(filename, plot = p, width = 7.5, height = 6.5, units = "in", dpi = 300,
           limitsize = FALSE)
  }
  
  message("PCA-Plot gespeichert unter: ", filename)
  
  return(pca_result)
}
