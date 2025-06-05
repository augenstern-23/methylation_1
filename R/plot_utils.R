# R/plot_utils.R

#' Hole alle CpGs aus meth_df, die innerhalb einer gegebenen Region liegen
#'
#' @param meth_df DataFrame mit Methylierungsdaten (inkl. chr, start, end, *_frac-Spalten)
#' @param region_row Ein einzelnes DataFrame-Row-Objekt (z. B. region_df[i, ])
#' @return Ein DataFrame mit allen CpGs in dieser Region

library(pheatmap)

get_region_cpgs <- function(meth_df, region_row) {
  meth_df %>%
    filter(
      chr == region_row$chr,
      start >= region_row$start,
      end <= region_row$end
    )
}

#' Erstelle Boxplot der Methylierungswerte für eine Region
#'
#' @param matches Ein DataFrame mit CpGs (z. B. Rückgabe von get_region_cpgs)
#' @param region Die Region (einzelne Zeile aus region_df)
#' @param sample_name_map Ein Named Vector mit Umbenennungen: names = Spaltennamen, values = Klartext
#' @return Ein ggplot2-Plot-Objekt oder NULL, falls keine Daten vorhanden sind
plot_region_methylation <- function(matches, region, sample_name_map) {
  frac_cols <- names(matches)[grepl("_frac$", names(matches))]
  
  if (length(frac_cols) == 0 || nrow(matches) == 0) {
    message("Keine CpGs oder *_frac-Spalten vorhanden – übersprungen")
    return(NULL)
  }
  
  # Long-Format für ggplot
  frac_df_long <- matches[, frac_cols, drop = FALSE] %>%
    tidyr::pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = "methylation"
    ) %>%
    dplyr::mutate(
      sample = dplyr::recode(sample, !!!sample_name_map)
    )
  
  # Boxplot erstellen
  p <- ggplot2::ggplot(frac_df_long, ggplot2::aes(x = sample, y = methylation)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.5) +
    ggplot2::labs(
      title = sprintf("Region %s:%d-%d", region$chr, region$start, region$end),
      x = "Sample",
      y = "Methylierungsgrad"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  
  return(p)
}


prepare_methylation_matrix <- function(df) {
  meth_mat <- as.matrix(df[ , grep("_frac_mean$", colnames(df))])
  meth_mat <- meth_mat[complete.cases(meth_mat), ]
  row_vars <- apply(meth_mat, 1, var)
  meth_mat <- meth_mat[row_vars > 0, ]
  return(meth_mat)
}

rename_sample_columns <- function(meth_mat, sample_name_map) {
  colnames(meth_mat) <- sample_name_map[colnames(meth_mat)]
  return(meth_mat)
}

library(pheatmap)

plot_methylation_heatmap <- function(meth_mat, filename = NULL) {
  if (!is.null(filename)) {
    dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
    pheatmap(meth_mat,
             scale = "row",
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "correlation",
             show_rownames = FALSE,
             fontsize_col = 10,
             filename = filename)  # <- direkte Speicherung durch pheatmap
    cat("✅ Heatmap gespeichert unter:", filename, "\n")
  } else {
    # interaktive Anzeige
    pheatmap(meth_mat,
             scale = "row",
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "correlation",
             show_rownames = FALSE,
             fontsize_col = 10)
  }
}

prepare_cpg_matrix <- function(df) {
  # Annahme: chr, start, end = Spalten 1–3, danach Methylierungswerte
  meth_mat <- as.matrix(df[ , grep("_frac$", colnames(df))])
  meth_mat <- meth_mat[complete.cases(meth_mat), ]
  row_vars <- apply(meth_mat, 1, var)
  meth_mat <- meth_mat[row_vars > 0, ]
  return(meth_mat)
}


