#' Berechnet mittlere Methylierung/Coverage je Region und Sample
#'
#' @param meth_df CpG-Tabelle mit chr, start, end, danach alle Sample-Spalten
#' @param region_df Regionstabelle mit chr, start, end
#' @return DataFrame mit allen Regionen + Mittelwerten pro Sample
#' 
summarise_region_methylation <- function(meth_df, region_df) {
  library(dplyr)
  
  # Standardisiere chr-Spalten (wichtig für sicheren Vergleich)
  meth_df$chr <- tolower(trimws(as.character(meth_df$chr)))
  region_df$chr <- tolower(trimws(as.character(region_df$chr)))
  
  result_list <- list()
  
  for (i in seq_len(nrow(region_df))) {
    region <- region_df[i, ]
    
    # Finde alle CpGs innerhalb der Region
    matches <- meth_df %>%
      filter(
        chr == region$chr,
        start >= region$start,
        end <= region$end
      )
    
    if (nrow(matches) > 0) {
      # Dynamisch: alle Spalten nach "end"
      end_index <- which(names(matches) == "end")
      sample_cols <- names(matches)[(end_index + 1):ncol(matches)]
      
      avg_values <- matches %>%
        summarise(
          across(
            .cols = all_of(sample_cols),
            .fns  = ~mean(.x, na.rm = TRUE),
            .names = "{.col}_mean"
          )
        )
    } else {
      # Wenn keine CpGs gefunden wurden → NA-Zeile mit allen Mittelwert-Spalten
      end_index <- which(names(meth_df) == "end")
      sample_cols <- names(meth_df)[(end_index + 1):ncol(meth_df)]
      avg_colnames <- paste0(sample_cols, "_mean")
      avg_values <- as.data.frame(matrix(NA, nrow = 1, ncol = length(avg_colnames)))
      colnames(avg_values) <- avg_colnames
    }
    
    # Kombiniere Region mit berechneten Mittelwerten
    result_list[[i]] <- bind_cols(region, avg_values)
  }
  
  # Endgültiges Ergebnis zusammenfügen
  result_df <- bind_rows(result_list)
  # Prüfen, ob alle Werte in den Methylierungsspalten (ab Spalte 4) NA sind
  keep_rows <- !apply(is.na(result_df[, 4:ncol(result_df)]), 1, all)
  # Nur die Zeilen behalten, die mindestens einen Wert haben
  clean_df <- result_df[keep_rows, ]
  
  return(clean_df)
}


filter_regions <- function(meth_df, region_df) {
  # Standardisiere chr-Spalten
  meth_df <- meth_df %>%
    mutate(chr = tolower(trimws(chr)))
  
  region_df <- region_df %>%
    mutate(chr = tolower(trimws(chr)))
  
  # Finde passende Zeilen in meth_df basierend auf region_df
  filtered_df <- semi_join(meth_df, region_df, by = c("chr", "start", "end"))
  
  return(filtered_df)
}

