#' Durchschnittliche Methylierung/Coverage pro Region (GRanges-Version)
#'
#' @param meth_df CpG-Methylierungsdaten mit Spalten chr, start, end, sampleX_frac, sampleX_cov
#' @param region_df Regionendaten mit Spalten chr, start, end
#' @return DataFrame mit Durchschnittswerten je Region
#' 
compute_region_averages_gr <- function(meth_df, region_df) {
  library(GenomicRanges)
  library(dplyr)
  
  # Vereinheitlichung chr-Spalte
  meth_df$chr <- tolower(trimws(as.character(meth_df$chr)))
  region_df$chr <- tolower(trimws(as.character(region_df$chr)))
  
  # Optional: Standardchromosomen filtern
  allowed_chr <- paste0("chr", c(1:22, "x", "y"))
  meth_df <- meth_df %>% filter(chr %in% allowed_chr)
  
  # GRanges-Objekte erzeugen
  meth_gr <- GRanges(
    seqnames = meth_df$chr,
    ranges = IRanges(start = meth_df$start, end = meth_df$end),
    meth_df[, !(names(meth_df) %in% c("chr", "start", "end"))]
  )
  
  region_gr <- GRanges(
    seqnames = region_df$chr,
    ranges = IRanges(start = region_df$start, end = region_df$end)
  )
  
  # Overlaps suchen
  hits <- findOverlaps(region_gr, meth_gr)
  hits_df <- as.data.frame(hits)  # Spalten: queryHits (region), subjectHits (CpG)
  
  # Verbinde CpG-Daten mit Regions-ID
  meth_hits <- as.data.frame(mcols(meth_gr[hits$subjectHits]))
  meth_hits$region_id <- hits$queryHits
  
  # Mittelwerte pro Region berechnen
  result <- meth_hits %>%
    group_by(region_id) %>%
    summarise(across(everything(), ~mean(.x, na.rm = TRUE))) %>%
    left_join(region_df %>% mutate(region_id = row_number()), by = "region_id") %>%
    select(chr, start, end, everything(), -region_id)
  
  return(result)
}
