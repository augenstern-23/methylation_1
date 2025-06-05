# R/utils.R

#### Funktion: Extrahiere chr, start, end aus einer Spalte im Format "chr1:567502-567715"
extract_hg19_coords <- function(hg19_vector) {
  # Erwartet: Vektor mit Strings wie "chr1:567502-567715"
  parts <- strcapture(
    pattern = "^([^:]+):([0-9]+)-([0-9]+)$",
    x = hg19_vector,
    proto = list(chr = character(), start = integer(), end = integer())
  )
  return(parts)
}

