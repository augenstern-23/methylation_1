
#awk 'NR==1 || $NF == "CpG"' GSE49828_WGBS_ICM_methylation_calling.bed.txt > ICM_filtered_CpG_only.txt
 
#### Datei auf CpGs in letzter Spalte filtern ####

# Pfade
input_file <- "data/raw/GSE49828_WGBS_ICM_methylation_calling.bed.txt"
output_file <- "data/processed/icm_cpg_only.txt"

# Bash-Befehl aufrufen
cmd <- paste0("awk 'NR==1 || $NF==\"CpG\"' ", input_file, " > ", output_file)
system(cmd)

message("CpG-Zeilen extrahiert und gespeichert in: ", output_file)



ICM_methylation <- read.csv("data/raw/ICM_filtered_CpG_only.txt", sep = "\t", header = TRUE)

# an chr1 ausprobieren
chr1_sample_table <- meth_df[meth_df[[1]] == "chr1", ]
chr1_icm <- ICM_methylation[ICM_methylation[[1]] == "chr1", ]


# Sicherstellen, dass Spaltennamen stimmen
colnames(chr1_icm) <- tolower(colnames(chr1_icm))
colnames(chr1_sample_table) <- tolower(colnames(chr1_sample_table))

# Neue Spalte vorbereiten
chr1_sample_table$icm_methylation <- NA_real_

# Loop über alle Zeilen in chr1_sample_table
for (i in seq_len(nrow(chr1_sample_table))) {
  region_chr   <- chr1_sample_table$chr[i]
  region_start <- chr1_sample_table$start[i]
  region_end   <- chr1_sample_table$end[i]
  
  # Filtere passende ICM-Positionen
  match_rows <- chr1_icm[
    chr1_icm$x.chr == region_chr &
      chr1_icm$pos >= region_start &
      chr1_icm$pos <= region_end,
  ]
  
  # Mittelwert berechnen, wenn Daten vorhanden sind
  if (nrow(match_rows) > 0) {
    chr1_sample_table$icm_methylation[i] <- mean(match_rows$metrate, na.rm = TRUE)
  }
}

# nur 23 matches --> falsches Referenzgenom?
