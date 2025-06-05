#!/usr/bin/env Rscript 

source("R/summarise_region_methylation.R")
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)



meth_df <- read.csv("data/raw/combined_coverage_methylation_matrix.tsv", sep = "\t", header = TRUE) %>%
  rename_with(tolower)

region_df <- read.csv("data/processed/hg19_coordinates.csv") %>%
  rename_with(tolower)

# df speichern mit Overlap zwischen ME und Methylation
overlap_df <- filter_regions(meth_df, region)

# Auswertung durchführen
result_df <- summarise_region_methylation(meth_df, region_df)


# Ausgabe speichern 
write.csv(result_df, "data/processed/methylation_summary.csv", row.names = FALSE)
write.csv(overlap_df, "data/processed/MEs_methylation.csv", row.names = FALSE)
