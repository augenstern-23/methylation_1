#!/usr/bin/env Rscript

# scripts/plot_methylation_boxplots.R

# Lade Pakete
library(dplyr)
library(tidyr)
library(ggplot2)

# Lade Funktionen
source("R/plot_utils.R")

# ---- Daten einlesen ----

meth_df <- read.csv("data/raw/combined_coverage_methylation_matrix.tsv", sep = "\t", header = TRUE)
region_df <- read.csv("data/processed/hg19_coordinates.csv")

# Optional: Spaltennamen vereinheitlichen (z.B. chr statt Chr)
names(meth_df) <- tolower(names(meth_df))
names(region_df) <- tolower(names(region_df))

# ---- Sample-Namen umbenennen ----

sample_name_map <- c(
  "x2_s13_frac"  = "naive",
  "x38_s14_frac" = "neuron_d26_primed",
  "x5_s15_frac"  = "primed",
  "x56_s16_frac" = "neuron_d26_reprimed",
  "x8_s17_frac"  = "reprimed_formative"
)

# ---- Boxplots für erste 5 Regionen ----

for (i in 1:5) {
  region <- region_df[i, ]
  matches <- get_region_cpgs(meth_df, region)
  p <- plot_region_methylation(matches, region, sample_name_map)
  
  if (!is.null(p)) print(p)
}
