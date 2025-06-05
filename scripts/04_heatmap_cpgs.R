#!/usr/bin/env Rscript

source("R/plot_utils.R")

# Datei mit CpG-Daten (jede Zeile = CpG)
cpg_df <- read.csv("data/processed/MEs_methylation.csv", header = TRUE)

meth_mat <- prepare_cpg_matrix(cpg_df)

# Spaltennamen umbenennen 
sample_name_map <- c(
  "x2_s13_frac"  = "naive",
  "x38_s14_frac" = "neuron_d26_primed",
  "x5_s15_frac"  = "primed",
  "x56_s16_frac" = "neuron_d26_reprimed",
  "x8_s17_frac"  = "reprimed_formative"
)
colnames(meth_mat) <- sample_name_map[colnames(meth_mat)]

# Heatmap speichern
plot_file <- "results/plots/cpg_methylation_heatmap.pdf"
plot_methylation_heatmap(meth_mat, filename = plot_file)
