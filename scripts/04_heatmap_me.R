#!/usr/bin/env Rscript

source("R/plot_utils.R")

meth_summary_df <- read.csv("data/processed/methylation_summary.csv", header = TRUE)
meth_mat <- prepare_methylation_matrix(meth_summary_df)

# Spaltennamen umbenennen 
sample_name_map <- c(
  "x2_s13_frac_mean"  = "naive",
  "x38_s14_frac_mean" = "neuron_d26_primed",
  "x5_s15_frac_mean"  = "primed",
  "x56_s16_frac_mean" = "neuron_d26_reprimed",
  "x8_s17_frac_mean"  = "reprimed_formative"
)

meth_mat <- rename_sample_columns(meth_mat, sample_name_map)

# Speichern 
#plot_file <- "results/plots/me_methylation_heatmap.pdf"
#plot_methylation_heatmap(meth_mat, filename = plot_file)

# PDF mit Titel
plot_methylation_heatmap(meth_mat, filename = "results/plots/heatmap_me.pdf", plot_title = "CpGs in MEs")

# PNG mit Titel
plot_methylation_heatmap(meth_mat, filename = "results/plots/heatmap_me.png", plot_title = "CpGs in MEs")
