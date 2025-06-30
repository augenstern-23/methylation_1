# PCA Workflow für alle CpGs und Top 10.000 variabelste CpGs

library(ggplot2)
source("R/plot_utils.R")  # enthält run_pca() und weitere Hilfsfunktionen

# 1. Daten einlesen
meth_df <- read_methylation_data("data/raw/combined_coverage_methylation_matrix.tsv")

# 2. Zeilen mit NA entfernen
meth_clean <- filter_complete_cases(meth_df)

# 3a. PCA auf Top 10.000 variable CpGs mit Mindest-Coverage
res_top <- select_top_variable_cpgs(meth_clean, n_top = 1000000, min_coverage = 10)
meth_mat_top <- res_top$meth
coord_info_top <- res_top$meta
message("🔢 Matrixgröße (Top 1M CpGs): ", dim(meth_mat_top)[1], " x ", dim(meth_mat_top)[2])
meth_top_filtered <- filter_zero_variance_samples(meth_mat_top)

# 4a. PCA auf Top CpGs ausführen
sample_mapping <- c(
  "x2_s13_frac"  = "naive",
  "x38_s14_frac" = "neuron_d26_primed",
  "x5_s15_frac"  = "primed",
  "x56_s16_frac" = "neuron_d26_reprimed",
  "x8_s17_frac"  = "reprimed_formative"
)

pca_top_png <- run_pca(
  meth_top_filtered,
  sample_mapping,
  filename = "results/plots/pca_top1m_cpgs.png",
  plot_title = "PCA Top 1,000,000 variable CpGs"
  )

pca_top <- run_pca(
  meth_top_filtered,
  sample_mapping = sample_mapping,
  filename = "results/plots/pca_top100k_cpgs.pdf",
  plot_title = "PCA Top 100,000 variable CpGs"
)




# 3b. PCA auf allen CpGs mit Mindest-Coverage
res_all <- select_all_cpgs_with_coverage(meth_clean, min_coverage = 10)
meth_mat_all <- res_all$meth
coord_info_all <- res_all$meta
message("🔢 Matrixgröße (alle CpGs): ", dim(meth_mat_all)[1], " x ", dim(meth_mat_all)[2])

# 4b. PCA auf allen CpGs ausführen
pca_all <- run_pca(
  meth_mat_all,
  sample_mapping = sample_mapping,
  filename = "results/plots/pca_all_cpgs.pdf", 
  plot_title = "PCA All CpGs"
)


