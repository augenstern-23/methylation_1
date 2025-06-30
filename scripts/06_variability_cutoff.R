source("R/plot_utils.R")  # enthält run_pca() und weitere Hilfsfunktionen


meth_df <- read_methylation_data("data/raw/combined_coverage_methylation_matrix.tsv")
# Zeilen mit NA entfernen
meth_clean <- filter_complete_cases(meth_df)



# 1. CpG-Matrix extrahieren
frac_cols <- grep("_frac$", colnames(meth_clean), value = TRUE)
meth_mat <- as.matrix(meth_clean[, frac_cols])

# 2. Nur CpGs mit vollständigen Werten verwenden
meth_mat <- meth_mat[complete.cases(meth_mat), ]

# 3. Zeilenvarianz berechnen (für jede CpG über alle Proben)
row_vars <- rowVars(meth_mat, na.rm = TRUE)

# 4. Sortieren der Varianzen (absteigend)
sorted_vars <- sort(row_vars, decreasing = TRUE)

# 5. Elbow-Plot zeichnen
plot(sorted_vars[1:50000], type = "l", log = "y",
     main = "Elbow Plot: CpG-Varianz (log-Skala)",
     xlab = "Top-N CpGs (nach Varianz)", ylab = "Varianz (log)")
abline(v = c(10000, 100000, 1000000), col = c("red", "blue", "green"), lty = 2)
legend("topright", legend = c("10k", "100k", "1M"),
       col = c("red", "blue", "green"), lty = 2)



pdf("results/plots/elbow_plot_cpg_variance.pdf", width = 7, height = 5)

plot(sorted_vars[1:50000], type = "l", log = "y",
     main = "Elbow Plot: CpG-Varianz (log-Skala)",
     xlab = "Top-N CpGs (nach Varianz)", ylab = "Varianz (log)")
abline(v = c(10000, 100000, 1000000), col = c("red", "blue", "green"), lty = 2)
legend("topright", legend = c("10k", "100k", "1M"),
       col = c("red", "blue", "green"), lty = 2)

dev.off()

Also, I played around a little and thought it might be better to use the CpGs with the highest variance in methylation between cell states for PCA, instead of including all sites.

To test this, I compared different cutoffs, meaning how many of the most variable CpGs we include in the analysis – I tried the top 10,000, 100,000, and 1 million CpGs (ranked by variance).
The idea is that CpGs with very low variance across samples are less informative for separating cell types and may just introduce noise.

Here’s what I found:
  
  With 10,000 CpGs, the PCA shows a clear and compact separation between cell types, with minimal noise.
With 1 million CpGs, the structure is still visible but the scale is much more stretched – probably due to many low-variance sites contributing mostly noise.
The 100,000-CpG version looked almost identical to the 10,000 one – until I realized I had accidentally reused the 10k matrix 😅 (fixed now).

I also plotted the cumulative variance across all CpGs and noticed that it rises very steeply in the beginning and then flattens out after ~10,000 CpGs. That supports the idea that most of the signal is concentrated in a relatively small subset of CpGs.

So for now, I’m using the top 10,000 most variable CpGs as input for the PCA and downstream plots.

What do you think about that cutoff? Would you also go with ~10k, or would you prefer including more sites?