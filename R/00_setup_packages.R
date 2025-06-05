# scripts/00_setup_packages.R

# Basispakete aus CRAN
packages_cran <- c("dplyr", "readr", "ggplot2", "tidyr")

# Bioconductor-Pakete
packages_bioc <- c("GenomicRanges", "SummarizedExperiment")

# Installiere CRAN-Pakete, falls nicht vorhanden
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

lapply(packages_cran, install_if_missing)

# Bioconductor-Manager installieren
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Bioconductor-Pakete installieren
for (pkg in packages_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}
