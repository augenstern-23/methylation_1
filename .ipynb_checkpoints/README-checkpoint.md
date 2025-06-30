# methylation_1

A reproducible R-based analysis pipeline for DNA methylation data.  
This project uses modular R scripts and a Conda environment to perform key analyses such as heatmaps, PCA, regional averaging, and CpG methylation visualization.

## Features

- Preprocessing of methylation data (e.g., hg19 table creation)
- Averaging methylation levels over genomic regions
- PCA of most variable methylation sites
- Generation of heatmaps and methylation plots
- Modular script structure in `scripts/`
- Reusable utility functions in `R/`
- Reproducible environment via Conda

## Installation

1. **Install Conda** (Miniconda or Anaconda)
2. **Create environment:**
   ```bash
   mamba env create -f env_methylation.yml -n env_methylation
   conda activate methylation
   ```

3. **(Optional)** Open the project in RStudio:
   - `methylation_1.Rproj`

## Usage

Run the scripts in `scripts/` for each step of the analysis:

```r
# Example: Compute average methylation
source("scripts/02_compute_average_methylation.R")

# Example: Create CpG heatmap
source("scripts/04_heatmap_cpgs.R")
```

The `R/` folder contains utility functions for plots, summaries, and region-level operations.

## Project Structure

```
.
├── R/                       # Utility functions
├── scripts/                # Analysis scripts
├── env_methylation.yml     # Conda environment file
├── methylation_1.Rproj     # RStudio project file
└── LICENSE
```

## License

This project is licensed under the MIT License – see [LICENSE](./LICENSE).

## Author

Developed by [augenstern-23](https://github.com/augenstern-23)
