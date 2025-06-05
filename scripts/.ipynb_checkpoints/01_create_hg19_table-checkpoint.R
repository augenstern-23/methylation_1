#!/usr/bin/env Rscript

# scripts/01_create_hg19_table.R

# Pakete laden
library(dplyr)
source("R/utils.R")

metastable_epiallele_list<- read.csv("data/processed/metastable_epiallele_list.csv")

# hg19-Spalte extrahieren, Neue Tabelle erstellen + speichern
hg19_coords <- extract_hg19_coords(metastable_epiallele_list$`Hg19 data`)
hg19_table <- bind_cols(hg19_coords)
print(head(hg19_table))

write.csv(hg19_table, "data/processed/hg19_coordinates.csv", row.names = FALSE)

