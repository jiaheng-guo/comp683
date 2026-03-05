# Optimized MERFISH Spatial Transcriptomics Reproduction

This repo contains a refactored, reproducible implementation of the analysis pipeline from **Ayata et al. (2025)** using MERFISH spatial transcriptomics data. The original code was provided as a single long, hard-coded R script with manual cross-language steps. Here we reorganize the workflow into modular functions, make the R→Python handoff explicit, and document a clear run procedure.

## What’s in this repo

- `MERFISH_analysis.r` — main end-to-end pipeline (R + Seurat)
- `find_nearest_neighbors2.py` — nearest-neighbor distance (NND) computation (Python)
- `counts_and_metadata_*.csv` — MERFISH per-sample input tables (not included here if large; see Data section)

## Key improvements vs. original release

- Modularized repetitive per-sample preprocessing into functions + loops (much shorter and easier to maintain)
- Removed unnecessary heavy exploratory steps (e.g., repeated JackStraw runs) in the default pipeline
- Increased memory limits explicitly to prevent `future.globals.maxSize` crashes in Seurat
- Automated the R→Python→R handoff for NND computation (no manual “run Python somewhere in the middle”)
- Added clear documentation for reproducibility

---

## Requirements

### R
Tested with:
- R **4.2.2**
- Seurat **5.0.3**

Required R packages:
- `Seurat`, `dplyr`, `Matrix`, `reticulate`, `scCustomize`, `ggplot2`, `sp`

Install in R:
```r
install.packages(c("dplyr","Matrix","reticulate","ggplot2","sp"))
install.packages("scCustomize")
install.packages("Seurat")

>`find_nearest_neighbors2.py` takes as input a .csv file named `all_cells.csv` that contains a list of cells, their cluster identity, and their X and Y positions. The script analyzes all cluster pairings and calculates the nearest neighbor distance for each cell in a cluster, for all pairings of clusters contained within all_cells.csv. Data is output to a folder named `python_nearest_neighbor_outputs` as .csv files for each cluster-cluster pairing.

>In Ayata_et_al_2025, `all_cells.csv` is generated from the metadata of a Seurat object using RStudio:

>**Example R code to generate "all.cells.csv"**
```R
Idents(FADPU538_cx) <- FADPU538_cx$celltype
table(Idents(FADPU538_cx))
output_dataframe=data.frame(FADPU538_cx@active.ident)
output_dataframe=cbind(output_dataframe,FADPU538_cx@meta.data$center.x)
output_dataframe=cbind(output_dataframe,FADPU538_cx@meta.data$center.y)
output_dataframe <- cbind(rownames(output_dataframe), data.frame(output_dataframe, row.names=NULL))
colnames(output_dataframe) <- c('cell', 'cluster', 'x', 'y')
write.csv(output_dataframe, 'all_cells.csv', row.names=FALSE)
```