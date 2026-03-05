#### Code for MERFISH Analysis #########
### Analysis run using R v4.2.2 and Seurat v5.0.3 ######

library(Seurat)
library(dplyr)
library(Matrix)
library(reticulate)
library(scCustomize)
library(ggplot2)
library(sp)

options(future.globals.maxSize = 8000 * 1024^2)  # 8 GB
options(scipen = 9999)

N_GENES <- 398
ROW_VOLUME <- 464
ROW_CENTER_X <- 465
ROW_CENTER_Y <- 466
ROW_TOTAL_TRANSCRIPTS <- 472  # overwritten with colSums of genes

# FADPU668 artifact polygon (reused in multiple sections)
FADPU668_ARTIFACT_X <- c(4432,4388,4201,4185,3993,3985,3797,3789,
                         3594,3586,2975,2971,2779,2763,2388,2364,1353,1353)
FADPU668_ARTIFACT_Y <- c(6328,6184,6186,6008,5992,5793,5789,5585,
                         5585,5178,5170,4870,4966,4774,4782,4379,4371,6603)

#' Load a MERFISH CSV and create a Seurat object with metadata.
load_merfish_sample <- function(filename, sample_name, genotype) {
  dat <- t(read.table(file = filename, sep = ",", header = TRUE, row.names = 1))

  # Compute total transcript count per cell from gene rows
  genes <- dat[1:N_GENES, ]
  dat[ROW_TOTAL_TRANSCRIPTS, ] <- colSums(genes)

  # Filter: transcripts > 40, volume > 100 um^3
  dat <- dat[, dat[ROW_TOTAL_TRANSCRIPTS, ] > 40]
  dat <- dat[, dat[ROW_VOLUME, ] > 100]

  # Extract metadata before subsetting to gene rows
  volume   <- dat[ROW_VOLUME, ]
  center.x <- dat[ROW_CENTER_X, ]
  center.y <- dat[ROW_CENTER_Y, ]
  mean.RNA <- mean(dat[ROW_TOTAL_TRANSCRIPTS, ])

  # Keep gene rows only and normalize
  dat <- dat[1:N_GENES, ]
  dat <- dat / mean.RNA
  dat <- dat / volume

  # Create Seurat object
  s <- CreateSeuratObject(counts = dat)
  s@meta.data <- cbind(s@meta.data,
                       data.frame(volume),
                       data.frame(center.x),
                       data.frame(center.y))
  s$sample   <- sample_name
  s$genotype <- genotype
  return(s)
}

run_clustering <- function(obj, n_dims, resolution, jackstraw_dims = NULL, do_umap = TRUE) {
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = N_GENES, verbose = FALSE)
  obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 50, verbose = FALSE)

  if (!is.null(jackstraw_dims)) {
    obj <- JackStraw(obj, num.replicate = 100, dims = jackstraw_dims)
    obj <- ScoreJackStraw(obj, dims = 1:jackstraw_dims)
    print(JackStrawPlot(obj, dims = 1:jackstraw_dims))
  }

  obj <- FindNeighbors(obj, dims = 1:n_dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)

  if (isTRUE(do_umap)) {
    # UMAP is for visualization; skipping it does NOT affect clustering/subsetting.
    obj <- RunUMAP(obj, dims = 1:n_dims, verbose = TRUE)
  }

  return(obj)
}

#' Annotate clusters with cell-type labels and store in a metadata column.
annotate_clusters <- function(obj, labels, col_name = "celltype") {
  Idents(obj) <- obj$seurat_clusters
  k <- length(levels(obj))

  # Make labels length match the number of clusters (Seurat levels).
  if (length(labels) != k) {
    warning(sprintf(
      "annotate_clusters(): label length (%d) != number of clusters (%d). Padding/truncating with 'UNLABELED'.",
      length(labels), k
    ))
    if (length(labels) < k) {
      labels <- c(labels, rep("UNLABELED", k - length(labels)))
    } else {
      labels <- labels[1:k]
    }
  }

  names(labels) <- levels(obj)
  obj <- RenameIdents(obj, labels)
  obj[[col_name]] <- Idents(obj)
  return(obj)
}

#' Remove FADPU668 cells located within the artifact polygon region.
remove_fadpu668_artifact <- function(obj) {
  Idents(obj) <- obj$sample
  fadpu668 <- subset(obj, idents = "FADPU668")
  others   <- subset(obj, idents = "FADPU668", invert = TRUE)

  in_artifact <- point.in.polygon(fadpu668$center.x, fadpu668$center.y,
                                  FADPU668_ARTIFACT_X, FADPU668_ARTIFACT_Y)
  fadpu668$Region <- ifelse(in_artifact > 0, "Exclude", "Include")
  fadpu668_clean <- subset(fadpu668, subset = Region == "Include")

  result <- merge(fadpu668_clean, y = others)
  result <- JoinLayers(result)
  return(result)
}

#' Subset cells within a spatial polygon.
subset_by_polygon <- function(obj, poly_x, poly_y) {
  in_poly <- point.in.polygon(obj$center.x, obj$center.y, poly_x, poly_y)
  obj$in_polygon <- in_poly
  return(subset(obj, subset = in_polygon == 1))
}

#' Run NND analysis for one sample: export microglia + plaque positions,
#' call the Python NND script, read back results, annotate plaque-associated.
run_nnd_for_sample <- function(microglia_cx, plaque_cx, sample_name) {
  # Build CSV with microglia labeled by sample name and plaques labeled "Plaques"
  mic_df <- data.frame(
    cell    = colnames(microglia_cx),
    cluster = sample_name,
    x       = microglia_cx$center.x,
    y       = microglia_cx$center.y
  )
  plq_df <- data.frame(
    cell    = colnames(plaque_cx),
    cluster = "Plaques",
    x       = plaque_cx$center.x,
    y       = plaque_cx$center.y
  )
  write.csv(rbind(mic_df, plq_df), "all_cells.csv", row.names = FALSE)

  # Run Python nearest-neighbor script
  system("python find_nearest_neighbors2.py")

  # Read NND results: each row is cell, nearest_plaque, distance
  nnd_file <- file.path("python_nearest_neighbor_outputs",
                        paste0(sample_name, "_Plaques_nearest_neighbors.csv"))
  nnd_data <- read.csv(nnd_file, header = FALSE)

  # Match distances back to microglia by cell name (safe regardless of order)
  nnd_distances <- setNames(nnd_data$V3, nnd_data$V1)
  microglia_cx$NND <- as.numeric(nnd_distances[colnames(microglia_cx)])

  # Annotate: plaque-associated (NND < 15) vs distal (NND >= 15)
  microglia_cx$JP <- microglia_cx$NND < 15
  Idents(microglia_cx) <- "JP"

  return(microglia_cx)
}
sample_info <- data.frame(
  filename = c(
    "counts_and_metadata_FADTV_411.csv",
    "counts_and_metadata_FADTV_348.csv",
    "counts_and_metadata_FADPU_538.csv",
    "counts_and_metadata_FADTV_635.csv",
    "counts_and_metadata_FADPU_521.csv",
    "counts_and_metadata_FADPU_663.csv",
    "counts_and_metadata_FADPU_662.csv",
    "counts_and_metadata_FADPU_668.csv"
  ),
  sample_name = c("FADTV411","FADTV348","FADPU538","FADTV635",
                  "FADPU521","FADPU663","FADPU662","FADPU668"),
  genotype = c("5xFAD_PUhigh","5xFAD","5xFAD","Ctrl",
               "5xFAD_PUlow","Ctrl","5xFAD","5xFAD_PUlow"),
  stringsAsFactors = FALSE
)

seurat_list <- mapply(
  load_merfish_sample,
  sample_info$filename,
  sample_info$sample_name,
  sample_info$genotype,
  SIMPLIFY = FALSE
)

all.combined <- merge(seurat_list[[1]], y = seurat_list[-1],
                      merge.data = TRUE,
                      add.cell.ids = as.character(1:8))
all.combined <- JoinLayers(all.combined)
all.combined <- subset(all.combined, nFeature_RNA > 10)

# ============================================================
# SECTION 4: Initial Clustering and Annotation (Round 0)
# ============================================================

all.combined <- run_clustering(all.combined, n_dims = 39, resolution = 1.8,
                               jackstraw_dims = NULL)

round0_labels <- c(
  "Excitatory Neurons","Inhibitory Neurons","Endothelial Cells","Oligodendrocytes","Astrocytes",
  "Excitatory Neurons","Excitatory Neurons","Excitatory Neurons","Oligodendrocytes","Excitatory Neurons",
  "Inhibitory Neurons","Oligodendrocytes","Inhibitory Neurons","Excitatory Neurons","Astrocytes",
  "Microglia","Inhibitory Neurons","OPCs","Oligodendrocyte/Neuron Hybrids","Microglia",
  "Astrocyte/Neuron Hybrids","Meningeal Cells","Excitatory Neurons","Pericytes","Oligodendrocyte/Neuron Hybrids",
  "Inhibitory Neurons","Inhibitory Neurons","Endothelial Cells","Oligodendrocytes","Oligodendrocyte/Astrocyte Hybrids",
  "Oligodendrocytes","Border-Associated Macrophages","Excitatory Neurons","Microglia/Oligodendrocyte Hybrids","Endothelial Cell/Oligodendrocyte Hybrids",
  "Excitatory Neurons","Oligodendrocytes","Meningeal Cells","Choroid Plexus","Oligodendrocytes",
  "Inhibitory Neurons","Inhibitory Neurons","Inhibitory Neurons","Oligodendrocyte/Neuron Hybrids","Ependymal Cells",
  "Astrocyte/Neuron Hybrids","Astrocyte/Neuron Hybrids","Astrocyte/Neuron Hybrids","Excitatory Neurons","Microglia/Astrocyte Hybrids",
  "Microglia/Endothelial Cell Hybrids","Astrocyte/Endothelial Cell Hybrids","FADPU521 Artifacts","Inhibitory Neurons","Oligodendrocytes",
  "Lymphoid Cells","OPCs","OPC/Endothelial Cell Hybrids","Microglia","Microglia/OPC Hybrids",
  "Ifng Plaques","OPC/Astrocyte/Neuron Hybrids","Excitatory Neurons","Inhibitory Neurons","Oligodendrocytes",
  "Inhibitory Neurons","Oligodendrocyte/Astrocyte Hybrids","Oligodendrocytes","Excitatory Neurons","Neuron/Oligodendocyte Hybrids",
  "Neuron/OPC Hybrids","Broad Hybrids","Broad Hybrids","Broad Hybrids","Excitatory Neurons",
  "Neuron/Oligodendrocyte Hybrids","Oligodendrocytes","Inhibitory Neurons"
)
all.combined <- annotate_clusters(all.combined, round0_labels, "celltype")
saveRDS(all.combined, "All_celltypes_annotated.rds")

### Round 1: remove hybrids and artifacts ----
Idents(all.combined) <- all.combined$celltype
clean1 <- subset(all.combined, idents = c(
  "Excitatory Neurons","Inhibitory Neurons","Endothelial Cells","Oligodendrocytes","Astrocytes",
  "Microglia","OPCs","Meningeal Cells","Pericytes","Border-Associated Macrophages",
  "Choroid Plexus","Ependymal Cells","Lymphoid Cells","Ifng Plaques"))

clean1 <- run_clustering(clean1, n_dims = 39, resolution = 1.8, jackstraw_dims = NULL, do_umap = FALSE)

round1_labels <- c(
  "Inhibitory Neurons","Excitatory Neurons","Oligodendrocytes","Inhibitory Neurons","Astrocytes",
  "Excitatory Neurons","Excitatory Neurons","Excitatory Neurons","Oligodendrocytes","Excitatory Neurons",
  "Excitatory Neurons","Oligodendrocytes","Microglia","Endothelial Cells","Astrocytes",
  "Inhibitory Neurons","Endothelial Cells","OPCs","Excitatory Neurons","Pericytes",
  "Excitatory Neurons","Pericyte/Endothelial/Meninges Hybrids","Microglia","Meningeal Cells","Inhibitory Neurons",
  "Inhibitory Neurons","Oligodendrocytes","Inhibitory Neurons","Excitatory Neurons","Border-Associated Macrophages",
  "Oligodendrocytes","Excitatory Neurons","Oligodendrocytes","Oligodendrocytes","Excitatory Neurons",
  "Choroid Plexus","Inhibitory Neurons","Meningeal Cells","Inhibitory Neurons","Inhibitory Neurons",
  "Ependymal Cells","Excitatory Neurons","Antigen-Presenting Cells","Inhibitory Neurons","Oligodendrocytes",
  "OPCs","T Cells","Astrocytes","Microglia","Oligodendrocytes",
  "Endothelial Cells","Ifng Plaques","Excitatory Neurons","Inhibitory Neurons","Excitatory Neurons",
  "Inhibitory Neurons","Hybrid","Hybrid","Hybrid","Hybrid",
  "Hybrid","Hybrid","Hybrid","Hybrid","Hybrid",
  "Hybrid","Excitatory Neurons","Hybrid"
)
clean1 <- annotate_clusters(clean1, round1_labels, "celltype_clean1")

### Round 2: remove hybrids ----
Idents(clean1) <- clean1$celltype_clean1
clean2 <- subset(clean1, idents = c(
  "Excitatory Neurons","Inhibitory Neurons","Endothelial Cells","Oligodendrocytes","Astrocytes",
  "Microglia","OPCs","Meningeal Cells","Pericytes","Border-Associated Macrophages",
  "Choroid Plexus","Ependymal Cells","Antigen-Presenting Cells","Ifng Plaques","T Cells"))

clean2 <- run_clustering(clean2, n_dims = 39, resolution = 1.8, jackstraw_dims = NULL, do_umap = FALSE)

round2_labels <- c(
  "Excitatory Neurons","Excitatory Neurons","Oligodendrocytes","Excitatory Neurons","Astrocytes",
  "Inhibitory Neurons","Inhibitory Neurons","Excitatory Neurons","Oligodendrocytes","Oligodendrocytes",
  "Astrocytes","Excitatory Neurons","Excitatory Neurons","Endothelial Cells","Endothelial Cells",
  "Inhibitory Neurons","Microglia","Microglia","Inhibitory Neurons","OPCs",
  "Inhibitory Neurons","Inhibitory Neurons","Meningeal Cells","Pericytes","Inhibitory Neurons",
  "Inhibitory Neurons","Oligodendrocytes","Excitatory Neurons","Border-Associated Macrophages","Oligodendrocytes",
  "Excitatory Neurons","Oligodendrocytes","Inhibitory Neurons","Choroid Plexus","Oligodendrocytes",
  "Inhibitory Neurons","Inhibitory Neurons","Ependymal Cells","Excitatory Neurons","Antigen-Presenting Cells",
  "Inhibitory Neurons","OPC/Neurons","T Cells","Astrocytes","Inhibitory Neurons",
  "Microglia","Oligodendrocytes","Ifng Plaques","Inhibitory Neurons","Excitatory Neurons",
  "Inhibitory Neurons","Hybrids","Excitatory Neurons","Excitatory Neurons","Excitatory Neurons",
  "Inhibitory Neurons","Oligodendrocytes","Excitatory Neurons","Excitatory Neurons","Oligodendrocytes",
  "Hybrids","Oligodendrocytes","Excitatory Neurons","Excitatory Neurons","Excitatory Neurons"
)
clean2 <- annotate_clusters(clean2, round2_labels, "celltype_clean2")

### Round 3: remove hybrids ----
Idents(clean2) <- clean2$celltype_clean2
clean3 <- subset(clean2, idents = c(
  "Excitatory Neurons","Inhibitory Neurons","Endothelial Cells","Oligodendrocytes","Astrocytes",
  "Microglia","OPCs","Meningeal Cells","Pericytes","Border-Associated Macrophages",
  "Choroid Plexus","Ependymal Cells","Antigen-Presenting Cells","Ifng Plaques","T Cells"))

clean3 <- run_clustering(clean3, n_dims = 36, resolution = 1.8, jackstraw_dims = NULL, do_umap = FALSE)

round3_labels <- c(
  "Inhibitory Neurons","Excitatory Neurons","Oligodendrocytes","Astrocytes","Excitatory Neurons",
  "Oligodendrocytes","Inhibitory Neurons","Excitatory Neurons","Oligodendrocytes","Excitatory Neurons",
  "Excitatory Neurons","Excitatory Neurons","Astrocytes","Microglia","Endothelial Cells",
  "Inhibitory Neurons","Endothelial Cells","Inhibitory Neurons","OPCs","Excitatory Neurons",
  "Excitatory Neurons","Microglia","Pericytes","Meningeal Fibroblasts","Inhibitory Neurons",
  "Inhibitory Neurons","Excitatory Neurons","Oligodendrocytes","Excitatory Neurons","Border-Associated Macrophages",
  "Oligodendrocytes","Antigen-Presenting Cells","Choroid Plexus Cells","Oligodendrocytes","Inhibitory Neurons",
  "Inhibitory Neurons","Inhibitory Neurons","Ependymal Cells","Excitatory Neuron - Artifacts","Inhibitory Neurons",
  "Astrocyte/Neuron Hybrids","Oligodendrocytes","T Cells","Neuron/BAM Hybrids","Ifng Plaques",
  "Hybrids","Hybrids","Hybrids","Hybrids","Inhibitory Neurons",
  "Hybrids","Hybrids"
)
clean3 <- annotate_clusters(clean3, round3_labels, "celltype_clean3")

### Round 4: remove hybrids and artifacts ----
Idents(clean3) <- clean3$celltype_clean3
clean4 <- subset(clean3, idents = c(
  "Excitatory Neurons","Inhibitory Neurons","Endothelial Cells","Oligodendrocytes","Astrocytes",
  "Microglia","OPCs","Meningeal Fibroblasts","Pericytes","Border-Associated Macrophages",
  "Choroid Plexus Cells","Ependymal Cells","Antigen-Presenting Cells","Ifng Plaques","T Cells"))

clean4 <- run_clustering(clean4, n_dims = 38, resolution = 2.4, jackstraw_dims = NULL, do_umap = FALSE)

round4_labels <- c(
  "Oligodendrocytes","Excitatory Neurons","Inhibitory Neurons","Excitatory Neurons","Excitatory Neurons",
  "Inhibitory Neurons","Excitatory Neurons","Excitatory Neurons","Oligodendrocytes","Endothelial Cells",
  "Excitatory Neurons","Inhibitory Neurons","Astrocytes","Excitatory Neurons","Endothelial Cells",
  "Astrocytes","Oligodendrocytes","Oligodendrocytes","Inhibitory Neurons","Astrocytes",
  "OPCs","Excitatory Neurons","Pericytes","Microglia","Inhibitory Neurons",
  "Microglia","Excitatory Neurons","Inhibitory Neurons","Meningeal Fibroblasts","Microglia",
  "Astrocytes","Excitatory Neurons","Inhibitory Neurons","Oligodendrocytes","Inhibitory Neurons",
  "Border-Associated Macrophages","Excitatory Neurons","Oligodendrocytes","Inhibitory Neurons","Oligodendrocytes",
  "Excitatory Neurons","Choroid Plexus Cells","Oligodendrocytes","Inhibitory Neurons","Microglia/Neuron Hybrids",
  "Inhibitory Neurons","Meningeal Fibroblasts","Ependymal Cells","Excitatory Neurons","Excitatory Neurons",
  "Oligodendrocytes","Antigen-Presenting Cells","Inhibitory Neurons","Oligodendrocytes","Excitatory Neurons",
  "T Cells","Inhibitory Neurons","Excitatory/BAMs","OL/T Cells","Oligodendrocytes",
  "Plaques","Hybrids","Hybrids","Inhibitory Neurons","64",
  "Hybrids","Hybrids","Excitatory Neurons","Hybrids","Excitatory Neurons",
  "Hybrids"
)
clean4 <- annotate_clusters(clean4, round4_labels, "celltype_clean4")

### Remove FADPU668 artifact region ----
clean4 <- remove_fadpu668_artifact(clean4)

### Round 5: final clean ----
Idents(clean4) <- clean4$celltype_clean4
clean5 <- subset(clean4, idents = c(
  "Excitatory Neurons","Inhibitory Neurons","Endothelial Cells","Oligodendrocytes","Astrocytes",
  "Microglia","OPCs","Meningeal Fibroblasts","Pericytes","Border-Associated Macrophages",
  "Choroid Plexus Cells","Ependymal Cells","Antigen-Presenting Cells","T Cells"))

clean5 <- run_clustering(clean5, n_dims = 38, resolution = 2.4, jackstraw_dims = NULL, do_umap = FALSE)

round5_labels <- c(
  "Oligodendrocytes","Inhibitory Neurons","Excitatory Neurons","Excitatory Neurons","Inhibitory Neurons",
  "Excitatory Neurons","Excitatory Neurons","Excitatory Neurons","Oligodendrocytes","Astrocytes",
  "Endothelial Cells","Inhibitory Neurons","Endothelial Cells","Microglia","Excitatory Neurons",
  "Excitatory Neurons","Astrocytes","Oligodendrocytes","Oligodendrocytes","Excitatory Neurons",
  "Excitatory Neurons","OPCs","Astrocytes","Pericytes","Excitatory Neurons",
  "Inhibitory Neurons","Inhibitory Neurons","Astrocytes","Inhibitory Neurons","Meningeal Fibroblasts",
  "Excitatory Neurons","Microglia","Inhibitory Neurons","Oligodendrocytes","Microglia",
  "Inhibitory Neurons","Border-Associated Macrophages","Excitatory Neurons","Oligodendrocytes","Oligodendrocytes",
  "Choroid Plexus Cells","Inhibitory Neurons","Newly-Formed Oligodendrocytes","Inhibitory Neurons","Meningeal Fibroblasts",
  "Inhibitory Neurons","Ependymal Cells","Excitatory Neurons","Excitatory Neurons","Inhibitory Neurons",
  "Antigen-Presenting Cells","Oligodendrocytes","Oligodendrocytes","T Cells","Oligodendrocytes",
  "Oligodendrocytes","Excitatory Neurons","Inhibitory Neurons","Excitatory Neurons","Excitatory Neurons",
  "Oligodendrocytes","Excitatory Neurons","Inhibitory Neurons","Inhibitory Neurons","Inhibitory Neurons",
  "Inhibitory Neurons","Inhibitory Neurons","Oligodendrocytes","Excitatory Neurons","Inhibitory Neurons",
  "Excitatory Neurons","Excitatory Neurons","Excitatory Neurons"
)
# FIX: original code referenced "all.combined.clean6" here, which was never created.
# It should be clean5 (the object created in the lines above).
clean5 <- annotate_clusters(clean5, round5_labels, "celltype_clean")

saveRDS(clean5, "All_celltypes_no_hybrids_annotated.rds")

Idents(all.combined) <- all.combined$celltype
microglia_unclean <- subset(all.combined, idents = c(
  "Microglia","Microglia/Oligodendrocyte Hybrids","Microglia/Endothelial Cell Hybrids",
  "Microglia/Astrocyte Hybrids","Microglia/OPC Hybrids","Border-Associated Macrophages"))

# Remove FADPU668 artifact region
microglia_unclean <- remove_fadpu668_artifact(microglia_unclean)

### Microglia subclustering - round 1 ----
microglia_unclean <- run_clustering(microglia_unclean, n_dims = 29, resolution = 2.8,
                                    jackstraw_dims = NULL, do_umap = FALSE)

# Keep only microglia clusters (remove BAM/hybrid clusters)
microglia_unclean2 <- subset(microglia_unclean,
  idents = c("0","1","2","3","4","7","9","10","12","14","15","19","23"))

### Microglia subclustering - round 2 ----
microglia_unclean2 <- run_clustering(microglia_unclean2, n_dims = 18, resolution = 2.8,
                                     jackstraw_dims = NULL, do_umap = FALSE)

# Keep only microglia clusters
microglia_clean <- subset(microglia_unclean2,
  idents = c("0","1","2","3","4","5","6","7","8","9","10","11","12","14",
             "16","17","19","20","21","22","23","24","25","26","27","28","29","30"))

### Microglia subclustering - round 3 (final) ----
microglia_clean <- run_clustering(microglia_clean, n_dims = 17, resolution = 0.7,
                                  jackstraw_dims = NULL)

saveRDS(microglia_clean, "Microglia_all_samples_no_hybrids.rds")

Idents(all.combined) <- all.combined$celltype
plaque_cells <- subset(all.combined, idents = "Ifng Plaques")

# Cortex polygon definitions per sample
cortex_polygons <- list(
  FADPU538 = list(
    x = c(3573,3893,3770,4161,4493,4814,5080,5195,5219,5174,5084,5039,5035,5068,
          5019,4959,4757,4611,4439,4016,3700,3758,3936,4065,4446,6898,7171,4861,
          3609,3230,3089,2980,3054,3437),
    y = c(6953,7315,7350,7694,7239,6875,6416,5952,5737,5669,4865,4684,4037,3114,
          2344,2052,1744,1599,1593,1667,1659,1138,870,563,183,152,6813,9488,
          8418,7844,7618,7275,7022,7020)
  ),
  FADTV348 = list(
    x = c(1802,1866,2059,2254,2370,2452,2415,2248,1969,1077,266,-81,919,1247,
          2372,3139,3193,3195,3271,3391,3435,2720,2535,2237,2048,1833,1548,1417,
          1368,1367,1455,1755),
    y = c(4103,4382,4567,4659,4719,4702,4121,3555,2987,3041,4301,7507,10339,11143,
          11613,11296,11104,10952,10719,10394,9823,9719,9635,9060,8425,7613,6749,
          6320,5915,5485,5098,4069)
  )
)

### FADPU538 NND analysis ----
FADPU538_mic <- subset(microglia_clean, subset = sample == "FADPU538")
FADPU538_plq <- subset(plaque_cells, subset = sample == "FADPU538")
FADPU538_mic_cx <- subset_by_polygon(FADPU538_mic, cortex_polygons$FADPU538$x, cortex_polygons$FADPU538$y)
FADPU538_plq_cx <- subset_by_polygon(FADPU538_plq, cortex_polygons$FADPU538$x, cortex_polygons$FADPU538$y)
FADPU538_mic_cx$region <- "Cx"
FADPU538_cx <- run_nnd_for_sample(FADPU538_mic_cx, FADPU538_plq_cx, "FADPU538")

### FADTV348 NND analysis ----
FADTV348_mic <- subset(microglia_clean, subset = sample == "FADTV348")
FADTV348_plq <- subset(plaque_cells, subset = sample == "FADTV348")
FADTV348_mic_cx <- subset_by_polygon(FADTV348_mic, cortex_polygons$FADTV348$x, cortex_polygons$FADTV348$y)
FADTV348_plq_cx <- subset_by_polygon(FADTV348_plq, cortex_polygons$FADTV348$x, cortex_polygons$FADTV348$y)
FADTV348_mic_cx$region <- "Cx"
FADTV348_cx <- run_nnd_for_sample(FADTV348_mic_cx, FADTV348_plq_cx, "FADTV348")

### merge cortical microglia from both 5xFAD samples ----
Microglia_cx_5xFAD <- merge(FADTV348_cx, FADPU538_cx)
Microglia_cx_5xFAD <- JoinLayers(Microglia_cx_5xFAD)

##Data Output
saveRDS(Microglia_cx_5xFAD, "Microglia_5xFAD_cortex_with_plaque_NND.rds")
