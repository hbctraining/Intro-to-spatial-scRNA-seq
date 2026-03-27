############################################################
# Spatially Derived Clusters with BANKSY
# Visium HD spatial transcriptomics workshop
# Author: Harvard Chan Bioinformatics Core
############################################################

library(Seurat)
# remotes::install_github("prabhakarlab/Banksy")
library(Banksy)
library(SeuratWrappers)
library(tidyverse)
library(qs)

############################################################
# Load processed Seurat object
############################################################

seurat_processed <- qread("intermediates/04_seurat_processed.qs")
seurat_processed

############################################################
# Prepare coordinates for BANKSY (x, y per bin)
############################################################

seurat_processed$cell <- rownames(seurat_processed@meta.data)

# CRC and NAT sample coordinates
coords_crc <- GetTissueCoordinates(seurat_processed, image = "P5CRC.008um")
coords_nat <- GetTissueCoordinates(seurat_processed, image = "P5NAT.008um")
coords <- rbind(coords_crc, coords_nat)

# Join coords into metadata
seurat_processed@meta.data <- left_join(
  seurat_processed@meta.data,
  coords,
  by = "cell"
)
rownames(seurat_processed@meta.data) <- seurat_processed$cell

############################################################
# Run BANKSY (creates BANKSY assay)
############################################################

# Example full run (can take a while):
seurat_banksy <- RunBanksy(
  seurat_processed,
  lambda       = 0.6,
  verbose      = TRUE,
  assay        = "Spatial.008um",
  slot         = "data",
  k_geom       = 15,
  dimx         = "x",
  dimy         = "y",
  group        = "orig.ident",
  split.scale  = TRUE
)

# For convenience, load precomputed object instead:
qsave(seurat_banksy, "intermediates/06_seurat_banksy_1.qs")
# seurat_banksy <- qread("data/intermediates/06_seurat_banksy_1.qs")

# BANKSY is now the active assay
seurat_banksy

############################################################
# Inspect BANKSY features (original + neighborhood-weighted)
############################################################

# First few features (original genes)
Features(seurat_banksy) %>% head()

# Last few features (lambda-scaled neighborhood-weighted, e.g. *.m0)
Features(seurat_banksy) %>% tail()

############################################################
# PCA + kNN + clustering on BANKSY matrix
############################################################

# Example full pipeline:
options(future.globals.maxSize = 200000000000)

seurat_banksy <- RunPCA(
  seurat_banksy,
  assay          = "BANKSY",
  reduction.name = "pca.banksy",
  features       = rownames(seurat_banksy),
  npcs           = 30
)

seurat_banksy <- FindNeighbors(
  seurat_banksy,
  reduction = "pca.banksy",
  dims      = 1:30
)

seurat_banksy <- FindClusters(
  seurat_banksy,
  cluster.name = "banksy_cluster",
  resolution   = 0.8
)

# After this, the important outputs are:
# - banksy_cluster in @meta.data
# - pca.banksy in @reductions

# Once done, we can drop the BANKSY assay to save memory:
DefaultAssay(seurat_banksy) <- "Spatial.008um"
seurat_banksy[["BANKSY"]] <- NULL

# For convenience, load post-processing object:
qsave(seurat_banksy, "intermediates/06_seurat_banksy_2.qs")
# seurat_banksy <- readRDS("data/seurat_banksy_2.qs")

############################################################
# Visualize BANKSY clusters on tissue
############################################################

SpatialDimPlot(
  seurat_banksy,
  group.by      = "banksy_cluster",
  pt.size.factor = 13,
  image.alpha    = 0
)

############################################################
# Compare BANKSY clusters vs RNA-based (projected) clusters
############################################################

# Spatial overlay of both clusterings
SpatialDimPlot(
  seurat_banksy,
  group.by       = c("banksy_cluster", "seurat_cluster.projected"),
  pt.size.factor = 13,
  image.alpha    = 0
)

# UMAP comparison of both clusterings
DimPlot(
  seurat_banksy,
  group.by = c("banksy_cluster", "seurat_cluster.projected")
)

############################################################
# Batch composition per cluster: RNA-based vs BANKSY
############################################################

# Proportion of cells in each RNA-based projected cluster by sample
ggplot(seurat_banksy@meta.data) +
  geom_bar(
    aes(x = seurat_cluster.projected, fill = orig.ident),
    position = position_fill()
  ) +
  theme_classic()

# Proportion of cells in each BANKSY cluster by sample
ggplot(seurat_banksy@meta.data) +
  geom_bar(
    aes(x = banksy_cluster, fill = orig.ident),
    position = position_fill()
  ) +
  theme_classic()
