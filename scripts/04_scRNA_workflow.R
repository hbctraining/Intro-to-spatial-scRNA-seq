############################################################
# Identifying Unwanted Variation
# Visium HD spatial transcriptomics workshop
# Author: Harvard Chan Bioinformatics Core
############################################################

library(Seurat)
library(tidyverse)
library(qs)

############################################################
# Load sketched Seurat object
############################################################

seurat_sketch <- qread("intermediates/03_seurat_sketch.qs")
seurat_sketch

############################################################
# Highly Variable Features (HVGs) on full Spatial assay
############################################################

# Identify the 15 most highly variable genes on Spatial.008um
ranked_variable_genes <- VariableFeatures(seurat_sketch, assay = "Spatial.008um")
top_genes <- ranked_variable_genes[1:15]
top_genes

############################################################
# Recalculate HVGs on sketch assay (10,000 bins)
############################################################

seurat_processed <- FindVariableFeatures(
  seurat_sketch,
  selection.method = "vst",
  nfeatures = 2000,
  assay = "sketch"
)
seurat_processed

# Identify the 15 most highly variable genes on sketch assay
ranked_variable_genes <- VariableFeatures(seurat_processed, assay = "sketch")
ranked_variable_genes[1:15]

############################################################
# Visualize Highly Variable Genes
############################################################

# Top 15 variable genes
ranked_variable_genes <- VariableFeatures(seurat_processed, assay = "sketch")
top_genes <- ranked_variable_genes[1:15]

# Plot average expression vs variance, labeling top 15 HVGs
p <- VariableFeaturePlot(seurat_processed)
LabelPoints(plot = p, points = top_genes, repel = TRUE)

############################################################
# Scale data (required before PCA)
############################################################

seurat_processed <- ScaleData(seurat_processed)
seurat_processed

############################################################
# PCA on sketch assay
############################################################

seurat_processed <- RunPCA(
  seurat_processed,
  assay = "sketch",
  reduction.name = "pca.sketch"
)
seurat_processed

# PCA plots
DimPlot(
  seurat_processed,
  group.by = "orig.ident",
  reduction = "pca.sketch"
)

DimPlot(
  seurat_processed,
  reduction = "pca.sketch",
  group.by = "orig.ident",
  split.by = "orig.ident"
)

############################################################
# Selecting PC dimensions (Elbow plot)
############################################################

ElbowPlot(
  object = seurat_processed,
  reduction = "pca.sketch",
  ndims = 50
)

# We will use 50 PCs for downstream steps

############################################################
# k-Nearest Neighbors graph (on 50 PCs)
############################################################

seurat_processed <- FindNeighbors(
  seurat_processed,
  assay = "sketch",
  reduction = "pca.sketch",
  dims = 1:50
)

############################################################
# Clustering (graph-based, resolution = 0.65)
############################################################

seurat_processed <- FindClusters(
  seurat_processed,
  cluster.name = "seurat_cluster.sketched",
  resolution = 0.65
)

# Check new cluster column in metadata
seurat_processed@meta.data %>%
  head() %>%
  relocate("seurat_cluster.sketched")

# Count clusters including NA values
table(seurat_processed$seurat_cluster.sketched, useNA = "ifany")

# Idents should now be cluster labels
Idents(seurat_processed) %>% head()

############################################################
# UMAP on sketch assay (50 PCs)
############################################################

seurat_processed <- RunUMAP(
  seurat_processed,
  reduction = "pca.sketch",
  reduction.name = "umap.sketch",
  return.model = TRUE,
  dims = 1:50
)
seurat_processed

# UMAP colored by sample
DimPlot(
  seurat_processed,
  group.by = "orig.ident",
  reduction = "umap.sketch"
)

# UMAP with cluster labels
p <- DimPlot(
  seurat_processed,
  reduction = "umap.sketch",
  label = TRUE
) +
  ggtitle("Sketched clustering")

LabelClusters(
  p,
  id = "ident",
  fontface = "bold",
  size = 5,
  bg.colour = "white",
  bg.r = 0.2,
  force = 0
)

############################################################
# Project clusters and embeddings back to full dataset
############################################################

# Allow larger objects for future-based methods
options(future.globals.maxSize = 2000000000)

seurat_processed <- ProjectData(
  object            = seurat_processed,
  assay             = "Spatial.008um",
  full.reduction    = "full.pca.sketch",
  sketched.assay    = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model        = "umap.sketch",
  dims              = 1:50,
  refdata           = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

seurat_processed

# Metadata now includes projected cluster labels for all bins
head(seurat_processed@meta.data)

############################################################
# Visualize projected clusters (full dataset)
############################################################

# Use full Spatial assay
DefaultAssay(seurat_processed) <- "Spatial.008um"

# Set identities to projected cluster assignments
Idents(seurat_processed) <- "seurat_cluster.projected"

# UMAP with projected clusters
p <- DimPlot(
  seurat_processed,
  reduction = "full.umap.sketch",
  label = TRUE
) +
  ggtitle("Projected clustering")

LabelClusters(
  p,
  id = "ident",
  fontface = "bold",
  size = 5,
  bg.colour = "white",
  bg.r = 0.2,
  force = 0
)

############################################################
# Optional: ensure cluster labels are ordered numerically
############################################################

seurat_processed$seurat_cluster.projected <-
  seurat_processed$seurat_cluster.projected %>%
  as.numeric() %>%
  as.factor()

seurat_processed$seurat_cluster.projected %>% head()

# Reset identities to ordered projected clusters
Idents(seurat_processed) <- "seurat_cluster.projected"

# UMAP with ordered cluster labels
p <- DimPlot(
  seurat_processed,
  reduction = "full.umap.sketch",
  label = TRUE
) +
  ggtitle("Projected clustering")

LabelClusters(
  p,
  id = "ident",
  fontface = "bold",
  size = 5,
  bg.colour = "white",
  bg.r = 0.2,
  force = 0
)

############################################################
# Spatial visualization of projected clusters
############################################################

SpatialDimPlot(
  seurat_processed,
  pt.size.factor = 12,
  image.alpha = 0
)

############################################################
# Save processed Seurat object
############################################################

qsave(seurat_processed, "intermediates/04_seurat_processed.qs")