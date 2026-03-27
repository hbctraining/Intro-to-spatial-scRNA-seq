############################################################
# Integration exploration (no actual integration performed)
# Visium HD spatial transcriptomics workshop
# Author: Harvard Chan Bioinformatics Core
############################################################

library(Seurat)
library(tidyverse)
library(qs)

############################################################
# Load processed Seurat object
############################################################

seurat_processed <- qread("data/04_seurat_processed.qs")
seurat_processed

############################################################
# Assess batch / condition effects: PCA vs UMAP
############################################################

# Projected PCA and UMAP, colored by sample
p1 <- DimPlot(
  seurat_processed,
  group.by  = "orig.ident",
  reduction = "full.pca.sketch"
) +
  NoLegend() +
  ggtitle("PCA (Projected)")

p2 <- DimPlot(
  seurat_processed,
  group.by  = "orig.ident",
  reduction = "full.umap.sketch"
) +
  ggtitle("UMAP (Projected)")

p1 | p2

############################################################
# Alternative UMAP visualization with split.by
############################################################

DimPlot(
  seurat_processed,
  group.by = "orig.ident",
  split.by = "orig.ident"
)

############################################################
# Custom UMAP visualization: gray background + colored foreground
############################################################

df <- FetchData(
  seurat_processed,
  c("fullumapsketch_1", "fullumapsketch_2", "orig.ident")
)

p <- ggplot(df) +
  geom_point(
    aes(x = fullumapsketch_1, y = fullumapsketch_2),
    color = "lightgray",
    alpha = 0.5,
    size = 0.1
  ) +
  geom_point(
    aes(
      x     = fullumapsketch_1,
      y     = fullumapsketch_2,
      color = orig.ident
    ),
    size = 0.1
  ) +
  theme_void() +
  facet_wrap(~orig.ident)

p

############################################################
# Batch distribution in clusters
############################################################

# Barplot of proportion of cells in each projected cluster by sample
ggplot(seurat_processed@meta.data) +
  geom_bar(
    aes(x = seurat_cluster.projected, fill = orig.ident),
    position = position_fill()
  ) +
  theme_classic()

# UMAP with projected clusters labeled
p <- DimPlot(
  seurat_processed,
  reduction = "full.umap.sketch",
  label = TRUE
) +
  ggtitle("Projected clustering")

LabelClusters(
  p,
  id        = "ident",
  fontface  = "bold",
  size      = 5,
  bg.colour = "white",
  bg.r      = 0.2,
  force     = 0
)

############################################################
# Assessing tumor-related marker genes (CEACAM5 / CEACAM6)
############################################################

# CEACAM5
VlnPlot(seurat_processed, "CEACAM5", pt.size = 0) +
  NoLegend()

FeaturePlot(seurat_processed, "CEACAM5")

# CEACAM6
VlnPlot(seurat_processed, "CEACAM6", pt.size = 0) +
  NoLegend()

FeaturePlot(seurat_processed, "CEACAM6")

# DotPlot for both markers across clusters
DotPlot(seurat_processed, c("CEACAM5", "CEACAM6"))