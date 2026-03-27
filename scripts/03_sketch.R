############################################################
# Normalization and Sketch Downsampling
# Visium HD spatial transcriptomics workshop
# Author: Harvard Chan Bioinformatics Core
# Created: December 2025
############################################################

# Load libraries
library(Seurat)
library(tidyverse)
library(scales)
library(qs)

############################################################
# Load filtered Seurat object
############################################################

seurat_filtered <- qread("intermediates/02_seurat_filtered.qs")
seurat_filtered

############################################################
# Log-normalization
############################################################

# Normalize the counts (standard log-normalized library size)
seurat_sketch <- NormalizeData(seurat_filtered)
seurat_sketch

############################################################
# Highly Variable Genes (HVGs)
############################################################

# Identify the most variable genes
seurat_sketch <- FindVariableFeatures(
  seurat_sketch,
  selection.method = "vst",
  nfeatures = 2000
)
seurat_sketch

# Identify the 15 most highly variable genes
ranked_variable_genes <- VariableFeatures(seurat_sketch)
top_genes <- ranked_variable_genes[1:15]
top_genes

############################################################
# Sketch Downsampling
############################################################

# We select 10,000 cells and create a new 'sketch' assay
seurat_sketch <- SketchData(
  object         = seurat_sketch,
  assay          = "Spatial.008um",
  ncells         = 10000,
  method         = "LeverageScore",
  sketched.assay = "sketch"
)

# Inspect object after sketching
seurat_sketch

############################################################
# Inspect metadata (includes leverage.score)
############################################################

# View(seurat_sketch@meta.data)  # use interactively if desired
head(seurat_sketch@meta.data, 15)

############################################################
# Visualizing leverage scores
############################################################

ggplot(seurat_sketch@meta.data) +
  geom_histogram(aes(x = leverage.score, fill = orig.ident),
                 alpha = 0.5,
                 bins = 100) +
  theme_classic() +
  scale_x_log10(labels = scales::label_number())

SpatialFeaturePlot(
  seurat_sketch,
  "leverage.score",
  pt.size.factor = 11,
  image.alpha = 0,
  max.cutoff = 2
)

############################################################
# Save Seurat object with sketch assay
############################################################

qsave(seurat_sketch, "intermediates/03_seurat_sketch.qs")