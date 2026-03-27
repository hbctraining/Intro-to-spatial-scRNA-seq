############################################################
# Quality Control - Visium HD spatial transcriptomics
# Working with merged Seurat object
# Author: Harvard Chan Bioinformatics Core
############################################################

library(Seurat)
library(tidyverse)
library(qs)

############################################################
# Load merged Seurat object
############################################################

seurat_merged <- qread("intermediates/01_seurat_merged.qs")

############################################################
# Prepare metadata
############################################################

meta <- seurat_merged@meta.data

# View(meta)  # use interactively if desired

############################################################
# Number of cells per sample
############################################################

ggplot(meta) +
  geom_bar(aes(x = orig.ident, fill = orig.ident),
           color = "black") +
  geom_text(aes(x = orig.ident, label = after_stat(count)),
            stat = "count", vjust = -1) +
  theme_classic()

############################################################
# UMI counts (transcripts) per cell
############################################################

# Spatial overlay of nCount_Spatial.008um
SpatialFeaturePlot(
  seurat_merged,
  "nCount_Spatial.008um",
  pt.size.factor = 12,
  image.alpha = 0,
  max.cutoff = "q90"
)

# Before filtration density of nCount_Spatial.008um
ggplot(meta) +
  geom_density(aes(x = nCount_Spatial.008um, fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 30, color = "pink") +
  geom_vline(xintercept = 50, color = "lightblue") +
  scale_x_log10() +
  theme_classic()

# After filtration density of nCount_Spatial.008um
meta_filt <- subset(
  meta,
  ((orig.ident == "P5CRC") & (nCount_Spatial.008um > 30)) |
    ((orig.ident == "P5NAT") & (nCount_Spatial.008um > 50))
)

ggplot(meta_filt) +
  geom_density(aes(x = nCount_Spatial.008um, fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 30, color = "pink") +
  geom_vline(xintercept = 50, color = "lightblue") +
  scale_x_log10() +
  theme_classic()

############################################################
# Genes detected per cell
############################################################

# Spatial overlay of nFeature_Spatial.008um
SpatialFeaturePlot(
  seurat_merged,
  "nFeature_Spatial.008um",
  pt.size.factor = 12,
  image.alpha = 0,
  max.cutoff = "q90"
)

# Before filtration density of nFeature_Spatial.008um
ggplot(meta) +
  geom_density(aes(x = nFeature_Spatial.008um, fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 30, color = "pink") +
  geom_vline(xintercept = 50, color = "lightblue") +
  scale_x_log10() +
  theme_classic()

# After filtration density of nFeature_Spatial.008um
meta_filt <- subset(
  meta,
  ((orig.ident == "P5CRC") & (nFeature_Spatial.008um > 30)) |
    ((orig.ident == "P5NAT") & (nFeature_Spatial.008um > 50))
)

ggplot(meta_filt) +
  geom_density(aes(x = nFeature_Spatial.008um, fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 30, color = "pink") +
  geom_vline(xintercept = 50, color = "lightblue") +
  scale_x_log10() +
  theme_classic()

############################################################
# Complexity (novelty) score: log10GenesPerUMI
############################################################

# Add number of genes per UMI for each cell to metadata in object
seurat_merged$log10GenesPerUMI <-
  log10(seurat_merged$nFeature_Spatial.008um) /
  log10(seurat_merged$nCount_Spatial.008um)
seurat_merged@meta.data <- seurat_merged@meta.data %>%
  mutate(log10GenesPerUMI = if_else(is.na(log10GenesPerUMI), 0, log10GenesPerUMI))


# Spatial overlay of complexity score
SpatialFeaturePlot(
  seurat_merged,
  "log10GenesPerUMI",
  pt.size.factor = 14,
  image.alpha = 0,
  min.cutoff = "q10"
)

# Before filtration density of complexity score
meta <- seurat_merged@meta.data
ggplot(meta) +
  geom_density(aes(x = log10GenesPerUMI, fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 0.80) +
  theme_classic()

# After filtration density of complexity score
meta_filt <- subset(meta, log10GenesPerUMI > 0.80)
ggplot(meta_filt) +
  geom_density(aes(x = log10GenesPerUMI, fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 0.80) +
  theme_classic()

############################################################
# Mitochondrial counts ratio
############################################################

# Compute percent mito ratio
seurat_merged$mitoRatio <-
  PercentageFeatureSet(object = seurat_merged, pattern = "^MT-")
seurat_merged$mitoRatio <- seurat_merged@meta.data$mitoRatio / 100
seurat_merged@meta.data <- seurat_merged@meta.data %>%
  mutate(mitoRatio = if_else(is.na(mitoRatio), 0, mitoRatio))

# Spatial overlay of mitochondrial ratio
SpatialFeaturePlot(
  seurat_merged,
  "mitoRatio",
  pt.size.factor = 12,
  image.alpha = 0,
  max.cutoff = "q90"
)

# Before filtration density of mitochondrial ratio
meta <- seurat_merged@meta.data
ggplot(meta) +
  geom_density(aes(x = mitoRatio, fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 0.25) +
  theme_classic()

# After filtration density of mitochondrial ratio
meta_filt <- subset(meta, mitoRatio < 0.25)
ggplot(meta_filt) +
  geom_density(aes(x = mitoRatio, fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 0.25) +
  theme_classic()

############################################################
# Example structure for Exercise 1 (left here as a template)
############################################################

# seurat_merged@meta.data %>%
#   arrange(mitoRatio) %>%
#   ggplot() +
#   geom_point(aes(x = ?, 
#                  y = ?,
#                  color = ?),
#              size = 0.5) +
#   ylim(0, 3500) + xlim(0, 3500) +
#   theme_bw()

############################################################
# Filter cells (very minimal filtering)
############################################################

# Filter on nCount_Spatial.008um
seurat_filtered <- subset(
  seurat_merged,
  ((orig.ident == "P5CRC") & (nCount_Spatial.008um > 30)) |
    ((orig.ident == "P5NAT") & (nCount_Spatial.008um > 50))
)

# Filter on nFeature_Spatial.008um
seurat_filtered <- subset(
  seurat_filtered,
  ((orig.ident == "P5CRC") & (nFeature_Spatial.008um > 30)) |
    ((orig.ident == "P5NAT") & (nFeature_Spatial.008um > 50))
)

# Filter on mitoRatio
seurat_filtered <- subset(seurat_filtered, mitoRatio < 0.25)

# Filter on log10GenesPerUMI
seurat_filtered <- subset(seurat_filtered, log10GenesPerUMI > 0.80)

seurat_filtered

############################################################
# How many cells were removed?
############################################################

n_removed <- ncol(seurat_merged) - ncol(seurat_filtered)
n_removed

############################################################
# Visualizing counts data after filtering
############################################################

# Violin plots of nCounts and nFeatures after filtration
p_ncount <- VlnPlot(
  seurat_filtered,
  features = "nCount_Spatial.008um",
  pt.size = 0,
  group.by = "orig.ident"
) +
  NoLegend()

p_nfeats <- VlnPlot(
  seurat_filtered,
  features = "nFeature_Spatial.008um",
  pt.size = 0,
  group.by = "orig.ident"
) +
  NoLegend()

p_ncount | p_nfeats

# Spatial overlay of nFeature and nCount after filtration
SpatialFeaturePlot(
  seurat_filtered,
  c("nFeature_Spatial.008um", "nCount_Spatial.008um"),
  pt.size.factor = 13,
  image.alpha = 0
)

############################################################
# Save filtered Seurat object
############################################################

# Save Seurat object (uncomment to run)
qsave(seurat_filtered, "intermediates/02_seurat_filtered.qs")