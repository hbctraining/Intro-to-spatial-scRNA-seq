############################################################
# Spatially Variable Genes
# Visium HD spatial transcriptomics workshop
# Author: Harvard Chan Bioinformatics Core
############################################################

library(Seurat)
library(tidyverse)
library(qs)
library(EnhancedVolcano)

############################################################
# Load processed Seurat object
############################################################

seurat_rctd <- qread("intermediates/07_seurat_rctd_unassigned.qs")
seurat_rctd <- subset(seurat_rctd, 
                      first_type == "unassigned",
                      invert = TRUE)


############################################################
# Spatially variable genes with Moran's I
############################################################
seurat_rctd <- ScaleData(seurat_rctd)

# seurat_rctd
# # S3 method for default
# seurat_rctd <- FindSpatiallyVariableFeatures(
#   seurat_rctd,
#   selection.method = "moransi",
#   verbose = TRUE)


# Split by FOV and run separately
fov_list <- c("P5CRC.008um", "P5NAT.008um")

svf_results <- lapply(fov_list, function(fov) {
  cells <- Cells(seurat_rctd[[fov]])  # get cells in this FOV
  sub <- subset(seurat_rctd, cells = cells)
  sub <- FindSpatiallyVariableFeatures(
    sub,
    selection.method = "moransi",
    verbose = TRUE
  )
  return(sub)
})