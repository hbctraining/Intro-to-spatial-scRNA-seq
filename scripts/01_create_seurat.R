############################################################
# Loading Spatial Data - Visium HD spatial transcriptomics
# Creating Seurat object(s) from SpaceRanger outputs
# Author: Harvard Chan Bioinformatics Core
############################################################

# Load libraries
library(fs)
library(tidyverse)
library(Seurat)
library(qs)

############################################################
# File structure (optional)
############################################################

# Print the first level in the downloaded directory
# dir_tree("data", recurse = 1)

# Look at the SpaceRanger output structure
dir_tree("data/P5CRC/", recurse = 1)

############################################################
# Create a Seurat object from SpaceRanger output
############################################################

# Example of looking up function arguments in an interactive session:
# ?Load10X_Spatial

# Create Seurat object for P5CRC with 8um and 16um bins
crc <- Load10X_Spatial(
  data.dir = "data/P5CRC/",
  bin.size = c(8, 16),
  slice    = "P5CRC"
)

# Print basic information
crc

############################################################
# Anatomy of a Seurat object
############################################################

# Assays
Assays(crc)

# Which assay is active by default?
DefaultAssay(crc)

# Change default assay to the 16um bins
DefaultAssay(crc) <- "Spatial.016um"
crc

############################################################
# Features and Cells
############################################################

# First features (genes) in object (from default assay)
Features(crc) %>% head()

# Number of genes in each assay
nrow(crc[["Spatial.008um"]])
nrow(crc[["Spatial.016um"]])

# First cells (columns)
Cells(crc) %>% head()

# Number of cells in each assay
ncol(crc[["Spatial.008um"]])
ncol(crc[["Spatial.016um"]])

############################################################
# Layers
############################################################

# List layers (count matrices) present
Layers(crc)

# Access a subset of the raw counts matrix
# Note: use "counts" as the layer name
LayerData(
  crc,
  assay = "Spatial.016um",
  layer = "counts"
)[1:5, 1:5]

############################################################
# Spatial fields
############################################################

# Get x,y coordinates for each bin
coords <- GetTissueCoordinates(crc)
head(coords)

# Spatial plot
SpatialDimPlot(crc)

############################################################
# Metadata
############################################################

# View metadata as a data.frame
head(crc@meta.data)

# Example of what some standard columns mean:
# - orig.ident   : Sample identity if known; defaults to "s"
# - nCount_*     : Number of UMIs per cell for that assay
# - nFeature_*   : Number of features (genes) detected per cell

# Update orig.ident to be sample name
crc@meta.data$orig.ident <- "P5CRC"
head(crc@meta.data)

# Access a single metadata column using dollar notation
head(crc$nCount_Spatial.008um)

############################################################
# Idents
############################################################

# Set cell identities to orig.ident
Idents(crc) <- "orig.ident"
head(Idents(crc))

############################################################
# Clean up single-sample object to save memory
############################################################

rm(crc)

############################################################
# Loading multiple samples using a for loop (16um bins)
############################################################

# We will use a bin size of 16um for the remainder of this script
samples <- c("P5CRC", "P5NAT")
list_seurat <- list()

for (sample in samples) {
  # Path to data directory
  data_dir <- paste0("data/", sample)
  
  # Create Seurat object and set orig.ident to be sample
  seurat <- Load10X_Spatial(
    data.dir = data_dir,
    bin.size = c(8),
    slice    = sample
  )
  seurat$orig.ident <- sample
  
  # Store Seurat object in our list
  list_seurat[[sample]] <- seurat
}

# Inspect list of Seurat objects
list_seurat

############################################################
# Merge datasets together
############################################################

# Merge P5CRC and P5NAT, add sample-specific prefixes
seurat_merged <- merge(
  x           = list_seurat[["P5CRC"]],
  y           = list_seurat[["P5NAT"]],
  add.cell.id = c("P5CRC", "P5NAT")
)

# Check merged object
seurat_merged

# At this stage you may see multiple layers, e.g., counts.1, counts.2.
# We want a single concatenated counts matrix across all cells.

# Join layers into a single "counts" layer
seurat_merged <- JoinLayers(seurat_merged)

# Confirm that there is a single counts layer
Layers(seurat_merged)

############################################################
# Evaluate merged object
############################################################

# Confirm total number of cells:
ncol(list_seurat[["P5CRC"]]) + ncol(list_seurat[["P5NAT"]])

# Should match:
ncol(seurat_merged)

# Inspect metadata to confirm cell name prefixes and orig.ident
head(seurat_merged@meta.data)
tail(seurat_merged@meta.data)

# Set Idents to orig.ident (sample IDs)
Idents(seurat_merged) <- "orig.ident"
head(Idents(seurat_merged))

############################################################
# Save merged Seurat object
############################################################

# Save integrated Seurat object (uncomment to run)
qsave(seurat_merged, "intermediates/01_seurat_merged.qs")