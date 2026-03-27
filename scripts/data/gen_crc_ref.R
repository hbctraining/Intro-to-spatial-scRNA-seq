library(Seurat)
library(tidyverse)
library(qs)

counts <- Read10X_h5("flex_reference/AggrOutput/HumanColonCancer_Flex_Multiplex_count_filtered_feature_bc_matrix.h5")
meta <- read.csv("flex_reference/SingleCell_MetaData.csv")
rownames(meta) <- meta$Barcode

seurat <- CreateSeuratObject(counts = counts,
                             meta.data = meta,
                             project = "CRC FLEX")

# Grab UMAP coordinates and put them in a reduction
umap_coords <- as.matrix(seurat@meta.data[, c("UMAP1", "UMAP2")])

# Create a DimReduc and store it in the object as "umap"
seurat[["umap"]] <- CreateDimReducObject(embeddings = umap_coords,
                                         key = "UMAP_",
                                         assay = DefaultAssay(seurat))


# Remove cells that did not pass QC
seurat <- subset(seurat, subset = (QCFilter == "Keep"))

DimPlot(seurat, group.by = "Level1", label = TRUE)
qsave(seurat, "crc_flex_ref.qs")