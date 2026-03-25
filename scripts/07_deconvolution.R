############################################################
# Spatially Derived Clusters with BANKSY
# Visium HD spatial transcriptomics workshop
# Author: Harvard Chan Bioinformatics Core
############################################################

library(Seurat)
library(tidyverse)
library(qs)

library(spacexr)

############################################################
# Load processed Seurat object
############################################################

seurat_banksy <- qread("intermediates/06_seurat_banksy_2.qs")

############################################################
# Load processed reference
############################################################

seurat_ref <- qread("data/crc_flex_ref.qs")

# UMAP
DimPlot(seurat_ref, group.by = "Level1", label = TRUE) +
  ggtitle("CRC single-cell FLEX reference dataset")

# Number of cells
ggplot(seurat_ref@meta.data,
       aes(x = Level1, fill = Level1)) +
  geom_bar() +
  theme_bw() + NoLegend() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("CRC single-cell FLEX reference dataset")


############################################################
# Run deconvolution on one sample
############################################################

# Create reference RCTD object
ref_counts <- seurat_ref[["RNA"]]$counts
ref_nUMI <- seurat_ref$nCount_RNA

# Requires celltypes to be factor
seurat_ref$Level1 <- factor(seurat_ref$Level1)
ref_cluster <- seurat_ref$Level1


# create the RCTD reference object
ref_rctd <- Reference(ref_counts,
                      ref_cluster,
                      ref_nUMI)


# P5NAT query object from sketch assay
nat <- subset(seurat_banksy,
              subset = (orig.ident == "P5NAT"))

DefaultAssay(nat) <- "sketch"

# create the RCTD query object
query_counts <- nat[["sketch"]]$counts
query_cells <- colnames(nat[["sketch"]])
query_coords <- GetTissueCoordinates(nat)[query_cells, 1:2]

# create the RCTD query object
query <- SpatialRNA(query_coords,
                    query_counts,
                    colSums(query_counts))

RCTD <- create.RCTD(query, 
                    ref_rctd, 
                    max_cores = 4)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet") # this command takes ~15 mins to run
RCTD_res <- RCTD@results$results_df

# Merge into metadata
colnames(RCTD_res) <- paste0(colnames(RCTD_res), "_sketch")
RCTD_res$cell <- rownames(RCTD_res)
meta <- nat@meta.data
meta <- left_join(meta, RCTD_res,
                  by = "cell")
rownames(meta) <- meta$cell
nat@meta.data <- meta

# Project onto larger dataset
nat <- ProjectData(
  object            = nat,
  assay             = "Spatial.008um",
  full.reduction    = "full.pca.sketch",
  sketched.assay    = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model        = "umap.sketch",
  dims              = 1:50,
  refdata           = list(spot_class = "spot_class_sketch",
                           first_type = "first_type_sketch",
                           second_type = "second_type_sketch"))


