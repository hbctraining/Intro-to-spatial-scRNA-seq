############################################################
# Deconvolution with RCTD
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
Idents(seurat_ref) <- "Level2"
seurat_ref <- subset(seurat_ref, downsample = 1000)

# UMAP
DimPlot(seurat_ref, group.by = "Level2", label = TRUE) +
  ggtitle("CRC single-cell FLEX reference dataset")

# Number of cells
ggplot(seurat_ref@meta.data,
       aes(x = Level2, fill = Level2)) +
  geom_bar() +
  theme_bw() + NoLegend() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("CRC single-cell FLEX reference dataset")

############################################################
# Run sketch on NAT sample
############################################################


# P5NAT query object from sketch assay
nat <- subset(seurat_banksy,
              subset = (orig.ident == "P5NAT"))

# # Remove previous results 
# nat[["sketch"]] <- NULL
# nat@reductions[["pca.sketch"]] <- NULL
# nat@reductions[["full.pca.sketch"]] <- NULL
# nat@graphs[["sketch_nn"]] <- NULL
# nat@graphs[["sketch_snn"]] <- NULL

nat <- FindVariableFeatures(
  nat,
  assay = "Spatial.008um",
  selection.method = "vst",
  nfeatures = 2000)


# We select 10,000 cells and create a new 'sketch' assay
nat$leverage.score <- NULL
nat <- SketchData(
  object         = nat,
  assay          = "Spatial.008um",
  ncells         = 10000,
  method         = "LeverageScore",
  sketched.assay = "sketch"
)
# DefaultAssay(nat) <- "sketch"


nat <- FindVariableFeatures(
  nat,
  selection.method = "vst",
  nfeatures = 2000,
  assay = "sketch"
)

nat <- ScaleData(nat)
nat <- RunPCA(
  nat,
  assay = "sketch",
  reduction.name = "pca.sketch")

# # ADD THIS: project PCA to full dataset to build full.pca.sketch
# nat <- ProjectIntegration(
#   object = nat,
#   sketched.assay = "sketch",
#   assay = "Spatial.008um",
#   reduction = "pca.sketch",
#   reduction.name = "pca.sketch.full"
# )

## Fix column names in pca.sketch.full to match exactly
# colnames(nat@reductions[["pca.sketch.full"]]@cell.embeddings) <- paste0("pca.sketch_", 1:50)


############################################################
# Run RCTD
############################################################

# Create reference RCTD object
ref_counts <- seurat_ref[["RNA"]]$counts
ref_nUMI <- seurat_ref$nCount_RNA

# Requires celltypes to be factor
seurat_ref$Level2 <- factor(seurat_ref$Level2)
ref_cluster <- seurat_ref$Level2


# create the RCTD reference object
ref_rctd <- Reference(ref_counts,
                      ref_cluster,
                      ref_nUMI)


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
                    max_cores = 6)
# 11:33
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet") # this command takes ~15 mins to run
RCTD_res <- RCTD@results$results_df
colnames(RCTD_res) <- paste0(colnames(RCTD_res), "_sketch")
write.csv(RCTD_res, "intermediates/07_RCTD_res_nat.csv")

# Merge into metadata
RCTD_res$cell <- rownames(RCTD_res)
meta <- nat@meta.data
meta <- left_join(meta, RCTD_res,
                  by = "cell")
rownames(meta) <- meta$cell
nat@meta.data <- meta

# Only fill NAs within sketch cells (RCTD low-quality drops)
# Leave non-sketch cells as NA - ProjectData will handle those
sketch_cells <- colnames(nat[["sketch"]])

nat$spot_class_sketch <- as.character(nat$spot_class_sketch)
nat$first_type_sketch  <- as.character(nat$first_type_sketch)
nat$second_type_sketch <- as.character(nat$second_type_sketch)

nat@meta.data[sketch_cells, "spot_class_sketch"][is.na(nat@meta.data[sketch_cells, "spot_class_sketch"])] <- "unassigned"
nat@meta.data[sketch_cells, "first_type_sketch"][is.na(nat@meta.data[sketch_cells, "first_type_sketch"])]   <- "unassigned"
nat@meta.data[sketch_cells, "second_type_sketch"][is.na(nat@meta.data[sketch_cells, "second_type_sketch"])] <- "unassigned"

# Project onto larger dataset
nat <- ProjectData(
  object            = nat,
  assay             = "Spatial.008um",
  full.reduction    = "pca.sketch.full",
  sketched.assay    = "sketch",
  sketched.reduction = "pca.sketch",
  dims              = 1:50,
  refdata           = list(spot_class = "spot_class_sketch",
                           first_type = "first_type_sketch",
                           second_type = "second_type_sketch"))

qsave(nat, "intermediates/07_nat_rctd.qs")
SpatialDimPlot(nat, group.by = "first_type", pt.size.factor = 7,
               image.alpha = 0)

############################################################
# Merge RCTD results
############################################################


crc_meta <- read.csv("data/DeconvolutionResults_P5CRC.csv") %>%
  mutate(barcode = paste0("P5CRC_", barcode))
colnames(crc_meta) <- c("cell", "spot_class", "first_type", "second_type")
nat_meta <- nat@meta.data %>% select(cell, spot_class, first_type, second_type) %>%
  remove_rownames()

rctd_all <- rbind(crc_meta, nat_meta)

meta <- seurat_banksy@meta.data
meta <- left_join(meta, rctd_all, by = "cell")
rownames(meta) <- meta$cell
seurat_banksy@meta.data <- meta


# Rename NA to unassigned
seurat_banksy@meta.data[, "spot_class"][is.na(seurat_banksy@meta.data[, "spot_class"])] <- "unassigned"
seurat_banksy@meta.data[, "first_type"][is.na(seurat_banksy@meta.data[, "first_type"])]   <- "unassigned"
seurat_banksy@meta.data[, "second_type"][is.na(seurat_banksy@meta.data[, "second_type"])] <- "unassigned"

SpatialDimPlot(seurat_banksy, group.by = "first_type", pt.size.factor = 7,
               image.alpha = 0)

qsave(seurat_banksy, "intermediates/07_seurat_rctd_unassigned.qs")
