############################################################
# SpatialCellChat
# Visium HD spatial transcriptomics workshop
# Author: Harvard Chan Bioinformatics Core
############################################################

library(Seurat)
library(tidyverse)
library(qs)
library(future)
library(jsonlite)

# devtools::install_github("jinworks/SpatialCellChat")
# BiocManager::install("BiocNeighbors")
# devtools::install_github("JEFworks-Lab/MERINGUE")
# devtools::install_github("KlugerLab/ALRA")
# devtools::install_github("zdebruine/RcppML")
# devtools::install_github("jinworks/SpatialCellChat")

# library(SpatialCellChat)
library(CellChat)

############################################################
# Load processed Seurat object
############################################################

seurat_rctd <- qread("intermediates/07_seurat_rctd_unassigned.qs")
seurat_rctd <- subset(seurat_rctd,
                      first_type %in% c("CD4 T cell",
                                        "CD8 T cell",
                                        "Macrophage",
                                        "Endothelial"))
Idents(seurat_rctd) <- "first_type"
seurat_rctd$samples <- factor(seurat_rctd$orig.ident)
seurat_rctd
# An object of class Seurat 
# 36170 features across 115578 samples within 2 assays 
# Active assay: Spatial.008um (18085 features, 2000 variable features)
# 2 layers present: counts, data
# 1 other assay present: sketch
# 5 dimensional reductions calculated: pca.sketch, umap.sketch, full.pca.sketch, full.umap.sketch, pca.banksy
# 2 spatial fields of view present: P5CRC.008um P5NAT.008um

SpatialDimPlot(seurat_rctd,
               image.alpha = 0.1,
               pt.size.factor = 12)

############################################################
# Cell chat (multiple sample)
############################################################

# Need normalized counts

# =============================================================================
# CellChat: Multi-Sample Spatial Transcriptomics Analysis
# =============================================================================
# Adapted from:
# https://github.com/jinworks/CellChat/blob/master/tutorial/
#   CellChat_analysis_of_multiple_spatial_transcriptomics_datasets.html
#
# Input object: seurat_rctd
#   - Two samples merged into one Seurat object
#   - Normalized counts: seurat_rctd[["Spatial.008um"]]@data
#   - seurat_rctd@meta.data$samples   -> sample labels (e.g. "Sample1", "Sample2")
#   - seurat_rctd@meta.data$first_type -> cell type labels from RCTD deconvolution
# =============================================================================

# ── 1. Extract inputs from the Seurat object ──────────────────────────────────

# Normalized gene-by-cell expression matrix from the 8 µm Spatial assay
data_input <- GetAssayData(seurat_rctd,
                           assay = "Spatial.008um",
                           layer = "data")

# Cell metadata: keep the two columns CellChat needs
#   'labels'  -> cell type per spot (from RCTD first_type)
#   'samples' -> which sample each spot belongs to
meta <- data.frame(
  labels  = seurat_rctd@meta.data$first_type,
  samples = seurat_rctd@meta.data$samples,
  row.names = colnames(seurat_rctd)
)

# Spatial coordinates – fetched separately per sample image then row-bound.
# Replace "Sample1" / "Sample2" with the exact image names stored in your
# Seurat object (check with: Images(seurat_rctd))
spatial_locs_crc <- GetTissueCoordinates(
  seurat_rctd,
  image = "P5CRC.008um",     # <-- replace with your first image name
  scale = NULL,          # use full-resolution (not scaled) coordinates
  cols  = c("imagerow", "imagecol")
) %>% dplyr::select(!c(cell))

spatial_locs_nat <- GetTissueCoordinates(
  seurat_rctd,
  image = "P5NAT.008um",     # <-- replace with your second image name
  scale = NULL,
  cols  = c("imagerow", "imagecol")
) %>% dplyr::select(!c(cell))

# Merge both coordinate data frames; row order must match colnames(seurat_rctd)
spatial_locs <- rbind(spatial_locs_crc, spatial_locs_nat)

# ── 1.1 Spatial scaling factors ────────────────────────────────────────────────
# CellChat v2 requires spatial.factors as a data.frame with two columns:
#
#   ratio : converts pixels → µm  =  spot_size_um / spot_diameter_fullres
#   tol   : the tolerance radius in µm used to define "neighbouring" spots
#           = spot_size_um / 2  (i.e. half the bin/spot diameter)
#
# For Visium HD 8 µm bins:
#   spot_size_um            = 8  (the bin size in microns)
#   spot_diameter_fullres   = from scalefactors_json.json ("spot_diameter_fullres")

sf_json     <- jsonlite::fromJSON("data/P5CRC/binned_outputs/square_008um/spatial/scalefactors_json.json")
spot_size_um <- sf_json$bin_size_um           # 8 for Visium HD 8 µm bins

conversion_factor <- spot_size_um / sf_json$spot_diameter_fullres

spatial_factors <- data.frame(
  ratio = conversion_factor,   # pixels-to-µm conversion factor
  tol   = spot_size_um / 2     # neighbour tolerance radius in µm
)

# ── 2. Create a single CellChat object (both samples together) ─────────────────
# Passing the 'samples' column lets CellChat aggregate communication across
# replicates while still using spatial coordinates for distance constraints.

cellchat <- createCellChat(
  object          = data_input,
  meta            = meta,
  group.by        = "labels",     # column that defines cell types
  datatype        = "spatial",
  coordinates     = spatial_locs,
  spatial.factors = spatial_factors
)

rm(data_input)

# Inspect the object
print(cellchat)


# ── 3. Set the ligand-receptor database ───────────────────────────────────────
# Use CellChatDB.human for human data, CellChatDB.mouse for mouse data.

cellchat_db <- CellChatDB.human   # <-- change to CellChatDB.mouse if needed

# Subset to secreted signaling only (fastest; omits Cell-Cell Contact & ECM).
# To use the full database instead, comment the next line and uncomment the one
# after it.
# cellchat_db_use <- subsetDB(cellchat_db, search = "Secreted Signaling")
cellchat_db_use <- subsetDB(cellchat_db, search = "Cell-Cell Contact")
# cellchat_db_use <- cellchat_db   # use the full database

cellchat@DB <- cellchat_db_use


# ── 4. Preprocess gene expression data ────────────────────────────────────────

# Subset to signaling genes only – reduces memory and speeds up computation.
# This step is required even when using the full database.
cellchat <- subsetData(cellchat)

# Run in parallel to speed up overexpression testing
plan("multisession", workers = 6)   # adjust 'workers' to available CPU cores

# Identify genes that are over-expressed relative to other cell types
cellchat <- identifyOverExpressedGenes(cellchat)

# Increases the size of the default vector
options(future.globals.maxSize= 200000000000)

# Flag ligand-receptor pairs where at least one partner is over-expressed
cellchat <- identifyOverExpressedInteractions(cellchat)

# (Optional) Smooth expression over the protein-protein interaction network to
# reduce dropout noise. If you run this, set raw.use = FALSE below.
# cellchat <- projectData(cellchat, PPI.human)
cellchat <- smoothData(cellchat, adj = PPI.human)



# ── 5. Infer cell-cell communication probabilities ────────────────────────────
# distance.use = TRUE  : spatial distances constrain which cell pairs can talk
# interaction.range    : max diffusion distance in µm (literature: ~250 µm)
# scale.distance       : scaling factor; tune if you get warnings about distances
# contact.dependent    : TRUE adds juxtacrine / cell-contact signalling
# contact.range        : max distance (µm) for contact-dependent signals;
#                        ~8-10 µm = cell diameter for single-cell res data,
#                        ~100 µm  = centre-to-centre distance for 10x Visium

cellchat <- computeCommunProb(
  cellchat,
  distance.use       = TRUE,
  interaction.range = 200,
  scale.distance     = 0.2,
  contact.dependent  = TRUE,
  contact.range      = 10,
  raw.use            = TRUE
)

# Remove any interactions supported by very few cells (avoids spurious signals)
cellchat <- filterCommunication(cellchat, min.cells = 10)


# ── 6. Aggregate to pathway level ─────────────────────────────────────────────
# Summarise L-R pair probabilities into signalling pathway probabilities
cellchat <- computeCommunProbPathway(cellchat)

# Count interactions and sum weights into a flat network matrix
cellchat <- aggregateNet(cellchat)

# Quick sanity-check plot: circle plots of total interactions and total strength
group_size <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(
  cellchat@net$count,
  vertex.weight = rowSums(cellchat@net$count),
  weight.scale  = TRUE,
  label.edge    = FALSE,
  title.name    = "Number of interactions"
)
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = rowSums(cellchat@net$weight),
  weight.scale  = TRUE,
  label.edge    = FALSE,
  title.name    = "Interaction strength"
)


# ── 7. Network centrality analysis ────────────────────────────────────────────
# Compute how much each cell type sends vs receives across all pathways
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Scatter plot: outgoing (x) vs incoming (y) signal strength per cell type
netAnalysis_signalingRole_scatter(cellchat)


# ── 8. Visualise a specific signalling pathway ────────────────────────────────
# List all pathways with significant communication:
print(cellchat@netP$pathways)

# Pick a pathway of interest (edit as needed)
pathways_show <- "PECAM1"   # <-- replace with a pathway from the list above

# -- Circle plot (aggregated across both samples) --
par(mfrow = c(1, 1))
netVisual_aggregate(cellchat, signaling = pathways_show, layout = "circle")

# -- Spatial plot for each sample separately --
# 'sample.use' selects which sample's coordinates to overlay
sample_names <- unique(meta$samples)   # e.g. c("Sample1", "Sample2")

for (samp in sample_names) {
  # par(mfrow = c(1, 1))
  netVisual_aggregate(
    cellchat,
    signaling        = pathways_show,
    sample.use       = samp,
    layout           = "spatial",
    edge.width.max   = 2,
    vertex.size.max  = 1,
    alpha.image      = 0.2,
    vertex.label.cex = 0,
    signaling.name       = paste(pathways_show, "-", samp)
  )
}

# Heatmap of signalling roles for this pathway
netAnalysis_signalingRole_network(
  cellchat,
  signaling = pathways_show,
  width     = 8,
  height    = 2.5,
  font.size = 10
)


# ── 9. Split object by sample for comparative analysis ────────────────────────
# Subset the unified CellChat object into one object per sample, then merge
# them as a named list so all comparison functions work correctly.

cellchat_list <- lapply(sample_names, function(s) {
  cells_use  <- rownames(meta)[meta$samples == s]
  
  # Get only the cell types that actually exist in this sample
  idents_use <- unique(meta$labels[meta$samples == s])
  
  subsetCellChat(cellchat, idents.use = idents_use, cells.use = cells_use)
})
names(cellchat_list) <- sample_names

# ── 10. Merge for comparative analysis ────────────────────────────────────────
# mergeCellChat combines the per-sample objects so comparison functions work
cellchat_merged <- mergeCellChat(
  cellchat_list,
  cell.prefix = TRUE,
  add.names = names(cellchat_list)
)


# ── 11. Compare overall communication between samples ─────────────────────────

# Bar plots: total interaction count and total strength per sample
gg1 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2))
gg2 <- compareInteractions(cellchat_merged, show.legend = FALSE,
                           group = c(1, 2), measure = "weight")
gg1 + gg2

# Differential circle plot: red = increased, blue = decreased in sample 2
par(mfrow = c(1, 1))
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE)

# Differential heatmap
netVisual_heatmap(cellchat_merged)


# ── 12. Pathway-level comparison ──────────────────────────────────────────────

# Rank pathways by change in information flow between samples
gg3 <- rankNet(cellchat_merged, mode = "comparison", measure = "weight",
               stacked = TRUE,  do.stat = TRUE)
gg4 <- rankNet(cellchat_merged, mode = "comparison", measure = "weight",
               stacked = FALSE, do.stat = TRUE)
gg3 + gg4

# Outgoing and incoming signalling role heatmaps per sample
pathway_union <- union(
  cellchat_list[[1]]@netP$pathways,
  cellchat_list[[2]]@netP$pathways
)

ht_out <- lapply(seq_along(cellchat_list), function(i) {
  netAnalysis_signalingRole_heatmap(
    cellchat_list[[i]],
    pattern   = "outgoing",
    signaling = pathway_union,
    title     = names(cellchat_list)[i],
    width     = 5,
    height    = 8
  )
})
ht_out[[1]] + ht_out[[2]]   # side-by-side outgoing heatmaps

ht_in <- lapply(seq_along(cellchat_list), function(i) {
  netAnalysis_signalingRole_heatmap(
    cellchat_list[[i]],
    pattern   = "incoming",
    signaling = pathway_union,
    title     = names(cellchat_list)[i],
    width     = 5,
    height    = 8
  )
})
ht_in[[1]] + ht_in[[2]]    # side-by-side incoming heatmaps


# ── 13. Network similarity (functional & structural) ─────────────────────────
# Embed signalling pathways in 2-D to see which pathways are shared or unique

# Functional similarity (similar senders/receivers)
cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "functional")
cellchat_merged <- netEmbedding(cellchat_merged,                 type = "functional")
cellchat_merged <- netClustering(cellchat_merged,                type = "functional")
netVisual_embeddingPairwise(cellchat_merged, type = "functional", label.size = 3.5)

# Rank pathways by their distance between the two samples (most changed first)
rankSimilarity(cellchat_merged, type = "functional")

# Structural similarity (similar network topology)
cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "structural")
cellchat_merged <- netEmbedding(cellchat_merged,                 type = "structural")
cellchat_merged <- netClustering(cellchat_merged,                type = "structural")
netVisual_embeddingPairwise(cellchat_merged, type = "structural", label.size = 3.5)
rankSimilarity(cellchat_merged, type = "structural")


# ── 14. L-R pair bubble plots ─────────────────────────────────────────────────
# Identify specific ligand-receptor pairs that are up- or down-regulated
# between the two samples.
# Edit sources.use / targets.use to focus on cell types of interest.

# Interactions INCREASED in sample 2
gg5 <- netVisual_bubble(
  cellchat_merged,
  sources.use    = NULL,   # NULL = all sources; or e.g. c(1, 2)
  targets.use    = NULL,
  comparison     = c(1, 2),
  max.dataset    = 2,      # interactions dominant in dataset 2
  title.name     = paste0("Increased in ", sample_names[2]),
  angle.x        = 45,
  remove.isolate = TRUE
)

# Interactions DECREASED in sample 2
gg6 <- netVisual_bubble(
  cellchat_merged,
  sources.use    = NULL,
  targets.use    = NULL,
  comparison     = c(1, 2),
  max.dataset    = 1,      # interactions dominant in dataset 1
  title.name     = paste0("Decreased in ", sample_names[2]),
  angle.x        = 45,
  remove.isolate = TRUE
)

gg5 + gg6


# ── 15. DEG-based signalling analysis ─────────────────────────────────────────
# Map differentially expressed genes onto L-R pairs to find mechanistic drivers

pos_dataset   <- sample_names[2]                   # define "positive" (e.g. treated) sample
features_name <- paste0(pos_dataset, "_merged")

# Compute pseudo-bulk DE between samples within the merged object
cellchat_merged <- identifyOverExpressedGenes(
  cellchat_merged,
  group.dataset = "datasets",
  pos.dataset   = pos_dataset,
  features.name = features_name,
  only.pos      = FALSE,
  thresh.pc     = 0.1,
  thresh.fc     = 0.1
)

# Map DE results back to cell-cell communication pairs
net <- netMappingDEG(cellchat_merged, features.name = features_name)

# Up-regulated L-R pairs in sample 2 (ligand log2FC > 0.05)
net_up <- subsetCommunication(
  cellchat_merged,
  net            = net,
  datasets       = pos_dataset,
  ligand.logFC   = 0.05,
  receptor.logFC = NULL
)

# Down-regulated L-R pairs in sample 2 (ligand log2FC < -0.05)
net_down <- subsetCommunication(
  cellchat_merged,
  net            = net,
  datasets       = sample_names[1],
  ligand.logFC   = -0.05,
  receptor.logFC = NULL
)

# Extract gene sets to drive chord diagrams
gene_up   <- extractGeneSubsetFromPair(net_up,   cellchat_merged)
gene_down <- extractGeneSubsetFromPair(net_down, cellchat_merged)

# Chord diagrams: up- and down-regulated L-R pairs side by side
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_chord_gene(
  cellchat_list[[2]],
  sources.use = NULL,
  targets.use = NULL,
  slot.name   = "net",
  net         = net_up,
  lab.cex     = 0.8,
  small.gap   = 3.5,
  title.name  = paste0("Up-regulated in ", sample_names[2])
)
netVisual_chord_gene(
  cellchat_list[[1]],
  sources.use = NULL,
  targets.use = NULL,
  slot.name   = "net",
  net         = net_down,
  lab.cex     = 0.8,
  small.gap   = 3.5,
  title.name  = paste0("Down-regulated in ", sample_names[2])
)


# ── 16. Spatial feature plots ─────────────────────────────────────────────────
# Visualise the expression of specific genes on the tissue slide per sample

genes_to_plot <- c("CXCL12", "CXCR4")   # <-- replace with genes of interest

for (samp in sample_names) {
  spatialFeaturePlot(
    cellchat,
    features      = genes_to_plot,
    sample.use    = samp,
    point.size    = 0.8,
    color.heatmap = "Reds",
    direction     = 1
  )
}


# ── 17. Save outputs ──────────────────────────────────────────────────────────

# Unified CellChat object (both samples together)
saveRDS(cellchat,        file = "cellchat_spatial_merged.rds")

# Per-sample list and the pairwise comparison object
saveRDS(cellchat_list,   file = "cellchat_spatial_list.rds")
saveRDS(cellchat_merged, file = "cellchat_spatial_comparison.rds")

message("Done. CellChat objects saved.")







# ###########################################################
# Idents(seurat_rctd) <- "first_type"
# # CellChat expects a column called "samples" that is a factor
# seurat_rctd$samples <- factor(seurat_rctd$orig.ident)
# 
# cellchat <- createCellChat(object = seurat_rctd, 
#                            group.by = "ident", 
#                            assay = "Spatial.008um")
# 
# # Load database
# CellChatDB <- CellChatDB.human 
# CellChatDB.use <- subsetDB(CellChatDB)
# cellchat@DB <- CellChatDB.use
# 
# showDatabaseCategory(CellChatDB)
# 
# 
# # Increases the size of the default vector
# options(future.globals.maxSize= 200000000000)
# 
# cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
# cellchat <- updateCellChat(cellchat)
# future::plan("multisession", workers = 8) # recommend running with at 8-16 cores
# 
# 
# cellchat <- identifyOverExpressedGenes(cellchat) # may take a couple minutes
# cellchat <- identifyOverExpressedInteractions(cellchat) # may take a couple minutes
# 
# 
# # this next command takes 0.5-2+ hours
# # can choose various methods for calculating average gene exp per group, 
# # 'triMean' allegedly produces fewer but stronger interactions
# cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
#                               distance.use = FALSE, interaction.range = 250, 
#                               scale.distance = NULL,
#                               contact.dependent = TRUE, contact.range = 100)
# 
# 
# #  filter out the cell-cell communication if < 10 cells per group
# cellchat <- filterCommunication(cellchat, min.cells = 10)
# 
# df.net <- subsetCommunication(cellchat) %>% dplyr::arrange(pval) 
# df.net %>% View()
# 
# groupSize <- as.numeric(table(cellchat@idents))
# par(mfrow = c(1,2), xpd=TRUE)
# netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), 
#                  weight.scale = T, label.edge= F, title.name = "Number of interactions")
# netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), 
#                  weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# 
# 
# 
# 
# ############################################################
# # Cell chat (split)
# ############################################################
# 
# # Split the seurat_rctd object by your experimental condition variable
# seurat_rctd.list <- SplitObject(seurat_rctd, split.by = "health")
# 
# # Prepare an empty list to store CellChat objects
# cellchat.list <- vector("list", length = length(seurat_rctd.list))
# names(cellchat.list) <- names(seurat_rctd.list)
# 
# for (i in seq_along(seurat_rctd.list)) {
#   # Step 1: Create a CellChat object for this condition
#   cellchat.list[[i]] <- createCellChat(object = seurat_rctd.list[[i]], 
#                                        group.by = "ident", assay = "RNA")
#   # Step 2: Assign the interaction database
#   cellchat.list[[i]]@DB <- CellChatDB.human
#   
#   # Step 3: Subset data (required step by CellChat)
#   cellchat.list[[i]] <- subsetData(cellchat.list[[i]])
#   
#   # Step 4: Identify overexpressed genes and interactions
#   cellchat.list[[i]] <- identifyOverExpressedGenes(cellchat.list[[i]])
#   cellchat.list[[i]] <- identifyOverExpressedInteractions(cellchat.list[[i]])
#   
#   # Step 5: Compute communication probabilities
#   cellchat.list[[i]] <- computeCommunProb(cellchat.list[[i]], type = "truncatedMean", trim = 0.1, 
#                                           distance.use = FALSE, interaction.range = 250, scale.distance = NULL, 
#                                           contact.dependent = TRUE, contact.range = 100)
#   
#   # Step 6: Filter communication
#   cellchat.list[[i]] <- filterCommunication(cellchat.list[[i]], min.cells = 10)
#   
#   subsetCommunication(cellchat.list[[i]]) %>% dplyr::arrange(pval) 
# }
# 
# # Merge all CellChat objects for comparison
# cellchat_merged <- mergeCellChat(cellchat.list, add.names = names(cellchat.list))
# 
# # Now you can run the comparison visualization
# netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE)