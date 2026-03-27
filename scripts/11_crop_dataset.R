library(Seurat)
library(tidyverse)
library(qs)

seurat_rctd <- qread("intermediates/07_seurat_rctd_unassigned.qs")

Idents(seurat_rctd) <- "first_type"
SpatialDimPlot(seurat_rctd,
               image.alpha = 0.1,
               pt.size.factor = 12)


df <- rbind(
  GetTissueCoordinates(seurat_rctd,
                       image = "P5CRC.008um"),
  GetTissueCoordinates(seurat_rctd,
                       image = "P5NAT.008um"))

df <- FetchData(seurat_rctd, 
          c("cell", "orig.ident",
            "first_type", 
            "seurat_cluster.projected",
            "banksy_cluster",
            "SELENOP", "SPP1")) %>%
  left_join(df, by = "cell") %>%
  arrange("SELENOP")

########################
# NAT
########################
SpatialFeaturePlot(seurat_rctd,
                   features = "SELENOP",
                   images = "P5NAT.008um")

# Split in half vertically, right-hand side
# Coordinates appear rotated
df_nat <- df %>% 
  subset(orig.ident == "P5NAT") %>%
  subset(y > median(y)) %>%
  subset(x < 40000) %>%
  subset(x > 29000) %>%
  subset(y < 62000) %>%
  subset(y > 54000)
ggplot(df_nat,
       aes(x = x, y = y,
           color = SELENOP)) +
  geom_point(size = 0.4) +
  theme_bw()

nrow(df_nat)

ggplot(df_nat,
       aes(x = x, y = y,
           color = first_type)) +
  geom_point(size = 0.4) +
  theme_bw() +
  facet_wrap(~first_type) +
  NoLegend()

ggplot(df_nat,
       aes(x = first_type,
           fill = first_type)) +
  geom_bar() +
  theme_bw() + NoLegend() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

########################
# CRC
########################

df_crc <- df %>% 
  subset(orig.ident == "P5CRC") %>%
  subset(x > median(x)) %>%
  subset(y > 52000) %>%
  subset(x > 57500)

ggplot(df_crc,
       aes(x = x, y = y,
           color = SELENOP)) +
  geom_point(size = 0.4) +
  theme_bw()

nrow(df_crc)

ggplot(df_crc,
       aes(x = x, y = y,
           color = first_type)) +
  geom_point(size = 0.4) +
  theme_bw() +
  facet_wrap(~first_type) +
  NoLegend()

ggplot(df_crc,
       aes(x = first_type,
           fill = first_type)) +
  geom_bar() +
  theme_bw() + NoLegend() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


########################
# Cropping
########################

cropped_cells <- c(df_crc$cell,
                   df_nat$cell)
seurat_cropped <- seurat_rctd[, cropped_cells]

SpatialFeaturePlot(seurat_cropped,
                   features = "SELENOP",
                   pt.size.factor = 15)


###############################
# Code to create cropped files
###############################

library(Seurat)
library(tidyverse)
library(DropletUtils)

# Creating cropped sections for several bin sizes
bins <- c("008", "016")

###############################
# CRC
###############################

# Manual
folder <- "P5CRC"
x_min <- 57500
x_max <- 64002
y_min <- 52000
y_max <- 64319

# AUTO
dir.create(paste0("data/", folder, "_cropped/"), 
           recursive = TRUE)
dir.create(paste0("data/", folder, "_cropped/binned_outputs/"), 
           recursive = TRUE)

# Create bin specific folder
for (bin in bins) {
  path_bin <- paste0("data/", folder, "_cropped/binned_outputs/square_", bin, "um")
  dir.create(path_bin, recursive = TRUE)
  
  # Copy spatial images from original folder
  dir.create(paste0(path_bin, "/spatial"), recursive = TRUE)
  
  path_json <- paste0("data/", folder, "/binned_outputs/square_", bin, 
                      "um/spatial/scalefactors_json.json")
  path_parq <- paste0("data/", folder, "/binned_outputs/square_", bin, 
                      "um/spatial/tissue_positions.parquet")
  path_png  <- paste0("data/", folder, "/binned_outputs/square_", bin, 
                      "um/spatial/tissue_lowres_image.png")
  files_spatial <- c(path_json, path_parq, path_png)
  file.copy(files_spatial,
            paste0(path_bin, "/spatial"))
  
  # Load object
  obj <- Load10X_Spatial(
    data.dir = paste0("data/", folder),
    bin.size = c(as.integer(bin)),
    slice = folder)
  
  # Cells to keep based on coordinates
  coords <- GetTissueCoordinates(obj)
  cells <- coords %>% 
    subset(x > x_min) %>%
    subset(x < x_max) %>%
    subset(y > y_min) %>%
    subset(y < y_max) %>% 
    row.names()
  
  # Subset object
  obj$filter <- TRUE
  obj@meta.data[cells, "filt"] <- FALSE
  obj_filt <- subset(obj, filt == FALSE)
  
  p <- SpatialFeaturePlot(obj, paste0("nCount_Spatial.", bin, "um"), 
                          pt.size.factor = 12) +
    SpatialFeaturePlot(obj_filt, paste0("nFeature_Spatial.", bin, "um"), 
                       pt.size.factor = 16)
  print(p)
  
  path_h5 <- paste0(path_bin, "/filtered_feature_bc_matrix.h5")
  write10xCounts(path_h5,
                 x=LayerData(obj_filt),
                 barcodes=colnames(obj_filt),
                 gene.id=rownames(obj_filt),
                 version="3",
                 type="HDF5",
                 overwrite = TRUE)
  
  rm(obj)
  rm(obj_filt)
}

###############################
# NAT
###############################


# Manual
folder <- "P5NAT"
x_min <- 29000
x_max <- 39987
y_min <- 54000
y_max <- 61991

# AUTO
dir.create(paste0("data/", folder, "_cropped/"), 
           recursive = TRUE)
dir.create(paste0("data/", folder, "_cropped/binned_outputs/"), 
           recursive = TRUE)

# Create bin specific folder
for (bin in bins) {
  path_bin <- paste0("data/", folder, "_cropped/binned_outputs/square_", bin, "um")
  dir.create(path_bin, recursive = TRUE)
  
  # Copy spatial images from original folder
  dir.create(paste0(path_bin, "/spatial"), recursive = TRUE)
  
  path_json <- paste0("data/", folder, "/binned_outputs/square_", bin, 
                      "um/spatial/scalefactors_json.json")
  path_parq <- paste0("data/", folder, "/binned_outputs/square_", bin, 
                      "um/spatial/tissue_positions.parquet")
  path_png  <- paste0("data/", folder, "/binned_outputs/square_", bin, 
                      "um/spatial/tissue_lowres_image.png")
  files_spatial <- c(path_json, path_parq, path_png)
  file.copy(files_spatial,
            paste0(path_bin, "/spatial"))
  
  # Load object
  obj <- Load10X_Spatial(
    data.dir = paste0("data/", folder),
    bin.size = c(as.integer(bin)),
    slice = folder)
  
  # Cells to keep based on coordinates
  coords <- GetTissueCoordinates(obj)
  cells <- coords %>% 
    subset(x > x_min) %>%
    subset(x < x_max) %>%
    subset(y > y_min) %>%
    subset(y < y_max) %>% 
    row.names()
  
  # Subset object
  obj$filter <- TRUE
  obj@meta.data[cells, "filt"] <- FALSE
  obj_filt <- subset(obj, filt == FALSE)
  
  p <- SpatialFeaturePlot(obj, paste0("nCount_Spatial.", bin, "um"), 
                          pt.size.factor = 12) +
    SpatialFeaturePlot(obj_filt, paste0("nFeature_Spatial.", bin, "um"), 
                       pt.size.factor = 16)
  print(p)
  
  path_h5 <- paste0(path_bin, "/filtered_feature_bc_matrix.h5")
  write10xCounts(path_h5,
                 x=LayerData(obj_filt),
                 barcodes=colnames(obj_filt),
                 gene.id=rownames(obj_filt),
                 version="3",
                 type="HDF5",
                 overwrite = TRUE)
  
  rm(obj)
  rm(obj_filt)
}