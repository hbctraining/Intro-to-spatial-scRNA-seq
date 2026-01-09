library(tidyverse)
library(Seurat)
library(DropletUtils)
library(cowplot)
library(png)
library(grid)
library(magick)
library(arrow)

folder <- "P5CRC"
bin <- "8"

# Load object
obj <- Load10X_Spatial(
  data.dir = paste0("data/", folder),
  bin.size = c(as.integer(bin)),
  slice = folder)

obj$cell <- Cells(obj)
coords <- GetTissueCoordinates(obj)
obj@meta.data <- left_join(obj@meta.data, 
                           coords, by = "cell")
rownames(obj@meta.data) <- obj$cell

x_min <- 55000
x_max <- Inf
y_min <- 52000
y_max <- Inf

ggplot(obj@meta.data) +
  geom_point(aes(x = x, y = y, 
                 color = log10(nCount_Spatial.008um)),
             size = 0.01) +
  theme_bw() + NoLegend() +
  scale_color_viridis_c() +
  geom_vline(xintercept = x_min) +
  geom_hline(yintercept = y_min)



summary(coords)
# x               y             cell          
# Min.   :39508   Min.   :44292   Length:541968     
# 1st Qu.:45159   1st Qu.:49228   Class :character  
# Median :51098   Median :54088   Mode  :character  
# Mean   :51243   Mean   :54157                     
# 3rd Qu.:57255   3rd Qu.:59044                     
# Max.   :64065   Max.   :68739                     

# Cells to keep based on coordinates
cells <- coords %>% 
  subset((x > y_min & x < y_max) & 
           (y > x_min & y < x_max)) %>% 
  row.names()

# Subset object
obj$filter <- TRUE
obj@meta.data[cells, "filt"] <- FALSE
obj_filt <- subset(obj, filt == FALSE)


p2 <- ggplot(obj_filt@meta.data) +
  geom_point(aes(x = x, y = y, 
                 color = log10(nCount_Spatial.008um)),
             size = 0.3) +
  theme_void() + NoLegend() +
  scale_color_viridis_c()



p1 + p2

