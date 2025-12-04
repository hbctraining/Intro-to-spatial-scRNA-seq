#!/bin/bash

###############3
# Mouse brain dataset from 10x
## Fresh Frozen
## https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-mouse-brain-fresh-frozen
#mkdir -p fresh_frozen
#curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fresh_Frozen/Visium_HD_Mouse_Brain_Fresh_Frozen_web_summary.html --output-dir fresh_frozen
#curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fresh_Frozen/Visium_HD_Mouse_Brain_Fresh_Frozen_cloupe_008um.cloupe --output-dir fresh_frozen
#curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fresh_Frozen/Visium_HD_Mouse_Brain_Fresh_Frozen_feature_slice.h5 --output-dir fresh_frozen
#curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fresh_Frozen/Visium_HD_Mouse_Brain_Fresh_Frozen_metrics_summary.csv --output-dir fresh_frozen
#curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fresh_Frozen/Visium_HD_Mouse_Brain_Fresh_Frozen_molecule_info.h5 --output-dir fresh_frozen
#curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fresh_Frozen/Visium_HD_Mouse_Brain_Fresh_Frozen_spatial.tar.gz --output-dir fresh_frozen
#curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fresh_Frozen/Visium_HD_Mouse_Brain_Fresh_Frozen_binned_outputs.tar.gz --output-dir fresh_frozen
#
## Unzip bin folder
#tar -xvf fresh_frozen/Visium_HD_Mouse_Brain_Fresh_Frozen_binned_outputs.tar.gz --directory fresh_frozen
#tar -xvf fresh_frozen/Visium_HD_Mouse_Brain_Fresh_Frozen_spatial.tar.gz --directory fresh_frozen
#
## Fixed Frozen
## https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-mouse-brain-fixed-frozen
#mkdir -p fixed_frozen
#curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_web_summary.html --output-dir fixed_frozen
#curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_cloupe_008um.cloupe --output-dir fixed_frozen
#curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_feature_slice.h5 --output-dir fixed_frozen
#curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_metrics_summary.csv --output-dir fixed_frozen
#curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_molecule_info.h5 --output-dir fixed_frozen
#curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_spatial.tar.gz --output-dir fixed_frozen
#curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_binned_outputs.tar.gz --output-dir fixed_frozen
#
## Unzip bin folder
#tar -xvf fixed_frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_binned_outputs.tar.gz --directory fixed_frozen
#tar -xvf fixed_frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_spatial.tar.gz --directory fixed_frozen


################
# CRC dataset from 10x
# https://www.10xgenomics.com/platforms/visium/product-family/dataset-human-crc

# Sample P5 CRC
# Output Files
#curl --output-dir P5CRC/ -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer_P5/Visium_HD_Human_Colon_Cancer_P5_binned_outputs.tar.gz
#curl --output-dir P5CRC/ -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer_P5/Visium_HD_Human_Colon_Cancer_P5_spatial.tar.gz
# curl --output-dir P5CRC/ -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer_P5/Visium_HD_Human_Colon_Cancer_P5_molecule_info.h5
# curl --output-dir P5CRC/ -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer_P5/Visium_HD_Human_Colon_Cancer_P5_cloupe_008um.cloupe
# curl --output-dir P5CRC/ -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer_P5/Visium_HD_Human_Colon_Cancer_P5_feature_slice.h5
# curl --output-dir P5CRC/ -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer_P5/Visium_HD_Human_Colon_Cancer_P5_metrics_summary.csv

# Sample P5 NAT
# Output Files
#curl --output-dir P5NAT -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Normal_P5/Visium_HD_Human_Colon_Normal_P5_binned_outputs.tar.gz
#curl --output-dir P5NAT -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Normal_P5/Visium_HD_Human_Colon_Normal_P5_spatial.tar.gz
# curl --output-dir P5NAT -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Normal_P5/Visium_HD_Human_Colon_Normal_P5_molecule_info.h5
# curl --output-dir P5NAT -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Normal_P5/Visium_HD_Human_Colon_Normal_P5_cloupe_008um.cloupe
# curl --output-dir P5NAT -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Normal_P5/Visium_HD_Human_Colon_Normal_P5_feature_slice.h5
# curl --output-dir P5NAT -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Normal_P5/Visium_HD_Human_Colon_Normal_P5_metrics_summary.csv


# Unzip P5 samples
#tar -xvf P5NAT/Visium_HD_Human_Colon_Normal_P5_binned_outputs.tar.gz -C P5NAT
#tar -xvf P5NAT/Visium_HD_Human_Colon_Normal_P5_spatial.tar.gz -C P5NAT

#tar -xvf P5CRC/Visium_HD_Human_Colon_Cancer_P5_binned_outputs.tar.gz -C P5CRC
#tar -xvf P5CRC/Visium_HD_Human_Colon_Cancer_P5_spatial.tar.gz -C P5CRC

# Published metadata for cells (barcodes are in completely different format)
#wget https://github.com/10XGenomics/HumanColonCancer_VisiumHD/raw/refs/heads/main/MetaData/SingleCell_MetaData.csv.gz
#gunzip SingleCell_MetaData.csv.gz

# Deconvolution results
#wget https://github.com/10XGenomics/HumanColonCancer_VisiumHD/raw/refs/heads/main/MetaData/DeconvolutionResults_P5CRC.csv.gz
#gunzip DeconvolutionResults_P5CRC.csv.gz
#wget https://github.com/10XGenomics/HumanColonCancer_VisiumHD/raw/refs/heads/main/MetaData/P5CRC_Metadata.parquet

