#!/bin/bash

# ############## MOUSE BRAIN
# # Mouse brain dataset from 10x
# # Fresh Frozen
# # https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-mouse-brain-fresh-frozen
# mkdir -p fresh_frozen
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fresh_Frozen/Visium_HD_Mouse_Brain_Fresh_Frozen_web_summary.html --output-dir fresh_frozen
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fresh_Frozen/Visium_HD_Mouse_Brain_Fresh_Frozen_cloupe_008um.cloupe --output-dir fresh_frozen
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fresh_Frozen/Visium_HD_Mouse_Brain_Fresh_Frozen_feature_slice.h5 --output-dir fresh_frozen
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fresh_Frozen/Visium_HD_Mouse_Brain_Fresh_Frozen_metrics_summary.csv --output-dir fresh_frozen
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fresh_Frozen/Visium_HD_Mouse_Brain_Fresh_Frozen_molecule_info.h5 --output-dir fresh_frozen
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fresh_Frozen/Visium_HD_Mouse_Brain_Fresh_Frozen_spatial.tar.gz --output-dir fresh_frozen
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fresh_Frozen/Visium_HD_Mouse_Brain_Fresh_Frozen_binned_outputs.tar.gz --output-dir fresh_frozen
# 
# # Unzip bin folder
# tar -xvf fresh_frozen/Visium_HD_Mouse_Brain_Fresh_Frozen_binned_outputs.tar.gz --directory fresh_frozen
# tar -xvf fresh_frozen/Visium_HD_Mouse_Brain_Fresh_Frozen_spatial.tar.gz --directory fresh_frozen
# 
# # Fixed Frozen
# # https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-mouse-brain-fixed-frozen
# mkdir -p fixed_frozen
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_web_summary.html --output-dir fixed_frozen
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_cloupe_008um.cloupe --output-dir fixed_frozen
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_feature_slice.h5 --output-dir fixed_frozen
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_metrics_summary.csv --output-dir fixed_frozen
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_molecule_info.h5 --output-dir fixed_frozen
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_spatial.tar.gz --output-dir fixed_frozen
# curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_binned_outputs.tar.gz --output-dir fixed_frozen
# 
# # Unzip bin folder
# tar -xvf fixed_frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_binned_outputs.tar.gz --directory fixed_frozen
# tar -xvf fixed_frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_spatial.tar.gz --directory fixed_frozen


################ HUMAN CRC
# CRC dataset from 10x
# https://www.10xgenomics.com/platforms/visium/product-family/dataset-human-crc

mkdir -p P5CRC
mkdir -p P5NAT

# Sample P5 CRC
# Output Files
cd P5CRC
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer_P5/Visium_HD_Human_Colon_Cancer_P5_binned_outputs.tar.gz
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer_P5/Visium_HD_Human_Colon_Cancer_P5_spatial.tar.gz
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer_P5/Visium_HD_Human_Colon_Cancer_P5_molecule_info.h5
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer_P5/Visium_HD_Human_Colon_Cancer_P5_cloupe_008um.cloupe
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer_P5/Visium_HD_Human_Colon_Cancer_P5_feature_slice.h5
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer_P5/Visium_HD_Human_Colon_Cancer_P5_metrics_summary.csv
cd ../

# Sample P5 NAT
# Output Files
cd P5NAT
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Normal_P5/Visium_HD_Human_Colon_Normal_P5_binned_outputs.tar.gz
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Normal_P5/Visium_HD_Human_Colon_Normal_P5_spatial.tar.gz
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Normal_P5/Visium_HD_Human_Colon_Normal_P5_molecule_info.h5
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Normal_P5/Visium_HD_Human_Colon_Normal_P5_cloupe_008um.cloupe
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Normal_P5/Visium_HD_Human_Colon_Normal_P5_feature_slice.h5
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Normal_P5/Visium_HD_Human_Colon_Normal_P5_metrics_summary.csv
cd ../

# Unzip CRC files
tar -xvf P5NAT/Visium_HD_Human_Colon_Normal_P5_binned_outputs.tar.gz -C P5NAT
tar -xvf P5NAT/Visium_HD_Human_Colon_Normal_P5_spatial.tar.gz -C P5NAT
tar -xvf P5CRC/Visium_HD_Human_Colon_Cancer_P5_binned_outputs.tar.gz -C P5CRC
tar -xvf P5CRC/Visium_HD_Human_Colon_Cancer_P5_spatial.tar.gz -C P5CRC

# Published metadata for cells (barcodes are in completely different format)
wget https://github.com/10XGenomics/HumanColonCancer_VisiumHD/raw/refs/heads/main/MetaData/SingleCell_MetaData.csv.gz
gunzip SingleCell_MetaData.csv.gz

# Deconvolution results
wget https://github.com/10XGenomics/HumanColonCancer_VisiumHD/raw/refs/heads/main/MetaData/DeconvolutionResults_P5CRC.csv.gz
gunzip DeconvolutionResults_P5CRC.csv.gz
wget https://github.com/10XGenomics/HumanColonCancer_VisiumHD/raw/refs/heads/main/MetaData/P5CRC_Metadata.parquet


mkdir flex_reference
mv DeconvolutionResults_P5CRC.csv flex_reference/
mv SingleCell_MetaData.csv flex_reference/
mv P5CRC_Metadata.parquet flex_reference/


### Download CRC flex counts dataset
cd flex_reference
# Input Files

curl -O https://cf.10xgenomics.com/samples/cell-exp/8.0.0/HumanColonCancer_Flex_Multiplex/HumanColonCancer_Flex_Multiplex_aggregation.csv
curl -O https://cf.10xgenomics.com/samples/cell-exp/8.0.0/HumanColonCancer_Flex_Multiplex/HumanColonCancer_Flex_Multiplex_count_feature_reference.csv

# Output Files
curl -O https://cf.10xgenomics.com/samples/cell-exp/8.0.0/HumanColonCancer_Flex_Multiplex/HumanColonCancer_Flex_Multiplex_count_analysis.tar.gz
curl -O https://cf.10xgenomics.com/samples/cell-exp/8.0.0/HumanColonCancer_Flex_Multiplex/HumanColonCancer_Flex_Multiplex_count_filtered_feature_bc_matrix.h5
curl -O https://cf.10xgenomics.com/samples/cell-exp/8.0.0/HumanColonCancer_Flex_Multiplex/HumanColonCancer_Flex_Multiplex_count_filtered_feature_bc_matrix.tar.gz
curl -O https://cf.10xgenomics.com/samples/cell-exp/8.0.0/HumanColonCancer_Flex_Multiplex/HumanColonCancer_Flex_Multiplex_count_cloupe.cloupe
curl -O https://cf.10xgenomics.com/samples/cell-exp/8.0.0/HumanColonCancer_Flex_Multiplex/HumanColonCancer_Flex_Multiplex_count_summary.json
