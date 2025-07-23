#!/bin/bash

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

# Fixed Frozen
# https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-mouse-brain-fixed-frozen
mkdir -p fixed_frozen
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_web_summary.html --output-dir fixed_frozen
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_cloupe_008um.cloupe --output-dir fixed_frozen
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_feature_slice.h5 --output-dir fixed_frozen
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_metrics_summary.csv --output-dir fixed_frozen
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_molecule_info.h5 --output-dir fixed_frozen
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_spatial.tar.gz --output-dir fixed_frozen
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Mouse_Brain_Fixed_Frozen/Visium_HD_Mouse_Brain_Fixed_Frozen_binned_outputs.tar.gz --output-dir fixed_frozen
