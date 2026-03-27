############################################################
# DE and Pathway Analysis
# Visium HD spatial transcriptomics workshop
# Author: Harvard Chan Bioinformatics Core
############################################################

library(Seurat)
library(tidyverse)
library(qs)
library(EnhancedVolcano)

############################################################
# Load processed Seurat object
############################################################

seurat_rctd <- qread("intermediates/07_seurat_rctd_unassigned.qs")

############################################################
# Macrophage findmarkers
############################################################

seurat_macro <- subset(seurat_rctd,
                       first_type == "Macrophage")
Idents(seurat_macro) <- "orig.ident"

# Calculate DGEs
dge <- FindMarkers(seurat_macro,
                   ident.1 = "P5CRC",
                   ident.2 = "P5NAT")
dge$gene <- rownames(dge)
View(dge)


############################################################
# Visualizations
############################################################

# Volcano plot
EnhancedVolcano(dge,
                        lab=dge$gene,
                        x="avg_log2FC",
                        y="p_val_adj",
                        title="FindMarkers Macrophage bins",
                        subtitle="CRC vs NAT")

# Filter significant genes
dge_sig <- dge %>% subset(p_val_adj < 0.05)


# Get the gene names and get the first 6 values
# Ignore ribosomal genes
genes <- dge_sig %>%
  filter(!str_detect(gene, "Rpl|Rps")) %>% 
  head(6)
genes <- genes$gene


VlnPlot(seurat_macro, genes,
        group.by = "orig.ident")
FeaturePlot(seurat_macro, genes, reduction = "full.umap.sketch")

SpatialFeaturePlot(seurat_macro, genes[1], 
                   image.alpha = 0.1, pt.size.factor = 15)


############################################################
# ORA
############################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)

# Create background dataset for hypergeometric testing using all tested 
# genes for significance in the results
all_genes <- as.character(dge$gene)

# Extract significant results for up-regulated
sigUp <- dplyr::filter(dge_sig, 
                       p_val_adj < 0.05, 
                       avg_log2FC > 0)
sigUp_genes <- as.character(sigUp$gene)

# Run GO enrichment analysis 
egoUp <- enrichGO(gene = sigUp_genes, 
                  universe = all_genes,
                  keyType = "SYMBOL",
                  OrgDb = org.Hs.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, 
                  readable = TRUE)

# Output results from GO analysis to a table
cluster_summaryUp <- data.frame(egoUp)
View(cluster_summaryUp)

# Dotplot 
dotplot(egoUp, showCategory=20)

# Add similarity matrix to the termsim slot of enrichment result
egoUp <- enrichplot::pairwise_termsim(egoUp)

# Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
emapplot(egoUp, showCategory = 50)


############################################################
# GSEA
############################################################

# Use a specific collection; C5 GO signatures
m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>%
  dplyr::select(gs_name, gene_symbol)

# Extract the foldchanges
foldchanges <- dge_sig$avg_log2FC

# Name each fold change with the corresponding gene symbol
names(foldchanges) <- dge_sig$gene
# Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)

# Run GSEA
msig_GSEA <- GSEA(foldchanges, TERM2GENE = m_t2g, verbose = FALSE)

# Extract the GSEA results
msigGSEA_results <- msig_GSEA@result

# Look at results ordered by NES
msigGSEA_results %>% arrange(-NES) %>% View()
