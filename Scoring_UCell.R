library(Seurat)
library(ggplot2)
library(UCell)
library(dplyr)
library(cowplot)
library(patchwork)

#------------------------ Signature scoring with UCell in Seurat ------------------------#
#
# Why shift to UCell scoring (rather than using the native AddModuleScore function from Seurat?
# UCell addresses several limitations associated with rank-based scoring, making it - at least for me - more robust and reliable:
# 1. Robust to expression scale and distribution:
# 	- AddModuleScore: calculates the average expression of a set of genes in each spot and subtracts the average expression of a random set of control genes. However, this approach can be sensitive to differences in gene expression levels, normalization, and sequencing depth across cells.
#   - UCell: uses rank-based statistics (Mann-Whitney U) rather than raw expression values, making it robust to variations in data scale and distribution. 
# 2. Normalization and interpretability:
# 	- AddModuleScore: produces scores that are not normalized and can be difficult to interpret directly, as they depend on the specific gene set and control genes used. This lack of standardization makes it harder to compare scores across different datasets.
# 	- UCell: produces normalized scores between 0 and 1, which are easier to interpret and compare across cells, samples, or datasets. The scores represent the relative enrichment of a gene signature within each cell.
# 3. Control Gene Selection:
# 	- AddModuleScore: relies on selecting random control genes of similar expression levels to the genes in the signature, which can introduce variability and noise into the scores due to arbitrary control gene choices.
# 	- UCell: does not require control genes, avoiding the randomness and variability associated with their selection, leading to more consistent and reproducible scoring.
# 4. Suitability for sparse and noisy data (applied to both single-cell and Visium data!):
# 	- AddModuleScore: can be affected by the sparsity data as many genes have zero expression in some cells/spots, potentially leading to misleading scores.
# 	- UCell: is more robust to sparse data since it focuses on ranks rather than raw expression values. It can effectively handle the zero-inflated nature of single-cell data without being heavily biased by low expression counts.
#
# ---------------------------------------------------------------------------------------#

#---
#NOT needed if using Seurat v5 (authors of UCell included a patch to connect directly with Seurat)
#source("helpers/AddModuleScore_with_UCell.R")
# source helper scripts if needed
# source("helpers/AddModuleScore_with_UCell.R")
source("helpers/SpatialFeaturePlotBlend.R")

# Function to load data and run the analysis pipeline
run_signatures <- function(visium_obj_path, signatures_path) {
  
  #--- Load Visium object
  visium_obj <- readRDS(visium_obj_path)

  #--- Load PMCB signatures
  sigs <- read.table(signatures_path, header = TRUE, sep = "\t")
  signatures <- lapply(sigs, as.character)
  names(signatures) <- colnames(sigs)

  #--- Scoring with UCell
  visium_obj <- AddModuleScore_UCell(visium_obj, features = signatures, name = NULL)

  #--- k-Nearest Neighbors smoothing to refine UCell scores
  visium_obj <- SmoothKNN(visium_obj,
                          signature.names = names(signatures),
                          reduction = "pca")

  #--- Visualization
  SpatialFeaturePlot(visium_obj, 
                     features = c("Type1_immuneResponse", "Type1_immuneResponse_kNN"), 
                     ncol = 2,
                     alpha = c(0.4, 1),
                     pt.size = 2.5)

  SpatialFeaturePlot(visium_obj, 
                     features = c("IL.17_Signaling", "IL.17_Signaling_kNN"), 
                     ncol = 2,
                     alpha = c(0.4, 1),
                     pt.size = 2.5)

  #--- Test with blend
  SpatialFeaturePlotBlend(visium_obj, features = c("IL.17_Signaling_kNN", "IL1B"), pt.size = 3)

  return(visium_obj)
}

# Example of calling the function with paths as arguments
# run_signatures("TIS4606_PP_preprocessed.rds", "Spatial_Signatures.txt")

