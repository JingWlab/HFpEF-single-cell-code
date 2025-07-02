##############################################
# Script information                                                      
# Title: Clustering and cell type annotation
# Date: 2023-06-13
# Description: None
##############################################

# Load required packages
library(Seurat)

## Load required data (from 02.Integration of Datasets.)
integrated <- readRDS("integrated.rds")

# Visualize clustering at different resolutions
resolutions <- c(0, 0.03, 0.06, 0.12, 0.2, 0.3, 0.4, 0.5, 0.8, 1, 2, 3)
for (res in resolutions) {
  integrated <- FindClusters(
    integrated, 
    resolution = res, 
    graph.name = paste0("integrated_snn_res.", res)
  )
}

# Switch to RNA assay for marker visualization
DefaultAssay(integrated) <- "RNA_cv"

# Manual annotation based on marker genes
#     "Cd79a", "Cd79b", "Ly6d",           # B-cell markers
#     "Tnnc1", "Actn2", "Myh6",           # Cardiomyocyte markers
#     "Fabp4", "Cdh5", "Pecam1", "Ece1",  # Endothelial markers
#     "Tm4sf1", "Vwf",                    # Andocardial cells markers
#     "Col1a2", "Dpt", "Col3a1",          # Fibroblast markers
#     "S100a8", "S100a9",                 # Granulocyte markers
#     "Hbb-bt", "Hba-a2", "Hbb-bs",       # Hemoglobin genes
#     "Ctss", "C1qa", "Cd68",             # Myeloid markers
#     "Cd3d", "Cd3e", "Cd3g"              # T-cell markers


integrated$annot <- "CM"  # Default to cardiomyocytes
integrated$annot[integrated$integrated_snn_res.0.8 %in% c(4,18)] <- "EC"      # Endothelial cells
integrated$annot[integrated$integrated_snn_res.0.8 == 11] <- "T"              # T-cells
integrated$annot[integrated$integrated_snn_res.0.8 == 14] <- "GN"             # Granulocytes
integrated$annot[integrated$integrated_snn_res.0.8 %in% c(3,15,17)] <- "MP"   # Macrophages
integrated$annot[integrated$integrated_snn_res.0.8 %in% c(2,10,19,6,9)] <- "FB" # Fibroblasts
integrated$annot[integrated$integrated_snn_res.0.8 == 12] <- "EcC"            # Endocardial cells
integrated$annot[integrated$integrated_snn_res.0.8 == 16] <- "B"              # B-cells

# Save fully processed and annotated data
saveRDS(integrated, file = "Integrated_HF_Annotated.rds")

##########################################################
### Visualization Functions for Key Figures ###
##########################################################

# Figure 1B: Cell type visualization
plot_Fig1B <- function(seurat_obj) {
  DimPlot(
    seurat_obj,
    group.by = "annot",
    label = TRUE,
    label.size = 6,
    repel = TRUE
  )
}

# Figure 1C: Marker gene dot plot
plot_Fig1C <- function(seurat_obj) {
  DotPlot(
    seurat_obj,
    features = c(
      "Cd79a", "Cd79b", "Ly6d",           # B-cell markers
      "Tnnc1", "Actn2", "Myh6",           # Cardiomyocyte markers
      "Fabp4", "Cdh5", "Pecam1", "Ece1",  # Endothelial markers
      "Tm4sf1", "Vwf",                    # Endocardial cells endothelial
      "Col1a2", "Dpt", "Col3a1",          # Fibroblast markers
      "S100a8", "S100a9",                 # Granulocyte markers
      "Hbb-bt", "Hba-a2", "Hbb-bs",       # Hemoglobin genes
      "Ctss", "C1qa", "Cd68",             # Myeloid markers
      "Cd3d", "Cd3e", "Cd3g"              # T-cell markers
    ),
    group.by = "annot"
  ) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Generate figures
plot_Fig1B(integrated)
plot_Fig1C(integrated)
