##############################################
# Script information                                                      
# Title: HFpEF and HFrEF data quality control and filtering
# Date: 2023-06-10
# Description: None
##############################################

# Load required packages
library(Seurat)

## Load required data (from 01.Quality control and filtering)
HFpEF <- readRDS("HFpEF.rds")
HFrEF <- readRDS("HFrEF.rds")

# Add metadata
HFpEF$Model <- "HFpEF"
HFrEF$Model <- "HFrEF"

# Integrate datasets
integration_features <- SelectIntegrationFeatures(object.list = list(HFpEF, HFrEF))
anchors <- FindIntegrationAnchors(
  object.list = list(HFpEF, HFrEF),
  dims = 1:30,
  anchor.features = integration_features
)
integrated <- IntegrateData(anchors, dims = 1:30)

# Set default assay for downstream analysis
DefaultAssay(integrated) <- "integrated"

# Calculate mitochondrial percentage post-integration
integrated[['percent.mt']] <- PercentageFeatureSet(integrated, pattern = '^mt-')

# Dimensionality reduction and clustering
integrated <- ScaleData(integrated) %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.1) %>%
  RunUMAP(dims = 1:20) %>%
  RunTSNE(dims = 1:15)

# Add metadata for disease models
integrated$Model <- ifelse(
  integrated$orig.ident == "SC", 
  "HFrEF", 
  "HFpEF"
)

# Save processed data
saveRDS(integrated, file = "integrated.rds")