##############################################
# Script information                                                      
# Title: Metabolic pathway analysis in cardiomyocyte
# Date: 2023-07-30
# Description: None
##############################################

# Install core packages
install.packages(c("devtools", "tidyverse", "data.table", "wesanderson", "Seurat"))
BiocManager::install(c("UCell", "AUCell", "GSEABase", "GSVA", "phangorn", "rsvd"))

# Install specialized packages
devtools::install_github("YosefLab/VISION@v2.1.0")        # Version-controlled installation
devtools::install_github("wu-yc/scMetabolism")           # Single-cell metabolism toolkit
devtools::install_github("saeyslab/nichenetr")           # Gene symbol conversion

# Load required packages
library(Seurat)              
library(scMetabolism)        # Metabolic pathway scoring
library(nichenetr)           # Ortholog conversion (mouse-human)
library(tidyverse)           # Data manipulation and visualization
library(AUCell)              # Gene set enrichment analysis
library(ComplexHeatmap)      # Advanced heatmap visualization

#' Convert mouse gene symbols to human for metabolic analysis
#' 
#' @param obj Seurat object with mouse gene symbols
#' @param metabolism.type Pathway database ("KEGG" or "REACTOME")
#' @return Seurat object with metabolic scores added to metadata
run_mouse_metabolism_analysis <- function(obj, metabolism.type = "REACTOME") {
  
  # Convert mouse to human gene symbols
  gene_trans <- rownames(obj) %>% 
    nichenetr::convert_mouse_to_human_symbols() %>% 
    as.data.frame()
  
  # Create mapping data frame
  gene_mouse <- data.frame(mouse = rownames(obj))
  gene_use <- cbind(gene_trans, gene_mouse)
  colnames(gene_use) <- c('human', 'mouse')
  
  # Remove genes without human orthologs
  gene_use <- na.omit(gene_use)
  
  # Subset to conserved genes
  mouse_data_trans <- subset(obj, features = gene_use$mouse)
  
  # Rename genes to human symbols
  mouse_data_trans <- RenameGenesSeurat(
    obj = mouse_data_trans,
    newnames = gene_use$human,
    gene.use = gene_use$mouse,
    de.assay = 'RNA'
  )
  
  # Calculate metabolic pathway scores using AUCell method
  mouse_metabolism <- sc.metabolism.Seurat(
    obj = mouse_data_trans,
    method = "AUCell",
    imputation = FALSE,        # No imputation for missing values
    ncores = 2,                # Use 2 cores for parallel processing
    metabolism.type = metabolism.type
  )
  
  return(mouse_metabolism)
}

#' Rename genes in all slots of a Seurat object
#' 
#' @param obj Seurat object
#' @param newnames New gene symbols (character vector)
#' @param gene.use Genes to keep (NULL keeps all)
#' @param de.assay Default assay for dimension reduction
RenameGenesSeurat <- function(obj, newnames, gene.use = NULL, de.assay) {
  # Validate input lengths
  if (is.null(gene.use)) {
    gene.use <- rownames(obj)
  } else {
    obj <- subset(obj, features = gene.use)
  }
  
  if (length(newnames) != length(gene.use)) {
    stop("Gene vectors must have the same length")
  }
  
  # Create gene mapping data frame
  gene_map <- data.frame(
    original = gene.use,
    new = make.unique(newnames)
  rownames(gene_map) <- gene_map$original
  
  # Process each assay in the object
  assays <- Seurat::Assays(obj)
  for (assay_name in assays) {
    assay_obj <- obj[[assay_name]]
    
    # Update feature names in count matrices
    matched_genes <- intersect(rownames(assay_obj), gene_map$original)
    if (length(matched_genes) > 0) {
      new_names <- gene_map[matched_genes, "new"]
      
      # Update counts and data slots
      rownames(assay_obj@counts) <- new_names
      rownames(assay_obj@data) <- new_names
      
      # Update variable features
      var_features <- VariableFeatures(obj, assay = assay_name)
      if (length(var_features) > 0) {
        new_var <- gene_map[var_features, "new"]
        VariableFeatures(obj, assay = assay_name) <- new_var
      }
      
      # Update scale.data if exists
      if (!is.null(assay_obj@scale.data)) {
        scaled_genes <- intersect(rownames(assay_obj@scale.data), gene_map$original)
        rownames(assay_obj@scale.data) <- gene_map[scaled_genes, "new"]
      }
    }
    obj[[assay_name]] <- assay_obj
  }
  
  # Update dimension reduction loadings
  if ("pca" %in% names(obj@reductions)) {
    pca_genes <- rownames(obj@reductions$pca@feature.loadings)
    matched_pca <- intersect(pca_genes, gene_map$original)
    if (length(matched_pca) > 0) {
      rownames(obj@reductions$pca@feature.loadings) <- gene_map[matched_pca, "new"]
    }
  }
  
  return(obj)
}

# MAIN ANALYSIS WORKFLOW
# ----------------------------
# Step 1: Subset HFpEF cardiomyocytes
Idents(S.CM) <- "stim"                          # Set identity to treatment groups
CM_HFpEF <- subset(S.CM, ident = "HFpEF")       # Select HFpEF cardiomyocytes

# Step 2: Perform metabolic pathway analysis
# Note: Uses REACTOME pathways (better coverage for metabolism)
mouse_CM_hfpef <- run_mouse_metabolism_analysis(
  obj = CM_HFpEF,
  metabolism.type = "REACTOME"
)

# Step 3: Add metabolic scores to original object
S.CM <- AddMetaData(S.CM, mouse_CM_hfpef@assays$METABOLISM@scale.data)

# ----------------------------
# 6. VISUALIZATION (FIGURE 2H)
# ----------------------------
# Create UMAP plot colored by fatty acid metabolism activity
fatty_acid_plot <- FeaturePlot(
  object = S.CM,
  features = "REACTOME_FATTY_ACID_METABOLISM",  # Pathway score column
  reduction = "umap",
  cols = c("lightgrey", "red"),                 # Low to high expression
  pt.size = 1,
  order = TRUE                                  # Plot high-expressing cells on top
) +
  ggtitle("Fatty Acid Metabolism in HFpEF Cardiomyocytes") +
  theme_minimal() +
  theme(legend.position = "right")

# Save plot
ggsave("Fig2H_FattyAcid_Metabolism.pdf", fatty_acid_plot, width = 8, height = 6)

# Note: Final figure was polished in Excel as per original workflow