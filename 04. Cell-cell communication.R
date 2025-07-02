##############################################
# Script information                                                      
# Title: Cell-cell communication analysis of HFpEF and HFrEF Datasets
# Date: 2023-06-20
# Description: None
##############################################

# Load required libraries
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)

## Load required data (from 03.Clustering and cell type annotation)
HF <- readRDS("Integrated_HF_Annotated.rds")

# ----------------------------
# 1. DATA PREPARATION
# ----------------------------
# Split dataset by condition (stim: HFpEF, HFpEF_Control, HFrEF, HFrEF_Control)
HF.split <- SplitObject(HF, split.by = "stim")

# Define cell types for analysis
cell_types <- c("CM", "EC", "FB", "EcC", "MP", "T", "Hbb_high")

# Subset relevant cell populations for each condition
HFrEF <- subset(HF.split$HFrEF, idents = cell_types)
HFrEF_Control <- subset(HF.split$HFrEF_Control, idents = cell_types)
HFpEF_Control <- subset(HF.split$HFpEF_Control, idents = cell_types)
HFpEF <- subset(HF.split$HFpEF, idents = cell_types)

# ----------------------------
# 2. CELLCHAT ANALYSIS FUNCTION
# ----------------------------
#' Perform CellChat analysis on a Seurat object
#' 
#' @param seurat_obj Seurat object
#' @return Processed CellChat object
run_cellchat_analysis <- function(seurat_obj) {
  # Create CellChat object
  cellchat <- createCellChat(
    object = seurat_obj@assays$RNA@data,
    meta = data.frame(labels = Idents(seurat_obj), 
    group.by = "labels"
  )
  
  # Use mouse signaling database
  CellChatDB <- CellChatDB.mouse
  
  # Select relevant database subsets
  CellChatDB.use <- subsetDB(CellChatDB, search = c(
    "Secreted Signaling", 
    "Cell-Cell Contact", 
    "ECM-Receptor"
  ))
  
  # Run analysis pipeline
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.mouse)
  cellchat <- computeCommunProb(cellchat, type = "triMean", trim = 0.1)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cellchat <- filterCommunication(cellchat, min.cells = 100)
  
  return(cellchat)
}

# ----------------------------
# 3. RUN ANALYSIS PER CONDITION
# ----------------------------
# Run CellChat for each condition
HFrEF.con <- run_cellchat_analysis(HFrEF_Control)
HFrEF.obj <- run_cellchat_analysis(HFrEF)
HFpEF.con <- run_cellchat_analysis(HFpEF_Control)
HFpEF.obj <- run_cellchat_analysis(HFpEF)

# ----------------------------
# 4. COMPARATIVE ANALYSIS
# ----------------------------
# Merge HFrEF vs Control
HFrEF_comparison <- mergeCellChat(
  list(Control = HFrEF.con, HFrEF = HFrEF.obj),
  add.names = c("Control", "HFrEF")
)

# Merge HFpEF vs Control
HFpEF_comparison <- mergeCellChat(
  list(Control = HFpEF.con, HFpEF = HFpEF.obj),
  add.names = c("Control", "HFpEF")
)

# ----------------------------
# 5. VISUALIZATION
# ----------------------------
## FIGURE 2A: HFrEF Signaling Network Comparison ##
# Count-based interaction heatmap
gg1 <- netVisual_heatmap(HFrEF_comparison, 
                         color.heatmap = "OrRd",
                         title.name = "Interaction Counts: HFrEF vs Control")

## FIGURE 2B: HFpEF Signaling Network Comparison ##
# Count-based interaction heatmap
gg2 <- netVisual_heatmap(HFpEF_comparison, 
                         color.heatmap = "OrRd",
                         title.name = "Interaction Counts: HFpEF vs Control")

# Interaction strength heatmap
gg3 <- netVisual_heatmap(HFpEF_comparison, 
                         measure = "weight",
                         color.heatmap = "Purples",
                         title.name = "Interaction Strength: HFpEF vs Control")

# Combine plots
combined_plot <- gg2 + gg3 + 
  plot_layout(ncol = 2) +
  plot_annotation(title = "HFpEF Cell-Cell Communication")

# Save plots
ggsave("Fig2A_HFrEF_heatmap.pdf", gg1, width = 10, height = 8)
ggsave("Fig2B_HFpEF_comparison.pdf", combined_plot, width = 16, height = 8)

# Note: Final figures were polished in Excel as per original workflow