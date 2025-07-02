##############################################
# Script information                                                      
# Title: HFpEF and HFrEF data quality control and filtering
# Date: 2023-06-10
# Description: None
##############################################

# Load required packages
library(dplyr)
library(Seurat)
library(stringr)

# Load gene annotation data
gene_info <- read.csv("/data/rawdata/HF/HFpEF/gene_info.csv", row.names = 1)
gene_info$Gene_new <- as.character(gene_info$Gene_Name)

# Handle duplicate gene names by appending gene length
gene_info$Gene_new[55494:55502] <- rownames(gene_info[55494:55502,])
duplicated_genes <- which(duplicated(gene_info$Gene_new))
gene_info$Gene_new[duplicated_genes] <- paste(
  gene_info$Gene_new[duplicated_genes], 
  gene_info$Gene_Length[duplicated_genes], 
  sep = "."
)

# Identify data files
exp_mtx_files <- list.files(
  path = '/data/rawdata/HF/HFpEF',
  pattern = "*_genematrix.csv", 
  full.names = TRUE
)
well_info_files <- list.files(
  path = '/data/rawdata/HF/HFpEF',
  pattern = "*_all.csv", 
  full.names = TRUE
)

# Initialize Seurat object
HFpEF <- NULL

# Process each sequencing chip dataset
for (i in seq_along(exp_mtx_files)) {
  # Load expression matrix
  exp_mtx <- as.matrix(read.csv(exp_mtx_files[i], header = TRUE, row.names = 1))
  
  # Extract chip identifier
  chip_idx <- str_sub(exp_mtx_files[i], 21, 26)
  
  # Load well information
  well_info <- read.csv(well_info_files[i], header = FALSE)
  colnames(well_info) <- c("Barcode", "Cell_info", "Counts")
  well_info <- well_info[1:ncol(exp_mtx), ]
  rownames(well_info) <- well_info$Barcode
  
  # Create unique cell identifiers
  well_info$Barcode_new <- paste(
    well_info$Barcode, 
    well_info$Cell_info, 
    chip_idx, 
    sep = "_"
  )
  
  # Update matrix column names
  colnames(exp_mtx) <- well_info[colnames(exp_mtx), "Barcode_new"]
  rownames(well_info) <- well_info$Barcode_new
  rownames(exp_mtx) <- gene_info[rownames(exp_mtx), "Gene_new"]
  
  # Remove control wells
  well_info <- well_info[!well_info$Cell_info %in% c("Pos_Ctrl", "Neg_Ctrl"), ]
  exp_mtx <- exp_mtx[, rownames(well_info)]
  
  # Create Seurat object
  seur <- CreateSeuratObject(
    counts = exp_mtx, 
    project = "HFpEF_CM", 
    meta.data = well_info
  )
  seur$chip_idx <- chip_idx
  
  # Merge datasets
  if (is.null(HFpEF)) {
    HFpEF <- seur
  } else {
    HFpEF <- merge(HFpEF, seur)
  }
}

# Basic QC filtering
HFpEF <- subset(HFpEF, nFeature_RNA > 500)  # Keep cells with >500 genes

# Load additional metadata
well_list_files <- list.files(
  path = '/data/rawdata/HF/HFpEF',
  pattern = "*_FD_WellList.TXT", 
  full.names = TRUE
)

combined_well_list <- NULL
for (file in well_list_files) {
  chip_idx <- str_sub(file, 21, 26)
  well_list <- read.table(file, sep = "\t", header = TRUE)
  well_list$Barcode1 <- paste0(
    str_sub(well_list$Barcode, 1, 8),
    str_sub(well_list$Barcode, 10, 17)
  )
  well_list$Barcode_new <- paste(
    well_list$Barcode1, 
    well_list$Sample, 
    chip_idx, 
    sep = "_"
  )
  combined_well_list <- rbind(combined_well_list, well_list)
}

# Add metadata to Seurat object
rownames(combined_well_list) <- combined_well_list$Barcode_new
HFpEF <- AddMetaData(
  HFpEF, 
  combined_well_list[colnames(HFpEF), ]
)

# Calculate mitochondrial and ribosomal percentages
HFpEF[['percent.mt']] <- PercentageFeatureSet(HFpEF, pattern = '^mt-')
HFpEF[['percent.ribo']] <- PercentageFeatureSet(HFpEF, pattern = '^Rp[sl]') + 1

# Apply QC filters based on mitochondrial percentage
mt_threshold_high <- mean(HFpEF$percent.mt) + 1.5 * sd(HFpEF$percent.mt)
mt_threshold_low <- mean(HFpEF$percent.mt) - 1.5 * sd(HFpEF$percent.mt)
HFpEF <- subset(
  HFpEF, 
  State == "Good" & 
  percent.mt > mt_threshold_low & 
  percent.mt < mt_threshold_high
)

# Remove mitochondrial genes
mt_genes <- grep("mt-", rownames(HFpEF))[-1]  # Keep first mt gene if present
HFpEF[["RNA_noMT"]] <- CreateAssayObject(
  counts = GetAssayData(HFpEF, assay = "RNA")[-mt_genes, ]
)
DefaultAssay(HFpEF) <- "RNA_noMT"

# Normalization and clustering
HFpEF <- NormalizeData(HFpEF, scale.factor = 10000)
HFpEF <- FindVariableFeatures(HFpEF)
HFpEF <- ScaleData(HFpEF)
HFpEF <- RunPCA(HFpEF, features = VariableFeatures(HFpEF))
HFpEF <- FindNeighbors(HFpEF, dims = 1:50)
HFpEF <- FindClusters(HFpEF, resolution = 0.3)
HFpEF <- RunUMAP(HFpEF, dims = 1:30)


# Load HFrEF dataset
HFrEF_counts <- read.csv("GSE120064_TAC_raw_umi_matrix.csv", row.names = 1)
HFrEF_meta <- read.table(
  "GSE120064_TAC_clean_cell_info_summary.txt", 
  row.names = 1, 
  header = TRUE
)

# Create HFrEF Seurat object
HFrEF <- CreateSeuratObject(
  counts = HFrEF_counts, 
  meta.data = HFrEF_meta, 
  project = "HFrEF_LW"
)

# Calculate mitochondrial percentage for HFrEF
HFrEF[['percent.mt']] <- PercentageFeatureSet(HFrEF, pattern = '^mt-')

# Apply QC filters to HFrEF
mt_percent <- HFrEF$percent.mt
hfref_mt_high <- mean(mt_percent) + 1.5 * sd(mt_percent)
hfref_mt_low <- mean(mt_percent) - 1.5 * sd(mt_percent)
HFrEF <- subset(
  HFrEF, 
  percent.mt > hfref_mt_low & 
  percent.mt < hfref_mt_high
)

# Create clean assays by removing unwanted genes
remove_genes <- function(object, patterns) {
  genes_to_remove <- unique(unlist(sapply(patterns, function(p) grep(p, rownames(object)))))
  CreateAssayObject(counts = GetAssayData(object)[-genes_to_remove, ])
}

# Process both datasets
HFrEF[["RNA_cv"]] <- remove_genes(HFrEF, c("Hba", "Hbb", "mt-"))
HFpEF[["RNA_cv"]] <- remove_genes(HFpEF, c("mt-", "Hba", "Hbb"))

# Set default assays
DefaultAssay(HFpEF) <- "RNA_cv"
DefaultAssay(HFrEF) <- "RNA_cv"

# Normalization and feature selection
HFpEF <- NormalizeData(HFpEF) %>% FindVariableFeatures()
HFrEF <- NormalizeData(HFrEF) %>% FindVariableFeatures()

# Save processed data
saveRDS(HFpEF, file = "HFpEF.rds")
saveRDS(HFrEF, file = "HFrEF.rds")