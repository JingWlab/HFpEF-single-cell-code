##############################################
# Script information                                                      
# Title: Sub-clustering analysis of cardiomyocyte
# Date: 2023-06-30
# Description: None
##############################################

# Load required packages
library(Seurat)
library(ggcorrplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

## Load required data (from 03.Clustering and cell type annotation)
S <- readRDS("Integrated_HF_Annotated.rds")

# Subset cardiomyocytes (CM) from the full dataset
S.CM <- subset(S, annot_corrected == "CM")

# Set the normalized count variance assay for feature selection
DefaultAssay(S.CM) <- "RNA_cv"  

# ----------------------------
# 3. DIMENSIONALITY REDUCTION
# ----------------------------
# Step 1: Find variable features
S.CM <- FindVariableFeatures(S.CM)

# Step 2: Switch to integrated assay for clustering
DefaultAssay(S.CM) <- "integrated"

# Step 3: Scale all genes and run PCA
all.genes <- rownames(S.CM)
S.CM <- ScaleData(S.CM, features = all.genes)
S.CM <- RunPCA(S.CM)

# Step 4: Determine significant PCs using elbow plot
ElbowPlot(S.CM)  # Select PCs where curve plateaus (e.g., 1-15)

# Step 5: Clustering and UMAP visualization
S.CM <- FindNeighbors(S.CM, dims = 1:50)      # Use sufficient PCs for neighborhood graph
S.CM <- FindClusters(S.CM, resolution = 0.6)  # Adjust resolution for cluster granularity
S.CM <- RunUMAP(S.CM, dims = 1:15)            # Use biologically relevant PCs

# Visualize clusters
DimPlot(S.CM, split.by = "stim", label = TRUE)  # By treatment condition
DimPlot(S.CM, group.by = "stim")               # Combined view by condition
DimPlot(S.CM, label = TRUE)                    # Default cluster view

# ----------------------------
# 4. SUBCLUSTER MERGING BY CORRELATION
# ----------------------------
# Identify significant marker genes for all clusters
CM.markers <- FindAllMarkers(S.CM, only.pos = TRUE)
sig.markers <- subset(CM.markers, 
                     p_val_adj < 1e-6 & avg_log2FC >= log2(1.5))
marker.list <- unique(sig.markers$gene)

# Create cluster IDs and set identities
S.CM$CM_num <- paste0("CM", S.CM$integrated_snn_res.0.6)
Idents(S.CM) <- "CM_num"

# Compute average expression matrix
avg.exp <- AverageExpression(S.CM, assays = "RNA_cv", 
                            features = marker.list)$RNA_cv
avg.exp <- avg.exp[, order(colnames(avg.exp))]  # Order clusters

# Calculate cluster correlations
corr.matrix <- cor(avg.exp)
p.matrix <- cor_pmat(avg.exp)

# FIGURE 2C PREVIEW: Cluster correlation heatmap
ggcorrplot(corr.matrix, 
           type = "lower",
           lab = TRUE, 
           p.mat = p.matrix,
           hc.order = TRUE,
           colors = c("#6D9EC1", "white", "#E46726"),
           ggtheme = theme_gray()) +
  ggtitle("Cardiomyocyte Subcluster Correlations")

# ----------------------------
# 5. FINAL SUBCLUSTER IDENTIFICATION
# ----------------------------
# MANUAL STEP: Merge highly correlated clusters based on heatmap
# Result: 8 distinct subclusters named CM1-CM8 (stored in metadata 'name')
# Example merging code (customize based on your correlation results):
# S.CM$name <- ifelse(S.CM$CM_num %in% c("CM1", "CM3"), "CM1", S.CM$CM_num)

# ----------------------------
# 6. VISUALIZATION OF FINAL CLUSTERS
# ----------------------------
## FIGURE 2C: Final UMAP ##
Idents(S.CM) <- "name"
DimPlot(S.CM, group.by = "name", label = TRUE, label.size = 6) +
  ggtitle("Cardiomyocyte Subclusters")

## FIGURE 2D: Marker Expression Dot Plot ##
marker_genes <- c("Pcdh7", "Ndufa4", "Nmrk2", "Tmem176b", 
                 "Dlc1", "Ndufa5", "Xirp2", "Nppa")

DotPlot(S.CM, features = marker_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Subcluster Marker Expression")

# ----------------------------
# 7. FUNCTIONAL ANALYSIS (FIG 2E-G)
# ----------------------------
# Set assay for differential expression
DefaultAssay(S.CM) <- "RNA_cv"
Idents(S.CM) <- "name"

# Find markers for final clusters
final.markers <- FindAllMarkers(S.CM, only.pos = TRUE)

# Function for enrichment analysis
run_cluster_enrichment <- function(cluster_id, marker_df = final.markers, top_n = 50) {
  # Get top marker genes
  cluster_genes <- subset(marker_df, cluster == cluster_id)$gene
  if (length(cluster_genes) > top_n) {
    cluster_genes <- cluster_genes[1:top_n]
  }
  
  # Convert to Entrez IDs
  gene_ids <- bitr(cluster_genes, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")
  
  # GO enrichment (Biological Process)
  go_bp <- enrichGO(gene = gene_ids$ENTREZID,
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    readable = TRUE)
  
  # KEGG pathway enrichment
  kegg <- enrichKEGG(gene = gene_ids$ENTREZID,
                     organism = "mmu",
                     pvalueCutoff = 0.05)
  
  # Save results
  write.csv(go_bp@result, file = paste0(cluster_id, "_GO_BP.csv"))
  write.csv(kegg@result, file = paste0(cluster_id, "_KEGG.csv"))
  
  # Generate plots
  go_plot <- dotplot(go_bp, showCategory = 15, title = paste(cluster_id, "GO-BP"))
  kegg_plot <- dotplot(kegg, showCategory = 15, title = paste(cluster_id, "KEGG"))
  
  return(list(go = go_bp, kegg = kegg, plots = list(go_plot, kegg_plot)))
}

# Run enrichment for all clusters
for (cluster in paste0("CM", 1:8)) {
  results <- run_cluster_enrichment(cluster)
  
  # Save plots
  ggsave(paste0(cluster, "_GO_plot.pdf"), results$plots[[1]], width = 10, height = 8)
  ggsave(paste0(cluster, "_KEGG_plot.pdf"), results$plots[[2]], width = 10, height = 8)
}

# Note: Top 4 pathways from each cluster were selected for Figure 2E-G
# Final figure assembly was completed in Excel