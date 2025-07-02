##############################################
# Script information                                                      
# Title: SCENIC_associated analysis & Figure 3 A-I & Figure S7-10
# Date: 2023-11-22
# Description: None
##############################################

# packages
library (SCENIC)
library (AUCell)
library (GENIE3)
library (NMF)
library (umap)
library (SeuratObject)
library (aertslab)
library (BiocGenerics)
library (DelayedArray)
library (SparseArray)
library (S4Vectors)

##### SCENIC_associated analysis #####
# Load sc-RNA data (Single-cell RNA-seq data)
X <- readRDS("Integrative_20230628.rds")

# Prepare cell meta information (metadata about the cells)
cellInfo <- data.frame(X@meta.data)
# Renaming columns in the metadata for clarity
colnames(cellInfo)[which(colnames(cellInfo)=="annot_corrected")] <- "celltype"  # Renaming 'annot_corrected' to 'celltype'
colnames(cellInfo)[which(colnames(cellInfo)=="stim")] <- "model"  # Renaming 'stim' to 'model'
colnames(cellInfo)[which(colnames(cellInfo)=="annot_stim")] <- "modelcell"  # Renaming 'annot_stim' to 'modelcell'

# Prepare expression matrix (raw gene expression counts for each cell)
exprMat <- as.matrix(X@assays$RNA@counts)

# Initialize SCENIC settings (setting up for the SCENIC analysis pipeline)
data(list="motifAnnotations_mgi_v9", package="RcisTarget")  # Load motif annotations for gene regulation
motifAnnotations_mgi <- motifAnnotations_mgi_v9
scenicOptions <- initializeScenic(org="mgi", dbDir="cisTarget_databases", nCores=10)  # Initialize SCENIC options with MGI dataset and set number of cores
saveRDS(scenicOptions, file="int/scenicOptions.Rds")  # Save the initialized SCENIC options

### Co-expression network analysis ###
# Log-transform the expression matrix for normalization
exprMat_log <- log2(exprMat + 1)

# Filter genes based on SCENIC options (removing genes that don't pass certain thresholds)
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]  # Create the filtered expression matrix based on the genes that passed filtering

# Run correlation analysis to identify gene co-expression
Correlation <- runCorrelation(exprMat_filtered, scenicOptions)

# Log-transform the filtered expression matrix for normalization
exprMat_filtered_log <- log2(exprMat_filtered + 1)

# Run GENIE3 algorithm to infer gene regulatory networks from the expression data
Genie3 <- runGenie3(exprMat_filtered_log, scenicOptions)

# Step 1: Run co-expression network analysis and create modules
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")  # Save the updated SCENIC options

# Step 2: Create regulons based on the co-expression network
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")  # Save the updated SCENIC options

# Step 3: Score cells based on the regulon activities
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

# Step 4: Binarize the AUC (Area Under Curve) values for the regulon activity
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

# Generate t-SNE plot of AUC values
tsneAUC(scenicOptions, aucType="AUC")

### Download the results and save them ###
# View the AUC value of each cell and save it to a CSV file
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
b <- aucell_regulonAUC@assays@data@listData$AUC  # Extract the AUC values
write.csv(b, file = "auc_eachregulon.csv")  # Save the AUC values to a CSV file

# Calculate the binary regulon activity (on/off for each cell's regulon)
minPerc <- .7  # Minimum percentage threshold for binary activity
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")  # Load binary regulon activity
write.csv(binaryRegulonActivity, file="onoff_eachcell_regulon.csv")  # Save the binary activity to a CSV file

### Calculate RSS (Residual Sum of Squares) to evaluate regulon activity per cell type ###
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation = cellInfo[colnames(regulonAUC), "celltype"])  # Calculate RSS for each cell type
write.csv(rss, file="celltypespec_reg.csv")  # Save the RSS values to a CSV file

### Subsequent processing was completed in Excel. ###


##### Figure 3 #####
##### Fig 3A  (input: data from above course -- "celltypespec_reg.csv" from "RSS calculation")
AA=readxl::read_xlsx("bubble.xlsx",sheet = 4)
AA$cellType= factor(AA$cellType,levels = c( "CM","EC","EcC","FB","GN","MP","B","T","Hbb_high"))
AA$Topic= factor(AA$Topic,levels = c("Hmgb3_extended (18g)",
                                     "Ppargc1a (78g)",
                                     "Rxrg_extended (512g)",
                                     "Nfe2l1_extended (107g)",
                                     "Mitf_extended (50g)",
                                     "Atf6 (157g)",
                                     "E2f6_extended (4759g)",
                                     "Rorc (33g)",
                                     "Elk3 (36g)",
                                     "Pparg (33g)",
                                     "Ets1 (307g)",
                                     "Klf7 (14g)",
                                     "Gata2_extended (631g)",
                                     "Foxc1_extended (12g)",
                                     "Plagl1_extended (16g)",
                                     "Nr2f2_extended (49g)",
                                     "Tal1_extended (10g)",
                                     "Tcf7l1 (55g)",
                                     "Tcf21 (95g)",
                                     "Zeb1 (16g)",
                                     "Creb3l2_extended (37g)",
                                     "Ar (12g)",
                                     "Mxd1 (43g)",
                                     "Nfe2 (41g)",
                                     "Nfil3_extended (77g)",
                                     "Runx2_extended (120g)",
                                     "Irf5 (66g)",
                                     "Spic (460g)",
                                     "Mafb (56g)",
                                     "Maf (73g)",
                                     "Zscan4c (18g)",
                                     "Pax5 (25g)",
                                     "Bcl11a (16g)",
                                     "Irf4 (45g)",
                                     "Eomes (18g)",
                                     "Stat4_extended (32g)",
                                     "Lef1 (16g)",
                                     "Runx3_extended (48g)",
                                     "Irf7 (53g)",
                                     "Gata1 (22g)",
                                     "Mxi1 (768g)",
                                     "Prdm16_extended (682g)"))
PP<-ggplot(data = AA, mapping = aes(x=Topic,y=cellType))+
  geom_point(aes(size=RSS,color=Z))+
  coord_flip()+ 
  scale_color_gradient(low="#DCDCDC",high ="red")+
  theme(legend.key.size=unit(0.5,'cm'))+
  theme(legend.text = element_text( size = 10,face = 'bold'))+
  theme(axis.text.x = element_text(size = 10, family = "myFont", hjust =0.5,face = 'bold'))+
  theme(axis.text.y = element_text(size = 0.5, family = "myFont", hjust =0.5,face = 'bold'))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour="#F5F5F5", size =0.05))+
  theme(panel.grid.minor = element_line(color = "green", size =0.05))
print(PP)


##### Fig 3B-C
 Data of Fig 3B-C were downloded during above "SCENIC_associated analysis", and drawn in excel.


##### Fig 3D
S <- readRDS("Integrative_20230628.rds")
Idents(S) <- S@meta.data$stim
X <-subset(S, subset=stim!= "TAC_progress")
DefaultAssay(X) <- "RNA"
### The expression of Ppargc1a in each cell type was extracted
VlnPlot(X,
        features = c("Ppargc1a"),
        split.by = "annot_corrected",
        pt.size = 0)
###Data was further sorted out and drawn in excel.


##### Fig 3E
N <- read.table("/NASdata/wuwq/HFpEF/ALL/HFp_CM_reg.txt",  
                header=T, 
                row.names=1,
                sep="\t") 
N<-as.matrix(N)
p <- pheatmap(N, scale="none",
              border="black", 
              cluster_cols = F, 
              cluster_rows = F,
              legend = T, 
              legend_breaks=c(1,3,5),
              legend_labels = c(1,3,5),
              fontsize_row = 5, 
              fontsize_col = 8,
              clustering_distance_rows = "euclidean",
              clustering_method="centroid",
              fontsize_number = 5,
              color =  colorRampPalette(c("white","pink","firebrick3"))(100))
print(p)


##### Fig 3F
#main color of the circle
top_languagesQ<-c(
  'Ppargc1a',
  'Atf6',
  'E2f6',
  'Mitf')
top_language_colorsQ <- c(
  '#efb306',
  '#e8351e',
  '#cd023d',
  '#0f8096',
  '#000000'
)
# name of the 4 circles
names(top_language_colorsQ) <- c(
  'Ppargc1a',
  'Atf6',
  'E2f6',
  'Mitf',
  'All'
)
#lines from inner to outer
edges1Q<-read.table("/NASdata/wuwq/HFpEF/ALL/edges1Q.txt", header=TRUE,sep = "")
edges1Q$total_code<-as.numeric(edges1Q$total_code)
edges2Q<-read.table("/NASdata/wuwq/HFpEF/ALL/edges2Q.txt", header=TRUE,sep = "")
edges2Q$from<- ""
edges2Q$total_code<- as.numeric(edges2Q$total_code)
edgesQ <- bind_rows(edges1Q, edges2Q)
vertices1Q<-read.table("/NASdata/wuwq/HFpEF/ALL/vertices1Q.txt", header=TRUE,sep = "")
vertices1Q$total_code<-as.numeric(vertices1Q$total_code)
vertices1Q$level<-as.numeric(vertices1Q$level)
vertices2Q<-read.table("/NASdata/wuwq/HFpEF/ALL/vertices2Q.txt", header=TRUE,sep = "")
vertices2Q$total_code<-as.numeric(vertices2Q$total_code)
vertices2Q$level<-as.numeric(vertices2Q$level)
vertices3Q <- tibble(
  node = '', language = NA, total_code = 0, level = 3
)
verticesQ <- bind_rows(vertices1Q, vertices2Q, vertices3Q) %>%
  mutate(
    radius = total_code**(1.8), # scaling circles
    language = factor(language, names(top_language_colorsQ))
  ) %>%
  arrange(level, language, node)

graphQ <- graph_from_data_frame(edgesQ, vertices = verticesQ)
# create custom layout by updating existing circle layout
layoutQ <- create_layout(graphQ, layout = 'circle')
outer_circleQ <- layoutQ %>%
  filter(level == 1) %>%
  mutate(language = factor(language, names(top_language_colorsQ))) %>%
  arrange(language, desc(name)) %>%
  mutate(
    x = cos((row_number() - 1) / 295 * 2 * pi),
    y = sin((row_number() - 1) / 295 * 2 * pi)
  )
# positioning circle centers manually by specifying polar coords
anglesQ <- c(30, 120, 210, 300, 0)
radiiQ <- c(0.5, 0.5, 0.5, 0.5, 0)
centersQ <- tibble(
  x = radiiQ * cos(anglesQ / 180 * pi),
  y = radiiQ * sin(anglesQ / 180 * pi)
)
inner_circleQ <- bind_cols(centersQ, select(filter(layoutQ, level != 1), -x, -y))
layoutQ[] <- bind_rows(outer_circleQ, inner_circleQ) %>% 
  arrange(.ggraph.index)
ggraph(layoutQ) +
  geom_edge_diagonal(
    aes(edge_color = node1.language, edge_alpha = as.factor(main)),
    edge_width = 0.2, show.legend = FALSE
  ) +
  geom_node_point(
    aes(size = radius*10000, color = language),
    alpha = 0.5, show.legend = FALSE
  ) +
  geom_node_text(
    aes(
      x = 1.0175* x,
      y = 1.0175 * y,
      label = name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = !(name %in% top_languagesQ)
    ),
    size = 2, hjust = 'outward', family = 'Oswald'
  ) +
  scale_edge_color_manual(values = top_language_colorsQ) +
  scale_color_manual(values = top_language_colorsQ) +
  scale_size_area(max_size = 70) +
  scale_edge_alpha_manual(values = c(0.6, 5)) +
  coord_fixed() +
  theme_void() +
  theme(
    text = element_text(family = 'Oswald'),
    legend.position = c(0.645, 0.51),
    plot.title = element_text(
      face = 'bold', hjust = 0.5, size = 5, margin = margin(t = 45, b = 3)
    ),
    plot.subtitle = element_text(
      face = 'plain', hjust = 0.5, size = 13, margin = margin(t = 5, b = 3)),
    plot.caption = element_text(
      face = 'plain', color = '#dedede', size = 8, hjust = 1,
      margin = margin(b = 20)
    )
  )


##### Fig 3G
N <- read.table("/NASdata/wuwq/HFpEF/ALL/heatmap.txt",  
                header=T, 
                row.names=1, 
                sep="\t") 
N<-as.matrix(N)
p <- pheatmap(N, scale="none",
              border="white", 
              cluster_cols = F, 
              cluster_rows = F,
              legend = T, 
              fontsize_row = F, 
              fontsize_col = 8,
              clustering_distance_rows = "euclidean", 
              clustering_method="centroid",
              fontsize_number = 5,
              color =  colorRampPalette(c("#4169e1","#87cefa","#e0ffff","white"))(100)) 
print(p)


##### Fig 3H
#main color of the circle
Rtop_languagesQ<-c(
  'Ppargc1a',
  'Rxrg',
  'Mitf',
  'Uqcrb'
)
#Rtop_language_colorsQ <- c(
#  '#efb306',
#  '#4e54ac',
#  '#0f8096',
#  '#17a769',
#  '#000000'
#)
Rtop_language_colorsQ <- c(
  '#e8351e',
  '#4e54ac',
  '#cd023d',
  '#17a769',
  '#000000'
)
# name of the 4 circles
names(Rtop_language_colorsQ) <- c(
  'Ppargc1a',
  'Rxrg',
  'Mitf',
  'Uqcrb',
  'All'
)
#lines from inner to outer
Redges1Q<-read.table("/NASdata/wuwq/HFpEF/ALL/edges1Q_HFr.txt", header=TRUE,sep = "")
Redges1Q$total_code<-as.numeric(Redges1Q$total_code)
Redges2Q<-read.table("/NASdata/wuwq/HFpEF/ALL/edges2Q_HFr.txt", header=TRUE,sep = "")
Redges2Q$from<- ""
Redges2Q$total_code<- as.numeric(Redges2Q$total_code)
RedgesQ <- bind_rows(Redges1Q, Redges2Q)
Rvertices1Q<-read.table("/NASdata/wuwq/HFpEF/ALL/vertices1Q_HFr.txt", header=TRUE,sep = "")
Rvertices1Q$total_code<-as.numeric(Rvertices1Q$total_code)
Rvertices1Q$level<-as.numeric(Rvertices1Q$level)
Rvertices2Q<-read.table("/NASdata/wuwq/HFpEF/ALL/vertices2Q_HFr.txt", header=TRUE,sep = "")
Rvertices2Q$total_code<-as.numeric(Rvertices2Q$total_code)
Rvertices2Q$level<-as.numeric(Rvertices2Q$level)
Rvertices3Q <- tibble(
  node = '', language = NA, total_code = 0, level = 3
)
RverticesQ <- bind_rows(Rvertices1Q, Rvertices2Q, Rvertices3Q) %>%
  mutate(
    radius = total_code**(1.8), # scaling circles
    language = factor(language, names(Rtop_language_colorsQ))
  ) %>%
  arrange(level, language, node)

graphQ <- graph_from_data_frame(RedgesQ, vertices = RverticesQ)


# create custom layout by updating existing circle layout
layoutQ <- create_layout(graphQ, layout = 'circle')
outer_circleQ <- layoutQ %>%
  filter(level == 1) %>%
  mutate(language = factor(language, names(Rtop_language_colorsQ))) %>%
  arrange(language, desc(name)) %>%
  mutate(
    x = cos((row_number() - 1) / 289 * 2 * pi),
    y = sin((row_number() - 1) / 289 * 2 * pi)
  )
# positioning circle centers manually by specifying polar coords
anglesQ <- c(10,100,190,280, 0)
radiiQ <- c(0.5, 0.5, 0.5, 0.5, 0)
centersQ <- tibble(
  x = radiiQ * cos(anglesQ / 180 * pi),
  y = radiiQ * sin(anglesQ / 180 * pi)
)
inner_circleQ <- bind_cols(centersQ,  select(filter(layoutQ, level != 1), -x, -y))
layoutQ[] <- bind_rows(outer_circleQ, inner_circleQ) %>% 
  arrange(.ggraph.index)
ggraph(layoutQ) +
  geom_edge_diagonal(
    aes(edge_color = node1.language, edge_alpha = as.factor(main)),
    edge_width = 0.2, show.legend = FALSE
  ) +
  geom_node_point(
    aes(size = radius*10000, color = language),
    alpha = 0.5, show.legend = FALSE
  ) +
  geom_node_text(
    aes(
      x = 1.0175* x,
      y = 1.0175 * y,
      label = name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = !(name %in% Rtop_languagesQ)
    ),
    size = 2, hjust = 'outward', family = 'Oswald'
  ) +
  scale_edge_color_manual(values = Rtop_language_colorsQ) +
  scale_color_manual(values = Rtop_language_colorsQ) +
  scale_size_area(max_size = 70) +
  scale_edge_alpha_manual(values = c(0.6, 5)) +
  coord_fixed() +
  theme_void() +
  theme(
    text = element_text(family = 'Oswald'),
    legend.position = c(0.645, 0.51),
    plot.title = element_text(
      face = 'bold', hjust = 0.5, size = 5, margin = margin(t = 45, b = 3)
    ),
    plot.subtitle = element_text(
      face = 'plain', hjust = 0.5, size = 13, margin = margin(t = 5, b = 3)),
    plot.caption = element_text(
      face = 'plain', color = '#dedede', size = 8, hjust = 1,
      margin = margin(b = 20)
    )
  )


##### Fig 3I
###Binary activity of single cells obtained based on SCENIC analysis, according to the barcode of CM1-8 subgroup, the binary activity of cells in CM3 subgroup is sorted in descending order in excel
### The figure was drawn in excel.


##### Fig S7-9 
Data of Fig S7-9 were downloded during above "SCENIC_associated analysis", and drawn in excel.


##### Fig S10
Data of Fig S10 were downloded during above "SCENIC_associated analysis".
N <- read.table("HFp_CM3_rss.txt",  
                header=T, 
                row.names=1, 
                sep="\t") 
N<-as.matrix(N)
p <- pheatmap(N, scale="none",
              border="black", 
              cluster_cols = F, 
              cluster_rows = F,
              legend = T, 
              legend_breaks=c(1,3,5),
              legend_labels = c(1,3,5),
              fontsize_row = 5, 
              fontsize_col = 8,
              clustering_distance_rows = "euclidean", 
              clustering_method="centroid",
              fontsize_number = 5,
              color =  colorRampPalette(c("white","pink","firebrick3"))(100))
print(p)

