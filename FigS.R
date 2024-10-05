##### Fig S1 基于PRISM 软件统计分析，于PPT中完成绘图

##### Fig S2A excel中完成绘图
##### Fig S2B https://satijalab.org/howmanycells 网站中完成绘图

##### Fig S2C
counts <- read.table("single_cell_genematrix.txt",sep = "\t")
counts <- counts[,-12577]
mat = matrix(0, nrow = 515, ncol = 12644)
rownames(mat)<- gene_info$ensamble

colnames(mat) <- nanogrid_nocontrol$sampleid
for (i in 1:515) {
  mat[i,]<- as.numeric(as.character(count1[gene_info$ensamble[i],]))
}

umi <- SingleCellExperiment(assays = list(counts = as.matrix(mat)),colData = nanogrid_nocontrol1)
altExp(umi,"ERCC") <- umi[grep("^ERCC-",rownames(umi)), ]
umi <- umi[grep("^ERCC-",rownames(umi),invert = T), ]
umi

library(readxl)
mito <- read_excel("mito.xlsx")
View(mito)                                                                                     
MT_names <- mito$gene_id
is_mito <- rownames(umi) %in% MT_names
table(is_mito)
umi_cell <- perCellQCMetrics(umi,subsets=list(Mito=is_mito))
head(umi_cell)[1:4,1:4]
umi_feature <- perFeatureQCMetrics(umi)
head(umi_feature)
umi <- addPerCellQC(umi, subsets=list(Mito=is_mito))
umi <- addPerFeatureQC(umi)
hist(
  umi$total,
  breaks = 100
)
qc.lib2 <- isOutlier(umi_cell$sum, 
                     nmads = 3,
                     log=TRUE, 
                     type="lower")
attr(qc.lib2, "thresholds")
abline(v = 25000, col = "red")


hist(
  umi_cell$detected,
  breaks = 100
)
qc.nexprs2 <- isOutlier(umi_cell$detected, 
                        nmads = 3,
                        log=TRUE,
                        type="lower")
attr(qc.nexprs2, "thresholds")
abline(v = 7000, col = "red")

qc.spike2 <- isOutlier(umi_cell$altexps_ERCC_percent, 
                       nmads = 3,
                       type="higher")
attr(qc.spike2, "thresholds")
qc.mito2 <- isOutlier(umi_cell$subsets_Mito_percent, 
                      nmads = 3,
                      type="higher")
attr(qc.mito2, "thresholds")
discard2 <- qc.lib2 | qc.nexprs2 | qc.spike2 | qc.mito2

DataFrame(LibSize=sum(qc.lib2), 
          NExprs=sum(qc.nexprs2), 
          SpikeProp=sum(qc.spike2), 
          MitoProp=sum(qc.mito2), 
          Total=sum(discard2))

reasons <- quickPerCellQC(umi_cell, 
                          sub.fields = c("subsets_Mito_percent", "altexps_ERCC_percent"))
colSums(as.matrix(reasons))
umi$discard <- reasons$discard
plotColData(umi, x="sum", y="subsets_Mito_percent", colour_by="discard")
plotColData(umi, x="sum", y="detected", colour_by="discard")
plotColData(umi, x="altexps_ERCC_percent", y="subsets_Mito_percent",colour_by="discard")
library(scales)
plotColData(umi, x="sum", y="detected", 
            colour_by="discard", other_fields = "Sample") + 
  facet_wrap(~Sample) + 
  scale_x_continuous(labels = unit_format(unit = "k", scale = 1e-3))

library(tidyverse)
library(SingleCellExperiment)
library(scater)

umi %>% 
  plotHighestExprs(exprs_values = "counts", 
                   feature_names_to_plot = "SYMBOL", 
                   colour_cells_by="detected")


##### Fig S2D
umi.qc.gene515 <- readRDS("~/umi.qc.gene515.rds")
plotRLE(umi.qc.gene515, exprs_values = "logcounts_raw",colour_by = "Stim") + ggtitle("RLE plot for log2(CPM) counts")

##### Fig S3 原始数据在supplyment table里，于于excel中完成绘图


##### Fig S4 ###承接Fig 2A-B
pathway.union <- union(HFrEF_Con@netP$pathways, HFrEF_HF@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(HFrEF_Con, pattern = "outgoing", signaling = pathway.union, title = "Control", width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(HFrEF_HF, pattern = "outgoing", signaling = pathway.union, title = "HFrEF", width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

HFrEF_HF <- netAnalysis_computeCentrality(HFrEF_HF, slot.name = "netP")
netAnalysis_signalingRole_network(HFrEF_HF, signaling = "COLLAGEN", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(HFrEF_HF, signaling = "FN1", width = 8, height = 2.5, font.size = 10)

##### Fig S5
pathway.union <- union(HFpEF_Con@netP$pathways, HFpEF_HF@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(HFpEF_Con, pattern = "outgoing", signaling = pathway.union, title = "Control", width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(HFpEF_HF, pattern = "outgoing", signaling = pathway.union, title = "HFrEF", width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

HFpEF_HF <- netAnalysis_computeCentrality(HFpEF_HF, slot.name = "netP")
netAnalysis_signalingRole_network(HFpEF_HF, signaling = "VCAM", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(HFpEF_HF, signaling = "HSPG", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(HFpEF_HF, signaling = "EPHA", width = 8, height = 2.5, font.size = 10)

##### Fig S6-8 基于SCENIC分析得到的每个细胞的binary activity,于excel中完成绘图


##### Fig S9
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


##### Fig S10A
#HFpEF_relavent
N <- read.table("/NASdata/wuwq/HFpEF/GWAS/DEGs_new/HFp-relevant.txt",  
                header=T, 
                row.names=1, 
                sep="\t") 
N<-as.matrix(N)
p <- pheatmap(N, scale="none",
              border="white",
              cluster_cols = F, 
              cluster_rows = F,
              legend = T, 
              legend_breaks=c(-2,-1,-0,1,2,3,4,5,6),
              legend_labels = c(-2,-1,-0,1,2,3,4,5,6),
              fontsize_row = 8, 
              fontsize_col = 8,
              clustering_distance_rows = "euclidean", 
              clustering_method="centroid",
              fontsize_number = 5,
              color =  colorRampPalette(c("#6495ED","#87CEFF","white","#FFC0CB","#F08080","#E34234","#D22B2B","#C41E3A","#800000"))(100))
print(p)


##### Fig S10B
M <- read.table("/NASdata/wuwq/HFpEF/GWAS/DEGs_new/TAC-relevant.txt",  
                header=T, 
                row.names=1, 
                sep="\t") 
M<-as.matrix(M)
p2 <- pheatmap(M, scale="none",
               border="white", 
               cluster_cols = F,
               cluster_rows = F,
               legend = T, 
               legend_breaks=c(-5,-4,-3,-2,-1,0,1,2),
               legend_labels = c(-5,-4,-3,-2,-1,0,1,2),
               fontsize_row = 8, 
               fontsize_col = 8,
               clustering_distance_rows = "euclidean", 
               clustering_method="centroid",
               fontsize_number = 5,
               color =  colorRampPalette(c("#000080","#104E8B","#4169E1","#6495ED","#87CEFF","#FFFFFF","#F08080","#CD4F39"))(100))
print(p2)


##### Fig S11
###load heart_altlas dataset
library(SeuratDisk)
Convert("/NASdata/wuwq/HFpEF/GWAS/global_raw.h5ad", dest="h5seurat",
        assay = "RNA",
        overwrite=F)
heart_altlas <- LoadH5Seurat('/data/wuwq/Figure_wwq/global_raw.h5seurat')
H = UpdateSeuratObject(heart_altlas)
unique(H@meta.data$type)
unique(H@meta.data$age_group)
unique(H@meta.data$cell_type)
### 获取FKBP7,RGS10在Ventricular_Cardiomyocyte中的表达量
DefaultAssay(H) <- "RNA"
a <- subset (H, subset= cell_type=="Ventricular_Cardiomyocyte")
FKBP7<- FetchData(a,vars = "FKBP7")
write.csv(FKBP7, file = "/data/wuwq/Figure_wwq/heart_altlas/FKBP7_a")
RGS10<- FetchData(a,vars = "RGS10")
write.csv(RGS10, file = "/data/wuwq/Figure_wwq/heart_altlas/RGS10_a")
###代码同上，依次获取FKBP7,RGS10在Atrial_Cardiomyocyte,Endothelial,Fibroblast,Lymphoid,Myeloid中的表达量
### 于excel完成后续作图
