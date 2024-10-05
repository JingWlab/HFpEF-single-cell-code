########## 心衰单细胞课题 ##########
########## SCENIC analysis
### load sc-RNA data
setwd("/")
S <- readRDS("Integrative_20230628.rds")
Idents(S) <- S@meta.data$stim
X <-subset(S, subset=stim!= "TAC_progress")
### correct annot_stim information
DefaultAssay(X)<- "RNA"
X$annot_stim <- paste(X$annot_corrected, X$stim, sep = "_")
unique(X$annot_stim)
#prepare cell meta information
setwd("/NASdata/wuwq/HFpEF/ALL")
cellInfo <- data.frame(X@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="annot_corrected")] <- "celltype"
colnames(cellInfo)[which(colnames(cellInfo)=="stim")] <- "model"
colnames(cellInfo)[which(colnames(cellInfo)=="annot_stim")] <- "modelcell"
saveRDS(cellInfo, file="int/cellInfo.Rds")
# prepare expression matrix
exprMat <- as.matrix(X@assays$RNA@counts)
saveRDS(exprMat, file="int/exprMat.Rds")
### Initialize settings
data(list="motifAnnotations_mgi_v9", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9
# ensure that cisTarget_databases folder has 2 files of 1G
scenicOptions <- initializeScenic(org="mgi", 
                                  dbDir="cisTarget_databases", nCores=10) 
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
### Co-expression network
exprMat_log <- log2(exprMat+1)
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
Correlation<-runCorrelation(exprMat_filtered, scenicOptions)
saveRDS(Correlation, file="int/Correlation.Rds") 
exprMat_filtered_log <- log2(exprMat_filtered+1) 
Genie3<-runGenie3(exprMat_filtered_log, scenicOptions)
saveRDS(Genie3, file="int/Genie3.Rds") 
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC")
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

### download the result
#AUC download
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
b<-aucell_regulonAUC@assays@data@listData$AUC
write.csv(b,file = "/NASdata/wuwq/HFpEF/ALL/auc_eachregulon.csv")
# binary regulon activity
minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
write.csv(binaryRegulonActivity,file="/NASdata/wuwq/HFpEF/ALL/onoff_eachcell_regulon.csv")
### RSS score to get cell-type specific regulon
regulonAUC<-loadInt(scenicOptions,"aucell_regulonAUC")
### After filtering RSS>0.01, 271 regulons were left
rss<-calcRSS(AUC=getAUC(regulonAUC),cellAnnotation = cellInfo[colnames(regulonAUC), "celltype"])
write.csv(rss,file="/NASdata/wuwq/HFpEF/ALL/celltypespec_reg.csv")
### 后续于excel进一步整理


##### Fig3A 
AA=readxl::read_xlsx("/NASdata/wuwq/HFpEF/ALL/bubble.xlsx",sheet = 4)
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


###于excel完成 Fig 3B-C 绘图

##### Fig 3D
S <- readRDS("Integrative_20230628.rds")
Idents(S) <- S@meta.data$stim
X <-subset(S, subset=stim!= "TAC_progress")
DefaultAssay(X) <- "RNA"
### 提取Ppargc1a在各个细胞类型的表达量信息
VlnPlot(X,
        features = c("Ppargc1a"),
        split.by = "annot_corrected",
        pt.size = 0)
###后续excel完成绘图
