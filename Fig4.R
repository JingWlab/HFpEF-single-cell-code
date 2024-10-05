########## GWAS & eQTL integration analysis
##### Fig 4A-D
###基于excel整理和在线网站LDLink的LD连锁不平衡计算，于ppt绘制完成


##### Fig 4E
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
### 获取ACTN2,HSPB7,TTN,FKBP7,HP,RGS10,YWHAE在Ventricular_Cardiomyocyte中的表达量
DefaultAssay(H) <- "RNA"
a <- subset (H, subset= cell_type=="Ventricular_Cardiomyocyte")
ACTN2<- FetchData(a,vars = "ACTN2")
write.csv(ACTN2, file = "/data/wuwq/Figure_wwq/heart_altlas/ACTN2_a")
HSPB7<- FetchData(a,vars = "HSPB7")
write.csv(HSPB7, file = "/data/wuwq/Figure_wwq/heart_altlas/HSPB7_a")
TTN<- FetchData(a,vars = "TTN")
write.csv(TTN, file = "/data/wuwq/Figure_wwq/heart_altlas/TTN_a")
FKBP7<- FetchData(a,vars = "FKBP7")
write.csv(FKBP7, file = "/data/wuwq/Figure_wwq/heart_altlas/FKBP7_a")
HP<- FetchData(a,vars = "HP")
write.csv(HP, file = "/data/wuwq/Figure_wwq/heart_altlas/HP_a")
RGS10<- FetchData(a,vars = "RGS10")
write.csv(RGS10, file = "/data/wuwq/Figure_wwq/heart_altlas/RGS10_a")
YWHAE<- FetchData(a,vars = "YWHAE")
write.csv(YWHAE, file = "/data/wuwq/Figure_wwq/heart_altlas/YWHAE_a")
###代码同上，依次获取ACTN2,HSPB7,TTN,FKBP7,HP,RGS10,YWHAE在Atrial_Cardiomyocyte,Endothelial,Fibroblast,Lymphoid,Myeloid中的表达量
### 于excel完成后续作图


##### Fig 4F
###作图依据：colocalization_GWAS_eQTL coloc共定位 (作图于excel中完成)
### E*_1~6 代表6个心脏相关组织，包括：Artery-Aorta, artery cornarory, heart atrial appendage, heart left ventricle, thyroid and whole blood
### 匹配MAF及得到共定位的PPH0-PPH4结果
#E1_1~E1_6  MAF
MAF_E1 <- read.table(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/MAF_E1.txt", fill=TRUE, header = T, as.is = T)
head(MAF_E1)
MAF_E1=MAF_E1[ , c(1,4)]
colnames(MAF_E1) <-c("rs_id", "maf")
head(MAF_E1)

#gwas
gwas_E1 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/gwaschr1_rs113151268.bed", fill=TRUE, header = T, as.is = T)

#E1_1
E1_1_HSPB7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E1_1_HSPB7.bed", fill=TRUE, header = T, as.is = T)
E1_1_HSPB7=E1_1_HSPB7[ , c(4,6,10)]
colnames(E1_1_HSPB7) <-c("rs_id", "BP","pval_nominal")
fE1_1_HSPB7 <- merge(MAF_E1,E1_1_HSPB7, by = "rs_id")
head(fE1_1_HSPB7)
dim(fE1_1_HSPB7)
write.table(fE1_1_HSPB7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE1_1_HSPB7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E1_1 <- merge(fE1_1_HSPB7, gwas_E1, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E1_1)
#coloc 
input_E1_1 <- input_E1_1[complete.cases(input_E1_1),]
coloc_E1_1 <- coloc.abf(
  dataset1=list(snp=input_E1_1$rs_id,pvalues=input_E1_1$P, type="cc", s=0.51, N=nrow(gwas_E1)),
  dataset2=list(snp=input_E1_1$rs_id,pvalues=input_E1_1$pval_nominal, type="quant", N=nrow(fE1_1_HSPB7)),
  MAF=input_E1_1$maf)

#E1_2
E1_2_HSPB7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E1_2_HSPB7.bed", fill=TRUE, header = T, as.is = T)
E1_2_HSPB7=E1_2_HSPB7[ , c(4,6,10)]
colnames(E1_2_HSPB7) <-c("rs_id", "BP","pval_nominal")
fE1_2_HSPB7 <- merge(MAF_E1,E1_2_HSPB7, by = "rs_id")
head(fE1_2_HSPB7)
dim(fE1_2_HSPB7)
write.table(fE1_2_HSPB7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE1_2_HSPB7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E1_2 <- merge(fE1_2_HSPB7, gwas_E1, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E1_2)
#coloc 
input_E1_2 <- input_E1_2[complete.cases(input_E1_2),]
coloc_E1_2 <- coloc.abf(
  dataset1=list(snp=input_E1_2$rs_id,pvalues=input_E1_2$P, type="cc", s=0.51, N=nrow(gwas_E1)),
  dataset2=list(snp=input_E1_2$rs_id,pvalues=input_E1_2$pval_nominal, type="quant", N=nrow(fE1_2_HSPB7)),
  MAF=input_E1_2$maf)

#E1_3
E1_3_HSPB7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E1_3_HSPB7.bed", fill=TRUE, header = T, as.is = T)
E1_3_HSPB7=E1_3_HSPB7[ , c(4,6,10)]
colnames(E1_3_HSPB7) <-c("rs_id", "BP","pval_nominal")
fE1_3_HSPB7 <- merge(MAF_E1,E1_3_HSPB7, by = "rs_id")
head(fE1_3_HSPB7)
dim(fE1_3_HSPB7)
write.table(fE1_3_HSPB7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE1_3_HSPB7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E1_3 <- merge(fE1_3_HSPB7, gwas_E1, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E1_3)
#coloc 
input_E1_3 <- input_E1_3[complete.cases(input_E1_3),]
coloc_E1_3 <- coloc.abf(
  dataset1=list(snp=input_E1_3$rs_id,pvalues=input_E1_3$P, type="cc", s=0.51, N=nrow(gwas_E1)),
  dataset2=list(snp=input_E1_3$rs_id,pvalues=input_E1_3$pval_nominal, type="quant", N=nrow(fE1_3_HSPB7)),
  MAF=input_E1_3$maf)

#E1_4
E1_4_HSPB7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E1_4_HSPB7.bed", fill=TRUE, header = T, as.is = T)
E1_4_HSPB7=E1_4_HSPB7[ , c(4,6,10)]
colnames(E1_4_HSPB7) <-c("rs_id", "BP","pval_nominal")
fE1_4_HSPB7 <- merge(MAF_E1,E1_4_HSPB7, by = "rs_id")
head(fE1_4_HSPB7)
dim(fE1_4_HSPB7)
write.table(fE1_4_HSPB7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE1_4_HSPB7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E1_4 <- merge(fE1_4_HSPB7, gwas_E1, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E1_4)
#coloc 
input_E1_4 <- input_E1_4[complete.cases(input_E1_4),]
coloc_E1_4 <- coloc.abf(
  dataset1=list(snp=input_E1_4$rs_id,pvalues=input_E1_4$P, type="cc", s=0.51, N=nrow(gwas_E1)),
  dataset2=list(snp=input_E1_4$rs_id,pvalues=input_E1_4$pval_nominal, type="quant", N=nrow(fE1_4_HSPB7)),
  MAF=input_E1_4$maf)

#E1_5
E1_5_HSPB7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E1_5_HSPB7.bed", fill=TRUE, header = T, as.is = T)
E1_5_HSPB7=E1_5_HSPB7[ , c(4,6,10)]
colnames(E1_5_HSPB7) <-c("rs_id", "BP","pval_nominal")
fE1_5_HSPB7 <- merge(MAF_E1,E1_5_HSPB7, by = "rs_id")
head(fE1_5_HSPB7)
dim(fE1_5_HSPB7)
write.table(fE1_5_HSPB7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE1_5_HSPB7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E1_5 <- merge(fE1_5_HSPB7, gwas_E1, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E1_5)
#coloc 
input_E1_5 <- input_E1_5[complete.cases(input_E1_5),]
coloc_E1_5 <- coloc.abf(
  dataset1=list(snp=input_E1_5$rs_id,pvalues=input_E1_5$P, type="cc", s=0.51, N=nrow(gwas_E1)),
  dataset2=list(snp=input_E1_5$rs_id,pvalues=input_E1_5$pval_nominal, type="quant", N=nrow(fE1_5_HSPB7)),
  MAF=input_E1_5$maf)

#E1_6
E1_6_HSPB7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E1_6_HSPB7.bed", fill=TRUE, header = T, as.is = T)
E1_6_HSPB7=E1_6_HSPB7[ , c(4,6,10)]
colnames(E1_6_HSPB7) <-c("rs_id", "BP","pval_nominal")
fE1_6_HSPB7 <- merge(MAF_E1,E1_6_HSPB7, by = "rs_id")
head(fE1_6_HSPB7)
dim(fE1_6_HSPB7)
write.table(fE1_6_HSPB7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE1_6_HSPB7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E1_6 <- merge(fE1_6_HSPB7, gwas_E1, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E1_6)
#coloc 
input_E1_6 <- input_E1_6[complete.cases(input_E1_6),]
coloc_E1_6 <- coloc.abf(
  dataset1=list(snp=input_E1_6$rs_id,pvalues=input_E1_6$P, type="cc", s=0.51, N=nrow(gwas_E1)),
  dataset2=list(snp=input_E1_6$rs_id,pvalues=input_E1_6$pval_nominal, type="quant", N=nrow(fE1_6_HSPB7)),
  MAF=input_E1_6$maf)


rm(list = ls())
###E2_1~E2_6  MAF
MAF_E2 <- read.table(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/MAF_E2.txt", fill=TRUE, header = T, as.is = T)
head(MAF_E2)
MAF_E2=MAF_E2[ , c(1,4)]
colnames(MAF_E2) <-c("rs_id", "maf")
head(MAF_E2)

#gwas
gwas_E2 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/norgwaschr1_rs28579893.bed", fill=TRUE, header = T, as.is = T)

#E2_1
E2_1_HSPB7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E2_1_HSPB7.bed", fill=TRUE, header = T, as.is = T)
E2_1_HSPB7=E2_1_HSPB7[ , c(4,6,10)]
colnames(E2_1_HSPB7) <-c("rs_id", "BP","pval_nominal")
fE2_1_HSPB7 <- merge(MAF_E2,E2_1_HSPB7, by = "rs_id")
head(fE2_1_HSPB7)
dim(fE2_1_HSPB7)
write.table(fE2_1_HSPB7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE2_1_HSPB7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E2_1 <- merge(fE2_1_HSPB7, gwas_E2, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E2_1)
#coloc 
input_E2_1 <- input_E2_1[complete.cases(input_E2_1),]
coloc_E2_1 <- coloc.abf(
  dataset1=list(snp=input_E2_1$rs_id,pvalues=input_E2_1$P, type="cc", s=0.51, N=nrow(gwas_E2)),
  dataset2=list(snp=input_E2_1$rs_id,pvalues=input_E2_1$pval_nominal, type="quant", N=nrow(fE2_1_HSPB7)),
  MAF=input_E2_1$maf)

#E2_2
E2_2_HSPB7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E2_2_HSPB7.bed", fill=TRUE, header = T, as.is = T)
E2_2_HSPB7=E2_2_HSPB7[ , c(4,6,10)]
colnames(E2_2_HSPB7) <-c("rs_id", "BP","pval_nominal")
fE2_2_HSPB7 <- merge(MAF_E2,E2_2_HSPB7, by = "rs_id")
head(fE2_2_HSPB7)
dim(fE2_2_HSPB7)
write.table(fE2_2_HSPB7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE2_2_HSPB7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E2_2 <- merge(fE2_2_HSPB7, gwas_E2, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E2_2)
#coloc 
input_E2_2 <- input_E2_2[complete.cases(input_E2_2),]
coloc_E2_2 <- coloc.abf(
  dataset1=list(snp=input_E2_2$rs_id,pvalues=input_E2_2$P, type="cc", s=0.51, N=nrow(gwas_E2)),
  dataset2=list(snp=input_E2_2$rs_id,pvalues=input_E2_2$pval_nominal, type="quant", N=nrow(fE2_2_HSPB7)),
  MAF=input_E2_2$maf)

#E2_3
E2_3_HSPB7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E2_3_HSPB7.bed", fill=TRUE, header = T, as.is = T)
E2_3_HSPB7=E2_3_HSPB7[ , c(4,6,10)]
colnames(E2_3_HSPB7) <-c("rs_id", "BP","pval_nominal")
fE2_3_HSPB7 <- merge(MAF_E2,E2_3_HSPB7, by = "rs_id")
head(fE2_3_HSPB7)
dim(fE2_3_HSPB7)
write.table(fE2_3_HSPB7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE2_3_HSPB7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E2_3 <- merge(fE2_3_HSPB7, gwas_E2, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E2_3)
#coloc 
input_E2_3 <- input_E2_3[complete.cases(input_E2_3),]
coloc_E2_3 <- coloc.abf(
  dataset1=list(snp=input_E2_3$rs_id,pvalues=input_E2_3$P, type="cc", s=0.51, N=nrow(gwas_E2)),
  dataset2=list(snp=input_E2_3$rs_id,pvalues=input_E2_3$pval_nominal, type="quant", N=nrow(fE2_3_HSPB7)),
  MAF=input_E2_3$maf)

#E2_4
E2_4_HSPB7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E2_4_HSPB7.bed", fill=TRUE, header = T, as.is = T)
E2_4_HSPB7=E2_4_HSPB7[ , c(4,6,10)]
colnames(E2_4_HSPB7) <-c("rs_id", "BP","pval_nominal")
fE2_4_HSPB7 <- merge(MAF_E2,E2_4_HSPB7, by = "rs_id")
head(fE2_4_HSPB7)
dim(fE2_4_HSPB7)
write.table(fE2_4_HSPB7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE2_4_HSPB7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E2_4 <- merge(fE2_4_HSPB7, gwas_E2, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E2_4)
#coloc 
input_E2_4 <- input_E2_4[complete.cases(input_E2_4),]
coloc_E2_4 <- coloc.abf(
  dataset1=list(snp=input_E2_4$rs_id,pvalues=input_E2_4$P, type="cc", s=0.51, N=nrow(gwas_E2)),
  dataset2=list(snp=input_E2_4$rs_id,pvalues=input_E2_4$pval_nominal, type="quant", N=nrow(fE2_4_HSPB7)),
  MAF=input_E2_4$maf)

#E2_5
E2_5_HSPB7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E2_5_HSPB7.bed", fill=TRUE, header = T, as.is = T)
E2_5_HSPB7=E2_5_HSPB7[ , c(4,6,10)]
colnames(E2_5_HSPB7) <-c("rs_id", "BP","pval_nominal")
fE2_5_HSPB7 <- merge(MAF_E2,E2_5_HSPB7, by = "rs_id")
head(fE2_5_HSPB7)
dim(fE2_5_HSPB7)
write.table(fE2_5_HSPB7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE2_5_HSPB7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E2_5 <- merge(fE2_5_HSPB7, gwas_E2, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E2_5)
#coloc 
input_E2_5 <- input_E2_5[complete.cases(input_E2_5),]
coloc_E2_5 <- coloc.abf(
  dataset1=list(snp=input_E2_5$rs_id,pvalues=input_E2_5$P, type="cc", s=0.51, N=nrow(gwas_E2)),
  dataset2=list(snp=input_E2_5$rs_id,pvalues=input_E2_5$pval_nominal, type="quant", N=nrow(fE2_5_HSPB7)),
  MAF=input_E2_5$maf)

#E2_6
E2_6_HSPB7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E2_6_HSPB7.bed", fill=TRUE, header = T, as.is = T)
E2_6_HSPB7=E2_6_HSPB7[ , c(4,6,10)]
colnames(E2_6_HSPB7) <-c("rs_id", "BP","pval_nominal")
fE2_6_HSPB7 <- merge(MAF_E2,E2_6_HSPB7, by = "rs_id")
head(fE2_6_HSPB7)
dim(fE2_6_HSPB7)
write.table(fE2_6_HSPB7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE2_6_HSPB7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E2_6 <- merge(fE2_6_HSPB7, gwas_E2, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E2_6)
#coloc 
input_E2_6 <- input_E2_6[complete.cases(input_E2_6),]
coloc_E2_6 <- coloc.abf(
  dataset1=list(snp=input_E2_6$rs_id,pvalues=input_E2_6$P, type="cc", s=0.51, N=nrow(gwas_E2)),
  dataset2=list(snp=input_E2_6$rs_id,pvalues=input_E2_6$pval_nominal, type="quant", N=nrow(fE2_6_HSPB7)),
  MAF=input_E2_6$maf)


rm(list = ls())
###E3_1~E3_6  MAF
MAF_E3 <- read.table(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/MAF_E3.txt", fill=TRUE, header = T, as.is = T)
head(MAF_E3)
MAF_E3=MAF_E3[ , c(1,4)]
colnames(MAF_E3) <-c("rs_id", "maf")
head(MAF_E3)

#gwas
gwas_E3 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/MTAGgwaschr6_rs6915002.bed", fill=TRUE, header = T, as.is = T)

#E3_1
E3_1_MLIP <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E3_1_MLIP.bed", fill=TRUE, header = T, as.is = T)
E3_1_MLIP=E3_1_MLIP[ , c(4,6,10)]
colnames(E3_1_MLIP) <-c("rs_id", "BP","pval_nominal")

fE3_1_MLIP <- merge(MAF_E3,E3_1_MLIP, by = "rs_id")
head(fE3_1_MLIP)
dim(fE3_1_MLIP)
write.table(fE3_1_MLIP,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE3_1_MLIP.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E3_1 <- merge(fE3_1_MLIP, gwas_E3, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E3_1)
#coloc 
input_E3_1 <- input_E3_1[complete.cases(input_E3_1),]
coloc_E3_1 <- coloc.abf(
  dataset1=list(snp=input_E3_1$rs_id,pvalues=input_E3_1$P, type="cc", s=0.51, N=nrow(gwas_E3)),
  dataset2=list(snp=input_E3_1$rs_id,pvalues=input_E3_1$pval_nominal, type="quant", N=nrow(fE3_1_MLIP)),
  MAF=input_E3_1$maf)

#E3_2
E3_2_MLIP <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E3_2_MLIP.bed", fill=TRUE, header = T, as.is = T)
E3_2_MLIP=E3_2_MLIP[ , c(4,6,10)]
colnames(E3_2_MLIP) <-c("rs_id", "BP","pval_nominal")
fE3_2_MLIP <- merge(MAF_E3,E3_2_MLIP, by = "rs_id")
head(fE3_2_MLIP)
dim(fE3_2_MLIP)
write.table(fE3_2_MLIP,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE3_2_MLIP.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E3_2 <- merge(fE3_2_MLIP, gwas_E3, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E3_2)
#coloc 
input_E3_2 <- input_E3_2[complete.cases(input_E3_2),]
coloc_E3_2 <- coloc.abf(
  dataset1=list(snp=input_E3_2$rs_id,pvalues=input_E3_2$P, type="cc", s=0.51, N=nrow(gwas_E3)),
  dataset2=list(snp=input_E3_2$rs_id,pvalues=input_E3_2$pval_nominal, type="quant", N=nrow(fE3_2_MLIP)),
  MAF=input_E3_2$maf)

#E3_3
E3_3_MLIP <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E3_3_MLIP.bed", fill=TRUE, header = T, as.is = T)
E3_3_MLIP=E3_3_MLIP[ , c(4,6,10)]
colnames(E3_3_MLIP) <-c("rs_id", "BP","pval_nominal")
fE3_3_MLIP <- merge(MAF_E3,E3_3_MLIP, by = "rs_id")
head(fE3_3_MLIP)
dim(fE3_3_MLIP)
write.table(fE3_3_MLIP,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE3_3_MLIP.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E3_3 <- merge(fE3_3_MLIP, gwas_E3, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E3_3)
#coloc 
input_E3_3 <- input_E3_3[complete.cases(input_E3_3),]
coloc_E3_3 <- coloc.abf(
  dataset1=list(snp=input_E3_3$rs_id,pvalues=input_E3_3$P, type="cc", s=0.51, N=nrow(gwas_E3)),
  dataset2=list(snp=input_E3_3$rs_id,pvalues=input_E3_3$pval_nominal, type="quant", N=nrow(fE3_3_MLIP)),
  MAF=input_E3_3$maf)

#E3_4
E3_4_MLIP <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E3_4_MLIP.bed", fill=TRUE, header = T, as.is = T)
E3_4_MLIP=E3_4_MLIP[ , c(4,6,10)]
colnames(E3_4_MLIP) <-c("rs_id", "BP","pval_nominal")
fE3_4_MLIP <- merge(MAF_E3,E3_4_MLIP, by = "rs_id")
head(fE3_4_MLIP)
dim(fE3_4_MLIP)
write.table(fE3_4_MLIP,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE3_4_MLIP.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E3_4 <- merge(fE3_4_MLIP, gwas_E3, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E3_4)
#coloc 
input_E3_4 <- input_E3_4[complete.cases(input_E3_4),]
coloc_E3_4 <- coloc.abf(
  dataset1=list(snp=input_E3_4$rs_id,pvalues=input_E3_4$P, type="cc", s=0.51, N=nrow(gwas_E3)),
  dataset2=list(snp=input_E3_4$rs_id,pvalues=input_E3_4$pval_nominal, type="quant", N=nrow(fE3_4_MLIP)),
  MAF=input_E3_4$maf)

#E3_5
E3_5_MLIP <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E3_5_MLIP.bed", fill=TRUE, header = T, as.is = T)
E3_5_MLIP=E3_5_MLIP[ , c(4,6,10)]
colnames(E3_5_MLIP) <-c("rs_id", "BP","pval_nominal")
fE3_5_MLIP <- merge(MAF_E3,E3_5_MLIP, by = "rs_id")
head(fE3_5_MLIP)
dim(fE3_5_MLIP)
write.table(fE3_5_MLIP,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE3_5_MLIP.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E3_5 <- merge(fE3_5_MLIP, gwas_E3, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E3_5)
#coloc 
input_E3_5 <- input_E3_5[complete.cases(input_E3_5),]
coloc_E3_5 <- coloc.abf(
  dataset1=list(snp=input_E3_5$rs_id,pvalues=input_E3_5$P, type="cc", s=0.51, N=nrow(gwas_E3)),
  dataset2=list(snp=input_E3_5$rs_id,pvalues=input_E3_5$pval_nominal, type="quant", N=nrow(fE3_5_MLIP)),
  MAF=input_E3_5$maf)

#E3_6
E3_6_MLIP <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E3_6_MLIP.bed", fill=TRUE, header = T, as.is = T)
E3_6_MLIP=E3_6_MLIP[ , c(4,6,10)]
colnames(E3_6_MLIP) <-c("rs_id", "BP","pval_nominal")
fE3_6_MLIP <- merge(MAF_E3,E3_6_MLIP, by = "rs_id")
head(fE3_6_MLIP)
dim(fE3_6_MLIP)
write.table(fE3_6_MLIP,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE3_6_MLIP.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E3_6 <- merge(fE3_6_MLIP, gwas_E3, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E3_6)
#coloc 
input_E3_6 <- input_E3_6[complete.cases(input_E3_6),]
coloc_E3_6 <- coloc.abf(
  dataset1=list(snp=input_E3_6$rs_id,pvalues=input_E3_6$P, type="cc", s=0.51, N=nrow(gwas_E3)),
  dataset2=list(snp=input_E3_6$rs_id,pvalues=input_E3_6$pval_nominal, type="quant", N=nrow(fE3_6_MLIP)),
  MAF=input_E3_6$maf)


rm(list = ls())
###E4E5_1~E4E5_6  MAF
MAF_E4E5 <- read.table(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/MAF_E4E5.txt", fill=TRUE, header = T, as.is = T)
head(MAF_E4E5)
MAF_E4E5=MAF_E4E5[ , c(1,4)]
colnames(MAF_E4E5) <-c("rs_id", "maf")
head(MAF_E4E5)

#gwas
gwas_E4E5 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/MTAGgwaschr1_rs12724121.bed", fill=TRUE, header = T, as.is = T)

#E4_1
E4_1_ACTN2 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E4_1_ACTN2.bed", fill=TRUE, header = T, as.is = T)
E4_1_ACTN2=E4_1_ACTN2[ , c(4,6,10)]
colnames(E4_1_ACTN2) <-c("rs_id", "BP","pval_nominal")
fE4_1_ACTN2 <- merge(MAF_E4E5,E4_1_ACTN2, by = "rs_id")
head(fE4_1_ACTN2)
dim(fE4_1_ACTN2)
write.table(fE4_1_ACTN2,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE4_1_ACTN2.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E4_1 <- merge(fE4_1_ACTN2, gwas_E4E5, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E4_1)
#coloc 
input_E4_1 <- input_E4_1[complete.cases(input_E4_1),]
coloc_E4_1 <- coloc.abf(
  dataset1=list(snp=input_E4_1$rs_id,pvalues=input_E4_1$P, type="cc", s=0.51, N=nrow(gwas_E4E5)),
  dataset2=list(snp=input_E4_1$rs_id,pvalues=input_E4_1$pval_nominal, type="quant", N=nrow(fE4_1_ACTN2)),
  MAF=input_E4_1$maf)

#E4_2
E4_2_ACTN2 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E4_2_ACTN2.bed", fill=TRUE, header = T, as.is = T)
E4_2_ACTN2=E4_2_ACTN2[ , c(4,6,10)]
colnames(E4_2_ACTN2) <-c("rs_id", "BP","pval_nominal")
fE4_2_ACTN2 <- merge(MAF_E4E5,E4_2_ACTN2, by = "rs_id")
head(fE4_2_ACTN2)
dim(fE4_2_ACTN2)
write.table(fE4_2_ACTN2,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE4_2_ACTN2.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E4_2 <- merge(fE4_2_ACTN2, gwas_E4E5, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E4_2)
#coloc 
input_E4_2 <- input_E4_2[complete.cases(input_E4_2),]
coloc_E4_2 <- coloc.abf(
  dataset1=list(snp=input_E4_2$rs_id,pvalues=input_E4_2$P, type="cc", s=0.51, N=nrow(gwas_E4E5)),
  dataset2=list(snp=input_E4_2$rs_id,pvalues=input_E4_2$pval_nominal, type="quant", N=nrow(fE4_2_ACTN2)),
  MAF=input_E4_2$maf)

#E4_3
E4_3_ACTN2 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E4_3_ACTN2.bed", fill=TRUE, header = T, as.is = T)
E4_3_ACTN2=E4_3_ACTN2[ , c(4,6,10)]
colnames(E4_3_ACTN2) <-c("rs_id", "BP","pval_nominal")
fE4_3_ACTN2 <- merge(MAF_E4E5,E4_3_ACTN2, by = "rs_id")
head(fE4_3_ACTN2)
dim(fE4_3_ACTN2)
write.table(fE4_3_ACTN2,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE4_3_ACTN2.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E4_3 <- merge(fE4_3_ACTN2, gwas_E4E5, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E4_3)
#coloc 
input_E4_3 <- input_E4_3[complete.cases(input_E4_3),]
coloc_E4_3 <- coloc.abf(
  dataset1=list(snp=input_E4_3$rs_id,pvalues=input_E4_3$P, type="cc", s=0.51, N=nrow(gwas_E4E5)),
  dataset2=list(snp=input_E4_3$rs_id,pvalues=input_E4_3$pval_nominal, type="quant", N=nrow(fE4_3_ACTN2)),
  MAF=input_E4_3$maf)

#E4_4
E4_4_ACTN2 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E4_4_ACTN2.bed", fill=TRUE, header = T, as.is = T)
E4_4_ACTN2=E4_4_ACTN2[ , c(4,6,10)]
colnames(E4_4_ACTN2) <-c("rs_id", "BP","pval_nominal")
fE4_4_ACTN2 <- merge(MAF_E4E5,E4_4_ACTN2, by = "rs_id")
head(fE4_4_ACTN2)
dim(fE4_4_ACTN2)
write.table(fE4_4_ACTN2,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE4_4_ACTN2.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E4_4 <- merge(fE4_4_ACTN2, gwas_E4E5, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E4_4)
#coloc 
input_E4_4 <- input_E4_4[complete.cases(input_E4_4),]
coloc_E4_4 <- coloc.abf(
  dataset1=list(snp=input_E4_4$rs_id,pvalues=input_E4_4$P, type="cc", s=0.51, N=nrow(gwas_E4E5)),
  dataset2=list(snp=input_E4_4$rs_id,pvalues=input_E4_4$pval_nominal, type="quant", N=nrow(fE4_4_ACTN2)),
  MAF=input_E4_4$maf)

#E4_5
E4_5_ACTN2 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E4_5_ACTN2.bed", fill=TRUE, header = T, as.is = T)
E4_5_ACTN2=E4_5_ACTN2[ , c(4,6,10)]
colnames(E4_5_ACTN2) <-c("rs_id", "BP","pval_nominal")
fE4_5_ACTN2 <- merge(MAF_E4E5,E4_5_ACTN2, by = "rs_id")
head(fE4_5_ACTN2)
dim(fE4_5_ACTN2)
write.table(fE4_5_ACTN2,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE4_5_ACTN2.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E4_5 <- merge(fE4_5_ACTN2, gwas_E4E5, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E4_5)
#coloc 
input_E4_5 <- input_E4_5[complete.cases(input_E4_5),]
coloc_E4_5 <- coloc.abf(
  dataset1=list(snp=input_E4_5$rs_id,pvalues=input_E4_5$P, type="cc", s=0.51, N=nrow(gwas_E4E5)),
  dataset2=list(snp=input_E4_5$rs_id,pvalues=input_E4_5$pval_nominal, type="quant", N=nrow(fE4_5_ACTN2)),
  MAF=input_E4_5$maf)

#E4_6
E4_6_ACTN2 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E4_6_ACTN2.bed", fill=TRUE, header = T, as.is = T)
E4_6_ACTN2=E4_6_ACTN2[ , c(4,6,10)]
colnames(E4_6_ACTN2) <-c("rs_id", "BP","pval_nominal")
fE4_6_ACTN2 <- merge(MAF_E4E5,E4_6_ACTN2, by = "rs_id")
head(fE4_6_ACTN2)
dim(fE4_6_ACTN2)
write.table(fE4_6_ACTN2,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE4_6_ACTN2.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E4_6 <- merge(fE4_6_ACTN2, gwas_E4E5, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E4_6)
#coloc 
input_E4_6 <- input_E4_6[complete.cases(input_E4_6),]
coloc_E4_6 <- coloc.abf(
  dataset1=list(snp=input_E4_6$rs_id,pvalues=input_E4_6$P, type="cc", s=0.51, N=nrow(gwas_E4E5)),
  dataset2=list(snp=input_E4_6$rs_id,pvalues=input_E4_6$pval_nominal, type="quant", N=nrow(fE4_6_ACTN2)),
  MAF=input_E4_6$maf)

#E5_1 ~ 6 have no MLF1_eQTL


rm(list = ls())
###E6_1~E6_6  MAF
MAF_E6 <- read.table(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/MAF_E6.txt", fill=TRUE, header = T, as.is = T)
head(MAF_E6)
MAF_E6=MAF_E6[ , c(1,4)]
colnames(MAF_E6) <-c("rs_id", "maf")
head(MAF_E6)

#gwas
gwas_E6 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/gwaschr2_rs1873164.bed", fill=TRUE, header = T, as.is = T)

#E6_1
E6_1_TTN <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E6_1_TTN.bed", fill=TRUE, header = T, as.is = T)
E6_1_TTN=E6_1_TTN[ , c(4,6,10)]
colnames(E6_1_TTN) <-c("rs_id", "BP","pval_nominal")
fE6_1_TTN <- merge(MAF_E6,E6_1_TTN, by = "rs_id")
head(fE6_1_TTN)
dim(fE6_1_TTN)
write.table(fE6_1_TTN,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE6_1_TTN.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E6_1 <- merge(fE6_1_TTN, gwas_E6, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E6_1)
#coloc 
input_E6_1 <- input_E6_1[complete.cases(input_E6_1),]
coloc_E6_1 <- coloc.abf(
  dataset1=list(snp=input_E6_1$rs_id,pvalues=input_E6_1$P, type="cc", s=0.51, N=nrow(gwas_E6)),
  dataset2=list(snp=input_E6_1$rs_id,pvalues=input_E6_1$pval_nominal, type="quant", N=nrow(fE6_1_TTN)),
  MAF=input_E6_1$maf)

#E6_2
E6_2_TTN <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E6_2_TTN.bed", fill=TRUE, header = T, as.is = T)
E6_2_TTN=E6_2_TTN[ , c(4,6,10)]
colnames(E6_2_TTN) <-c("rs_id", "BP","pval_nominal")
fE6_2_TTN <- merge(MAF_E6,E6_2_TTN, by = "rs_id")
head(fE6_2_TTN)
dim(fE6_2_TTN)
write.table(fE6_2_TTN,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE6_2_TTN.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E6_2 <- merge(fE6_2_TTN, gwas_E6, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E6_2)
#coloc 
input_E6_2 <- input_E6_2[complete.cases(input_E6_2),]
coloc_E6_2 <- coloc.abf(
  dataset1=list(snp=input_E6_2$rs_id,pvalues=input_E6_2$P, type="cc", s=0.51, N=nrow(gwas_E6)),
  dataset2=list(snp=input_E6_2$rs_id,pvalues=input_E6_2$pval_nominal, type="quant", N=nrow(fE6_2_TTN)),
  MAF=input_E6_2$maf)

#E6_3
E6_3_TTN <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E6_3_TTN.bed", fill=TRUE, header = T, as.is = T)
E6_3_TTN=E6_3_TTN[ , c(4,6,10)]
colnames(E6_3_TTN) <-c("rs_id", "BP","pval_nominal")
fE6_3_TTN <- merge(MAF_E6,E6_3_TTN, by = "rs_id")
head(fE6_3_TTN)
dim(fE6_3_TTN)
write.table(fE6_3_TTN,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE6_3_TTN.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E6_3 <- merge(fE6_3_TTN, gwas_E6, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E6_3)
#coloc 
input_E6_3 <- input_E6_3[complete.cases(input_E6_3),]
coloc_E6_3 <- coloc.abf(
  dataset1=list(snp=input_E6_3$rs_id,pvalues=input_E6_3$P, type="cc", s=0.51, N=nrow(gwas_E6)),
  dataset2=list(snp=input_E6_3$rs_id,pvalues=input_E6_3$pval_nominal, type="quant", N=nrow(fE6_3_TTN)),
  MAF=input_E6_3$maf)

#E6_4
E6_4_TTN <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E6_4_TTN.bed", fill=TRUE, header = T, as.is = T)
E6_4_TTN=E6_4_TTN[ , c(4,6,10)]
colnames(E6_4_TTN) <-c("rs_id", "BP","pval_nominal")
fE6_4_TTN <- merge(MAF_E6,E6_4_TTN, by = "rs_id")
head(fE6_4_TTN)
dim(fE6_4_TTN)
write.table(fE6_4_TTN,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE6_4_TTN.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E6_4 <- merge(fE6_4_TTN, gwas_E6, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E6_4)
#coloc 
input_E6_4 <- input_E6_4[complete.cases(input_E6_4),]
coloc_E6_4 <- coloc.abf(
  dataset1=list(snp=input_E6_4$rs_id,pvalues=input_E6_4$P, type="cc", s=0.51, N=nrow(gwas_E6)),
  dataset2=list(snp=input_E6_4$rs_id,pvalues=input_E6_4$pval_nominal, type="quant", N=nrow(fE6_4_TTN)),
  MAF=input_E6_4$maf)

#E6_5
E6_5_TTN <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E6_5_TTN.bed", fill=TRUE, header = T, as.is = T)
E6_5_TTN=E6_5_TTN[ , c(4,6,10)]
colnames(E6_5_TTN) <-c("rs_id", "BP","pval_nominal")
fE6_5_TTN <- merge(MAF_E6,E6_5_TTN, by = "rs_id")
head(fE6_5_TTN)
dim(fE6_5_TTN)
write.table(fE6_5_TTN,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE6_5_TTN.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E6_5 <- merge(fE6_5_TTN, gwas_E6, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E6_5)
#coloc 
input_E6_5 <- input_E6_5[complete.cases(input_E6_5),]
coloc_E6_5 <- coloc.abf(
  dataset1=list(snp=input_E6_5$rs_id,pvalues=input_E6_5$P, type="cc", s=0.51, N=nrow(gwas_E6)),
  dataset2=list(snp=input_E6_5$rs_id,pvalues=input_E6_5$pval_nominal, type="quant", N=nrow(fE6_5_TTN)),
  MAF=input_E6_5$maf)

#E6_6
E6_6_TTN <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E6_6_TTN.bed", fill=TRUE, header = T, as.is = T)
E6_6_TTN=E6_6_TTN[ , c(4,6,10)]
colnames(E6_6_TTN) <-c("rs_id", "BP","pval_nominal")
fE6_6_TTN <- merge(MAF_E6,E6_6_TTN, by = "rs_id")
head(fE6_6_TTN)
dim(fE6_6_TTN)
write.table(fE6_6_TTN,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE6_6_TTN.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E6_6 <- merge(fE6_6_TTN, gwas_E6, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E6_6)
#coloc 
input_E6_6 <- input_E6_6[complete.cases(input_E6_6),]
coloc_E6_6 <- coloc.abf(
  dataset1=list(snp=input_E6_6$rs_id,pvalues=input_E6_6$P, type="cc", s=0.51, N=nrow(gwas_E6)),
  dataset2=list(snp=input_E6_6$rs_id,pvalues=input_E6_6$pval_nominal, type="quant", N=nrow(fE6_6_TTN)),
  MAF=input_E6_6$maf)


rm(list = ls())
###E7E11_1~E7E11_6  MAF
MAF_E7E11 <- read.table(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/MAF_E7E11.txt", fill=TRUE, header = T, as.is = T)
head(MAF_E7E11)
MAF_E7E11=MAF_E7E11[ , c(1,4)]
colnames(MAF_E7E11) <-c("rs_id", "maf")
head(MAF_E7E11)

#gwas
gwas_E7E11 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/Comgwaschr2_rs2220127.bed", fill=TRUE, header = T, as.is = T)

#E7_1
E7_1_TTN <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E7_1_TTN.bed", fill=TRUE, header = T, as.is = T)
E7_1_TTN=E7_1_TTN[ , c(4,6,10)]
colnames(E7_1_TTN) <-c("rs_id", "BP","pval_nominal")
fE7_1_TTN <- merge(MAF_E7E11,E7_1_TTN, by = "rs_id")
head(fE7_1_TTN)
dim(fE7_1_TTN)
write.table(fE7_1_TTN,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE7_1_TTN.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E7_1 <- merge(fE7_1_TTN, gwas_E7E11, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E7_1)
#coloc 
input_E7_1 <- input_E7_1[complete.cases(input_E7_1),]
coloc_E7_1 <- coloc.abf(
  dataset1=list(snp=input_E7_1$rs_id,pvalues=input_E7_1$P, type="cc", s=0.51, N=nrow(gwas_E7E11)),
  dataset2=list(snp=input_E7_1$rs_id,pvalues=input_E7_1$pval_nominal, type="quant", N=nrow(fE7_1_TTN)),
  MAF=input_E7_1$maf)

#E7_2
E7_2_TTN <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E7_2_TTN.bed", fill=TRUE, header = T, as.is = T)
E7_2_TTN=E7_2_TTN[ , c(4,6,10)]
colnames(E7_2_TTN) <-c("rs_id", "BP","pval_nominal")
fE7_2_TTN <- merge(MAF_E7E11,E7_2_TTN, by = "rs_id")
head(fE7_2_TTN)
dim(fE7_2_TTN)
write.table(fE7_2_TTN,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE7_2_TTN.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E7_2 <- merge(fE7_2_TTN, gwas_E7E11, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E7_2)
#coloc 
input_E7_2 <- input_E7_2[complete.cases(input_E7_2),]
coloc_E7_2 <- coloc.abf(
  dataset1=list(snp=input_E7_2$rs_id,pvalues=input_E7_2$P, type="cc", s=0.51, N=nrow(gwas_E7E11)),
  dataset2=list(snp=input_E7_2$rs_id,pvalues=input_E7_2$pval_nominal, type="quant", N=nrow(fE7_2_TTN)),
  MAF=input_E7_2$maf)

#E7_3
E7_3_TTN <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E7_3_TTN.bed", fill=TRUE, header = T, as.is = T)
E7_3_TTN=E7_3_TTN[ , c(4,6,10)]
colnames(E7_3_TTN) <-c("rs_id", "BP","pval_nominal")
fE7_3_TTN <- merge(MAF_E7E11,E7_3_TTN, by = "rs_id")
head(fE7_3_TTN)
dim(fE7_3_TTN)
write.table(fE7_3_TTN,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE7_3_TTN.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E7_3 <- merge(fE7_3_TTN, gwas_E7E11, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E7_3)
#coloc 
input_E7_3 <- input_E7_3[complete.cases(input_E7_3),]
coloc_E7_3 <- coloc.abf(
  dataset1=list(snp=input_E7_3$rs_id,pvalues=input_E7_3$P, type="cc", s=0.51, N=nrow(gwas_E7E11)),
  dataset2=list(snp=input_E7_3$rs_id,pvalues=input_E7_3$pval_nominal, type="quant", N=nrow(fE7_3_TTN)),
  MAF=input_E7_3$maf)

#E7_4
E7_4_TTN <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E7_4_TTN.bed", fill=TRUE, header = T, as.is = T)
E7_4_TTN=E7_4_TTN[ , c(4,6,10)]
colnames(E7_4_TTN) <-c("rs_id", "BP","pval_nominal")
fE7_4_TTN <- merge(MAF_E7E11,E7_4_TTN, by = "rs_id")
head(fE7_4_TTN)
dim(fE7_4_TTN)
write.table(fE7_4_TTN,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE7_4_TTN.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E7_4 <- merge(fE7_4_TTN, gwas_E7E11, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E7_4)
#coloc 
input_E7_4 <- input_E7_4[complete.cases(input_E7_4),]
coloc_E7_4 <- coloc.abf(
  dataset1=list(snp=input_E7_4$rs_id,pvalues=input_E7_4$P, type="cc", s=0.51, N=nrow(gwas_E7E11)),
  dataset2=list(snp=input_E7_4$rs_id,pvalues=input_E7_4$pval_nominal, type="quant", N=nrow(fE7_4_TTN)),
  MAF=input_E7_4$maf)

#E7_5
E7_5_TTN <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E7_5_TTN.bed", fill=TRUE, header = T, as.is = T)
E7_5_TTN=E7_5_TTN[ , c(4,6,10)]
colnames(E7_5_TTN) <-c("rs_id", "BP","pval_nominal")
fE7_5_TTN <- merge(MAF_E7E11,E7_5_TTN, by = "rs_id")
head(fE7_5_TTN)
dim(fE7_5_TTN)
write.table(fE7_5_TTN,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE7_5_TTN.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E7_5 <- merge(fE7_5_TTN, gwas_E7E11, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E7_5)
#coloc 
input_E7_5 <- input_E7_5[complete.cases(input_E7_5),]
coloc_E7_5 <- coloc.abf(
  dataset1=list(snp=input_E7_5$rs_id,pvalues=input_E7_5$P, type="cc", s=0.51, N=nrow(gwas_E7E11)),
  dataset2=list(snp=input_E7_5$rs_id,pvalues=input_E7_5$pval_nominal, type="quant", N=nrow(fE7_5_TTN)),
  MAF=input_E7_5$maf)

#E7_6
E7_6_TTN <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E7_6_TTN.bed", fill=TRUE, header = T, as.is = T)
E7_6_TTN=E7_6_TTN[ , c(4,6,10)]
colnames(E7_6_TTN) <-c("rs_id", "BP","pval_nominal")
fE7_6_TTN <- merge(MAF_E7E11,E7_6_TTN, by = "rs_id")
head(fE7_6_TTN)
dim(fE7_6_TTN)
write.table(fE7_6_TTN,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE7_6_TTN.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E7_6 <- merge(fE7_6_TTN, gwas_E7E11, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E7_6)
#coloc 
input_E7_6 <- input_E7_6[complete.cases(input_E7_6),]
coloc_E7_6 <- coloc.abf(
  dataset1=list(snp=input_E7_6$rs_id,pvalues=input_E7_6$P, type="cc", s=0.51, N=nrow(gwas_E7E11)),
  dataset2=list(snp=input_E7_6$rs_id,pvalues=input_E7_6$pval_nominal, type="quant", N=nrow(fE7_6_TTN)),
  MAF=input_E7_6$maf)

#E11_1
E11_1_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E11_1_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E11_1_FKBP7=E11_1_FKBP7[ , c(4,6,10)]
colnames(E11_1_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE11_1_FKBP7 <- merge(MAF_E7E11,E11_1_FKBP7, by = "rs_id")
head(fE11_1_FKBP7)
dim(fE11_1_FKBP7)
write.table(fE11_1_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE11_1_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E11_1 <- merge(fE11_1_FKBP7, gwas_E7E11, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E11_1)
#coloc 
input_E11_1 <- input_E11_1[complete.cases(input_E11_1),]
coloc_E11_1 <- coloc.abf(
  dataset1=list(snp=input_E11_1$rs_id,pvalues=input_E11_1$P, type="cc", s=0.51, N=nrow(gwas_E7E11)),
  dataset2=list(snp=input_E11_1$rs_id,pvalues=input_E11_1$pval_nominal, type="quant", N=nrow(fE11_1_FKBP7)),
  MAF=input_E11_1$maf)

#E11_2
E11_2_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E11_2_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E11_2_FKBP7=E11_2_FKBP7[ , c(4,6,10)]
colnames(E11_2_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE11_2_FKBP7 <- merge(MAF_E7E11,E11_2_FKBP7, by = "rs_id")
head(fE11_2_FKBP7)
dim(fE11_2_FKBP7)
write.table(fE11_2_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE11_2_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E11_2 <- merge(fE11_2_FKBP7, gwas_E7E11, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E11_2)
#coloc 
input_E11_2 <- input_E11_2[complete.cases(input_E11_2),]
coloc_E11_2 <- coloc.abf(
  dataset1=list(snp=input_E11_2$rs_id,pvalues=input_E11_2$P, type="cc", s=0.51, N=nrow(gwas_E7E11)),
  dataset2=list(snp=input_E11_2$rs_id,pvalues=input_E11_2$pval_nominal, type="quant", N=nrow(fE11_2_FKBP7)),
  MAF=input_E11_2$maf)

#E11_3
E11_3_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E11_3_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E11_3_FKBP7=E11_3_FKBP7[ , c(4,6,10)]
colnames(E11_3_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE11_3_FKBP7 <- merge(MAF_E7E11,E11_3_FKBP7, by = "rs_id")
head(fE11_3_FKBP7)
dim(fE11_3_FKBP7)
write.table(fE11_3_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE11_3_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E11_3 <- merge(fE11_3_FKBP7, gwas_E7E11, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E11_3)
#coloc 
input_E11_3 <- input_E11_3[complete.cases(input_E11_3),]
coloc_E11_3 <- coloc.abf(
  dataset1=list(snp=input_E11_3$rs_id,pvalues=input_E11_3$P, type="cc", s=0.51, N=nrow(gwas_E7E11)),
  dataset2=list(snp=input_E11_3$rs_id,pvalues=input_E11_3$pval_nominal, type="quant", N=nrow(fE11_3_FKBP7)),
  MAF=input_E11_3$maf)

#E11_4
E11_4_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E11_4_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E11_4_FKBP7=E11_4_FKBP7[ , c(4,6,10)]
colnames(E11_4_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE11_4_FKBP7 <- merge(MAF_E7E11,E11_4_FKBP7, by = "rs_id")
head(fE11_4_FKBP7)
dim(fE11_4_FKBP7)
write.table(fE11_4_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE11_4_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E11_4 <- merge(fE11_4_FKBP7, gwas_E7E11, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E11_4)
#coloc 
input_E11_4 <- input_E11_4[complete.cases(input_E11_4),]
coloc_E11_4 <- coloc.abf(
  dataset1=list(snp=input_E11_4$rs_id,pvalues=input_E11_4$P, type="cc", s=0.51, N=nrow(gwas_E7E11)),
  dataset2=list(snp=input_E11_4$rs_id,pvalues=input_E11_4$pval_nominal, type="quant", N=nrow(fE11_4_FKBP7)),
  MAF=input_E11_4$maf)

#E11_5
E11_5_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E11_5_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E11_5_FKBP7=E11_5_FKBP7[ , c(4,6,10)]
colnames(E11_5_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE11_5_FKBP7 <- merge(MAF_E7E11,E11_5_FKBP7, by = "rs_id")
head(fE11_5_FKBP7)
dim(fE11_5_FKBP7)
write.table(fE11_5_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE11_5_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E11_5 <- merge(fE11_5_FKBP7, gwas_E7E11, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E11_5)
#coloc 
input_E11_5 <- input_E11_5[complete.cases(input_E11_5),]
coloc_E11_5 <- coloc.abf(
  dataset1=list(snp=input_E11_5$rs_id,pvalues=input_E11_5$P, type="cc", s=0.51, N=nrow(gwas_E7E11)),
  dataset2=list(snp=input_E11_5$rs_id,pvalues=input_E11_5$pval_nominal, type="quant", N=nrow(fE11_5_FKBP7)),
  MAF=input_E11_5$maf)

#E11_6
E11_6_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E11_6_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E11_6_FKBP7=E11_6_FKBP7[ , c(4,6,10)]
colnames(E11_6_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE11_6_FKBP7 <- merge(MAF_E7E11,E11_6_FKBP7, by = "rs_id")
head(fE11_6_FKBP7)
dim(fE11_6_FKBP7)
write.table(fE11_6_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE11_6_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E11_6 <- merge(fE11_6_FKBP7, gwas_E7E11, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E11_6)
#coloc 
input_E11_6 <- input_E11_6[complete.cases(input_E11_6),]
coloc_E11_6 <- coloc.abf(
  dataset1=list(snp=input_E11_6$rs_id,pvalues=input_E11_6$P, type="cc", s=0.51, N=nrow(gwas_E7E11)),
  dataset2=list(snp=input_E11_6$rs_id,pvalues=input_E11_6$pval_nominal, type="quant", N=nrow(fE11_6_FKBP7)),
  MAF=input_E11_6$maf)


rm(list = ls())
###E8E12_1~E8E12_6  MAF
MAF_E8E12 <- read.table(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/MAF_E8E12.txt", fill=TRUE, header = T, as.is = T)
head(MAF_E8E12)
MAF_E8E12=MAF_E8E12[ , c(1,4)]
colnames(MAF_E8E12) <-c("rs_id", "maf")
head(MAF_E8E12)

#gwas
gwas_E8E12 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/gwaschr2_rs2562845.bed", fill=TRUE, header = T, as.is = T)

#E8_1
E8_1_TTN <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E8_1_TTN.bed", fill=TRUE, header = T, as.is = T)
E8_1_TTN=E8_1_TTN[ , c(4,6,10)]
colnames(E8_1_TTN) <-c("rs_id", "BP","pval_nominal")
fE8_1_TTN <- merge(MAF_E8E12,E8_1_TTN, by = "rs_id")
head(fE8_1_TTN)
dim(fE8_1_TTN)
write.table(fE8_1_TTN,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE8_1_TTN.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E8_1 <- merge(fE8_1_TTN, gwas_E8E12, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E8_1)
#coloc 
input_E8_1 <- input_E8_1[complete.cases(input_E8_1),]
coloc_E8_1 <- coloc.abf(
  dataset1=list(snp=input_E8_1$rs_id,pvalues=input_E8_1$P, type="cc", s=0.51, N=nrow(gwas_E8E12)),
  dataset2=list(snp=input_E8_1$rs_id,pvalues=input_E8_1$pval_nominal, type="quant", N=nrow(fE8_1_TTN)),
  MAF=input_E8_1$maf)

#E8_2
E8_2_TTN <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E8_2_TTN.bed", fill=TRUE, header = T, as.is = T)
E8_2_TTN=E8_2_TTN[ , c(4,6,10)]
colnames(E8_2_TTN) <-c("rs_id", "BP","pval_nominal")
fE8_2_TTN <- merge(MAF_E8E12,E8_2_TTN, by = "rs_id")
head(fE8_2_TTN)
dim(fE8_2_TTN)
write.table(fE8_2_TTN,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE8_2_TTN.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E8_2 <- merge(fE8_2_TTN, gwas_E8E12, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E8_2)
#coloc 
input_E8_2 <- input_E8_2[complete.cases(input_E8_2),]
coloc_E8_2 <- coloc.abf(
  dataset1=list(snp=input_E8_2$rs_id,pvalues=input_E8_2$P, type="cc", s=0.51, N=nrow(gwas_E8E12)),
  dataset2=list(snp=input_E8_2$rs_id,pvalues=input_E8_2$pval_nominal, type="quant", N=nrow(fE8_2_TTN)),
  MAF=input_E8_2$maf)

#E8_3
E8_3_TTN <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E8_3_TTN.bed", fill=TRUE, header = T, as.is = T)
E8_3_TTN=E8_3_TTN[ , c(4,6,10)]
colnames(E8_3_TTN) <-c("rs_id", "BP","pval_nominal")
fE8_3_TTN <- merge(MAF_E8E12,E8_3_TTN, by = "rs_id")
head(fE8_3_TTN)
dim(fE8_3_TTN)
write.table(fE8_3_TTN,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE8_3_TTN.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E8_3 <- merge(fE8_3_TTN, gwas_E8E12, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E8_3)
#coloc 
input_E8_3 <- input_E8_3[complete.cases(input_E8_3),]
coloc_E8_3 <- coloc.abf(
  dataset1=list(snp=input_E8_3$rs_id,pvalues=input_E8_3$P, type="cc", s=0.51, N=nrow(gwas_E8E12)),
  dataset2=list(snp=input_E8_3$rs_id,pvalues=input_E8_3$pval_nominal, type="quant", N=nrow(fE8_3_TTN)),
  MAF=input_E8_3$maf)

#E8_4
E8_4_TTN <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E8_4_TTN.bed", fill=TRUE, header = T, as.is = T)
E8_4_TTN=E8_4_TTN[ , c(4,6,10)]
colnames(E8_4_TTN) <-c("rs_id", "BP","pval_nominal")
fE8_4_TTN <- merge(MAF_E8E12,E8_4_TTN, by = "rs_id")
head(fE8_4_TTN)
dim(fE8_4_TTN)
write.table(fE8_4_TTN,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE8_4_TTN.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E8_4 <- merge(fE8_4_TTN, gwas_E8E12, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E8_4)
#coloc 
input_E8_4 <- input_E8_4[complete.cases(input_E8_4),]
coloc_E8_4 <- coloc.abf(
  dataset1=list(snp=input_E8_4$rs_id,pvalues=input_E8_4$P, type="cc", s=0.51, N=nrow(gwas_E8E12)),
  dataset2=list(snp=input_E8_4$rs_id,pvalues=input_E8_4$pval_nominal, type="quant", N=nrow(fE8_4_TTN)),
  MAF=input_E8_4$maf)

#E8_5
E8_5_TTN <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E8_5_TTN.bed", fill=TRUE, header = T, as.is = T)
E8_5_TTN=E8_5_TTN[ , c(4,6,10)]
colnames(E8_5_TTN) <-c("rs_id", "BP","pval_nominal")
fE8_5_TTN <- merge(MAF_E8E12,E8_5_TTN, by = "rs_id")
head(fE8_5_TTN)
dim(fE8_5_TTN)
write.table(fE8_5_TTN,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE8_5_TTN.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E8_5 <- merge(fE8_5_TTN, gwas_E8E12, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E8_5)
#coloc 
input_E8_5 <- input_E8_5[complete.cases(input_E8_5),]
coloc_E8_5 <- coloc.abf(
  dataset1=list(snp=input_E8_5$rs_id,pvalues=input_E8_5$P, type="cc", s=0.51, N=nrow(gwas_E8E12)),
  dataset2=list(snp=input_E8_5$rs_id,pvalues=input_E8_5$pval_nominal, type="quant", N=nrow(fE8_5_TTN)),
  MAF=input_E8_5$maf)

#E8_6
E8_6_TTN <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E8_6_TTN.bed", fill=TRUE, header = T, as.is = T)
E8_6_TTN=E8_6_TTN[ , c(4,6,10)]
colnames(E8_6_TTN) <-c("rs_id", "BP","pval_nominal")
fE8_6_TTN <- merge(MAF_E8E12,E8_6_TTN, by = "rs_id")
head(fE8_6_TTN)
dim(fE8_6_TTN)
write.table(fE8_6_TTN,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE8_6_TTN.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E8_6 <- merge(fE8_6_TTN, gwas_E8E12, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E8_6)
#coloc 
input_E8_6 <- input_E8_6[complete.cases(input_E8_6),]
coloc_E8_6 <- coloc.abf(
  dataset1=list(snp=input_E8_6$rs_id,pvalues=input_E8_6$P, type="cc", s=0.51, N=nrow(gwas_E8E12)),
  dataset2=list(snp=input_E8_6$rs_id,pvalues=input_E8_6$pval_nominal, type="quant", N=nrow(fE8_6_TTN)),
  MAF=input_E8_6$maf)

#E12_1
E12_1_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E12_1_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E12_1_FKBP7=E12_1_FKBP7[ , c(4,6,10)]
colnames(E12_1_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE12_1_FKBP7 <- merge(MAF_E8E12,E12_1_FKBP7, by = "rs_id")
head(fE12_1_FKBP7)
dim(fE12_1_FKBP7)
write.table(fE12_1_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE12_1_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E12_1 <- merge(fE12_1_FKBP7, gwas_E8E12, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E12_1)
#coloc 
input_E12_1 <- input_E12_1[complete.cases(input_E12_1),]
coloc_E12_1 <- coloc.abf(
  dataset1=list(snp=input_E12_1$rs_id,pvalues=input_E12_1$P, type="cc", s=0.51, N=nrow(gwas_E8E12)),
  dataset2=list(snp=input_E12_1$rs_id,pvalues=input_E12_1$pval_nominal, type="quant", N=nrow(fE12_1_FKBP7)),
  MAF=input_E12_1$maf)

#E12_2
E12_2_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E12_2_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E12_2_FKBP7=E12_2_FKBP7[ , c(4,6,10)]
colnames(E12_2_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE12_2_FKBP7 <- merge(MAF_E8E12,E12_2_FKBP7, by = "rs_id")
head(fE12_2_FKBP7)
dim(fE12_2_FKBP7)
write.table(fE12_2_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE12_2_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E12_2 <- merge(fE12_2_FKBP7, gwas_E8E12, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E12_2)
#coloc 
input_E12_2 <- input_E12_2[complete.cases(input_E12_2),]
coloc_E12_2 <- coloc.abf(
  dataset1=list(snp=input_E12_2$rs_id,pvalues=input_E12_2$P, type="cc", s=0.51, N=nrow(gwas_E8E12)),
  dataset2=list(snp=input_E12_2$rs_id,pvalues=input_E12_2$pval_nominal, type="quant", N=nrow(fE12_2_FKBP7)),
  MAF=input_E12_2$maf)

#E12_3
E12_3_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E12_3_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E12_3_FKBP7=E12_3_FKBP7[ , c(4,6,10)]
colnames(E12_3_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE12_3_FKBP7 <- merge(MAF_E8E12,E12_3_FKBP7, by = "rs_id")
head(fE12_3_FKBP7)
dim(fE12_3_FKBP7)
write.table(fE12_3_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE12_3_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E12_3 <- merge(fE12_3_FKBP7, gwas_E8E12, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E12_3)
#coloc 
input_E12_3 <- input_E12_3[complete.cases(input_E12_3),]
coloc_E12_3 <- coloc.abf(
  dataset1=list(snp=input_E12_3$rs_id,pvalues=input_E12_3$P, type="cc", s=0.51, N=nrow(gwas_E8E12)),
  dataset2=list(snp=input_E12_3$rs_id,pvalues=input_E12_3$pval_nominal, type="quant", N=nrow(fE12_3_FKBP7)),
  MAF=input_E12_3$maf)

#E12_4
E12_4_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E12_4_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E12_4_FKBP7=E12_4_FKBP7[ , c(4,6,10)]
colnames(E12_4_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE12_4_FKBP7 <- merge(MAF_E8E12,E12_4_FKBP7, by = "rs_id")
head(fE12_4_FKBP7)
dim(fE12_4_FKBP7)
write.table(fE12_4_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE12_4_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E12_4 <- merge(fE12_4_FKBP7, gwas_E8E12, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E12_4)
#coloc 
input_E12_4 <- input_E12_4[complete.cases(input_E12_4),]
coloc_E12_4 <- coloc.abf(
  dataset1=list(snp=input_E12_4$rs_id,pvalues=input_E12_4$P, type="cc", s=0.51, N=nrow(gwas_E8E12)),
  dataset2=list(snp=input_E12_4$rs_id,pvalues=input_E12_4$pval_nominal, type="quant", N=nrow(fE12_4_FKBP7)),
  MAF=input_E12_4$maf)

#E12_5
E12_5_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E12_5_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E12_5_FKBP7=E12_5_FKBP7[ , c(4,6,10)]
colnames(E12_5_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE12_5_FKBP7 <- merge(MAF_E8E12,E12_5_FKBP7, by = "rs_id")
head(fE12_5_FKBP7)
dim(fE12_5_FKBP7)
write.table(fE12_5_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE12_5_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E12_5 <- merge(fE12_5_FKBP7, gwas_E8E12, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E12_5)
#coloc 
input_E12_5 <- input_E12_5[complete.cases(input_E12_5),]
coloc_E12_5 <- coloc.abf(
  dataset1=list(snp=input_E12_5$rs_id,pvalues=input_E12_5$P, type="cc", s=0.51, N=nrow(gwas_E8E12)),
  dataset2=list(snp=input_E12_5$rs_id,pvalues=input_E12_5$pval_nominal, type="quant", N=nrow(fE12_5_FKBP7)),
  MAF=input_E12_5$maf)

#E12_6
E12_6_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E12_6_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E12_6_FKBP7=E12_6_FKBP7[ , c(4,6,10)]
colnames(E12_6_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE12_6_FKBP7 <- merge(MAF_E8E12,E12_6_FKBP7, by = "rs_id")
head(fE12_6_FKBP7)
dim(fE12_6_FKBP7)
write.table(fE12_6_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE12_6_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E12_6 <- merge(fE12_6_FKBP7, gwas_E8E12, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E12_6)
#coloc 
input_E12_6 <- input_E12_6[complete.cases(input_E12_6),]
coloc_E12_6 <- coloc.abf(
  dataset1=list(snp=input_E12_6$rs_id,pvalues=input_E12_6$P, type="cc", s=0.51, N=nrow(gwas_E8E12)),
  dataset2=list(snp=input_E12_6$rs_id,pvalues=input_E12_6$pval_nominal, type="quant", N=nrow(fE12_6_FKBP7)),
  MAF=input_E12_6$maf)


rm(list = ls())
###E9_1~E9_6  MAF
MAF_E9 <- read.table(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/MAF_E9.txt", fill=TRUE, header = T, as.is = T)
head(MAF_E9)
MAF_E9=MAF_E9[ , c(1,4)]
colnames(MAF_E9) <-c("rs_id", "maf")
head(MAF_E9)

#gwas
gwas_E9 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/gwaschr2_rs10497529.bed", fill=TRUE, header = T, as.is = T)

#E9_1
E9_1_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E9_1_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E9_1_FKBP7=E9_1_FKBP7[ , c(4,6,10)]
colnames(E9_1_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE9_1_FKBP7 <- merge(MAF_E9,E9_1_FKBP7, by = "rs_id")
head(fE9_1_FKBP7)
dim(fE9_1_FKBP7)
write.table(fE9_1_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE9_1_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E9_1 <- merge(fE9_1_FKBP7, gwas_E9, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E9_1)
#coloc 
input_E9_1 <- input_E9_1[complete.cases(input_E9_1),]
coloc_E9_1 <- coloc.abf(
  dataset1=list(snp=input_E9_1$rs_id,pvalues=input_E9_1$P, type="cc", s=0.51, N=nrow(gwas_E9)),
  dataset2=list(snp=input_E9_1$rs_id,pvalues=input_E9_1$pval_nominal, type="quant", N=nrow(fE9_1_FKBP7)),
  MAF=input_E9_1$maf)

#E9_2
E9_2_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E9_2_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E9_2_FKBP7=E9_2_FKBP7[ , c(4,6,10)]
colnames(E9_2_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE9_2_FKBP7 <- merge(MAF_E9,E9_2_FKBP7, by = "rs_id")
head(fE9_2_FKBP7)
dim(fE9_2_FKBP7)
write.table(fE9_2_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE9_2_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E9_2 <- merge(fE9_2_FKBP7, gwas_E9, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E9_2)
#coloc 
input_E9_2 <- input_E9_2[complete.cases(input_E9_2),]
coloc_E9_2 <- coloc.abf(
  dataset1=list(snp=input_E9_2$rs_id,pvalues=input_E9_2$P, type="cc", s=0.51, N=nrow(gwas_E9)),
  dataset2=list(snp=input_E9_2$rs_id,pvalues=input_E9_2$pval_nominal, type="quant", N=nrow(fE9_2_FKBP7)),
  MAF=input_E9_2$maf)

#E9_3
E9_3_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E9_3_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E9_3_FKBP7=E9_3_FKBP7[ , c(4,6,10)]
colnames(E9_3_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE9_3_FKBP7 <- merge(MAF_E9,E9_3_FKBP7, by = "rs_id")
head(fE9_3_FKBP7)
dim(fE9_3_FKBP7)
write.table(fE9_3_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE9_3_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E9_3 <- merge(fE9_3_FKBP7, gwas_E9, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E9_3)
#coloc 
input_E9_3 <- input_E9_3[complete.cases(input_E9_3),]
coloc_E9_3 <- coloc.abf(
  dataset1=list(snp=input_E9_3$rs_id,pvalues=input_E9_3$P, type="cc", s=0.51, N=nrow(gwas_E9)),
  dataset2=list(snp=input_E9_3$rs_id,pvalues=input_E9_3$pval_nominal, type="quant", N=nrow(fE9_3_FKBP7)),
  MAF=input_E9_3$maf)

#E9_4
E9_4_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E9_4_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E9_4_FKBP7=E9_4_FKBP7[ , c(4,6,10)]
colnames(E9_4_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE9_4_FKBP7 <- merge(MAF_E9,E9_4_FKBP7, by = "rs_id")
head(fE9_4_FKBP7)
dim(fE9_4_FKBP7)
write.table(fE9_4_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE9_4_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E9_4 <- merge(fE9_4_FKBP7, gwas_E9, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E9_4)
#coloc 
input_E9_4 <- input_E9_4[complete.cases(input_E9_4),]
coloc_E9_4 <- coloc.abf(
  dataset1=list(snp=input_E9_4$rs_id,pvalues=input_E9_4$P, type="cc", s=0.51, N=nrow(gwas_E9)),
  dataset2=list(snp=input_E9_4$rs_id,pvalues=input_E9_4$pval_nominal, type="quant", N=nrow(fE9_4_FKBP7)),
  MAF=input_E9_4$maf)

#E9_5
E9_5_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E9_5_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E9_5_FKBP7=E9_5_FKBP7[ , c(4,6,10)]
colnames(E9_5_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE9_5_FKBP7 <- merge(MAF_E9,E9_5_FKBP7, by = "rs_id")
head(fE9_5_FKBP7)
dim(fE9_5_FKBP7)
write.table(fE9_5_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE9_5_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E9_5 <- merge(fE9_5_FKBP7, gwas_E9, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E9_5)
#coloc 
input_E9_5 <- input_E9_5[complete.cases(input_E9_5),]
coloc_E9_5 <- coloc.abf(
  dataset1=list(snp=input_E9_5$rs_id,pvalues=input_E9_5$P, type="cc", s=0.51, N=nrow(gwas_E9)),
  dataset2=list(snp=input_E9_5$rs_id,pvalues=input_E9_5$pval_nominal, type="quant", N=nrow(fE9_5_FKBP7)),
  MAF=input_E9_5$maf)

#E9_6
E9_6_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E9_6_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E9_6_FKBP7=E9_6_FKBP7[ , c(4,6,10)]
colnames(E9_6_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE9_6_FKBP7 <- merge(MAF_E9,E9_6_FKBP7, by = "rs_id")
head(fE9_6_FKBP7)
dim(fE9_6_FKBP7)
write.table(fE9_6_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE9_6_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E9_6 <- merge(fE9_6_FKBP7, gwas_E9, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E9_6)
#coloc 
input_E9_6 <- input_E9_6[complete.cases(input_E9_6),]
coloc_E9_6 <- coloc.abf(
  dataset1=list(snp=input_E9_6$rs_id,pvalues=input_E9_6$P, type="cc", s=0.51, N=nrow(gwas_E9)),
  dataset2=list(snp=input_E9_6$rs_id,pvalues=input_E9_6$pval_nominal, type="quant", N=nrow(fE9_6_FKBP7)),
  MAF=input_E9_6$maf)


rm(list = ls())
###E10_1~E10_6  MAF
MAF_E10 <- read.table(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/MAF_E10.txt", fill=TRUE, header = T, as.is = T)
head(MAF_E10)
MAF_E10=MAF_E10[ , c(1,4)]
colnames(MAF_E10) <-c("rs_id", "maf")
head(MAF_E10)

#gwas
gwas_E10 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/Comgwaschr2_rs142556838.bed", fill=TRUE, header = T, as.is = T)

#E10_1
E10_1_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E10_1_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E10_1_FKBP7=E10_1_FKBP7[ , c(4,6,10)]
colnames(E10_1_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE10_1_FKBP7 <- merge(MAF_E10,E10_1_FKBP7, by = "rs_id")
head(fE10_1_FKBP7)
dim(fE10_1_FKBP7)
write.table(fE10_1_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE10_1_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E10_1 <- merge(fE10_1_FKBP7, gwas_E10, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E10_1)
#coloc 
input_E10_1 <- input_E10_1[complete.cases(input_E10_1),]
coloc_E10_1 <- coloc.abf(
  dataset1=list(snp=input_E10_1$rs_id,pvalues=input_E10_1$P, type="cc", s=0.51, N=nrow(gwas_E10)),
  dataset2=list(snp=input_E10_1$rs_id,pvalues=input_E10_1$pval_nominal, type="quant", N=nrow(fE10_1_FKBP7)),
  MAF=input_E10_1$maf)

#E10_2
E10_2_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E10_2_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E10_2_FKBP7=E10_2_FKBP7[ , c(4,6,10)]
colnames(E10_2_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE10_2_FKBP7 <- merge(MAF_E10,E10_2_FKBP7, by = "rs_id")
head(fE10_2_FKBP7)
dim(fE10_2_FKBP7)
write.table(fE10_2_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE10_2_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E10_2 <- merge(fE10_2_FKBP7, gwas_E10, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E10_2)
#coloc 
input_E10_2 <- input_E10_2[complete.cases(input_E10_2),]
coloc_E10_2 <- coloc.abf(
  dataset1=list(snp=input_E10_2$rs_id,pvalues=input_E10_2$P, type="cc", s=0.51, N=nrow(gwas_E10)),
  dataset2=list(snp=input_E10_2$rs_id,pvalues=input_E10_2$pval_nominal, type="quant", N=nrow(fE10_2_FKBP7)),
  MAF=input_E10_2$maf)

#E10_3
E10_3_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E10_3_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E10_3_FKBP7=E10_3_FKBP7[ , c(4,6,10)]
colnames(E10_3_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE10_3_FKBP7 <- merge(MAF_E10,E10_3_FKBP7, by = "rs_id")
head(fE10_3_FKBP7)
dim(fE10_3_FKBP7)
write.table(fE10_3_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE10_3_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E10_3 <- merge(fE10_3_FKBP7, gwas_E10, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E10_3)
#coloc 
input_E10_3 <- input_E10_3[complete.cases(input_E10_3),]
coloc_E10_3 <- coloc.abf(
  dataset1=list(snp=input_E10_3$rs_id,pvalues=input_E10_3$P, type="cc", s=0.51, N=nrow(gwas_E10)),
  dataset2=list(snp=input_E10_3$rs_id,pvalues=input_E10_3$pval_nominal, type="quant", N=nrow(fE10_3_FKBP7)),
  MAF=input_E10_3$maf)

#E10_4
E10_4_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E10_4_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E10_4_FKBP7=E10_4_FKBP7[ , c(4,6,10)]
colnames(E10_4_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE10_4_FKBP7 <- merge(MAF_E10,E10_4_FKBP7, by = "rs_id")
head(fE10_4_FKBP7)
dim(fE10_4_FKBP7)
write.table(fE10_4_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE10_4_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E10_4 <- merge(fE10_4_FKBP7, gwas_E10, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E10_4)
#coloc 
input_E10_4 <- input_E10_4[complete.cases(input_E10_4),]
coloc_E10_4 <- coloc.abf(
  dataset1=list(snp=input_E10_4$rs_id,pvalues=input_E10_4$P, type="cc", s=0.51, N=nrow(gwas_E10)),
  dataset2=list(snp=input_E10_4$rs_id,pvalues=input_E10_4$pval_nominal, type="quant", N=nrow(fE10_4_FKBP7)),
  MAF=input_E10_4$maf)

#E10_5
E10_5_FKBP7 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E10_5_FKBP7.bed", fill=TRUE, header = T, as.is = T)
E10_5_FKBP7=E10_5_FKBP7[ , c(4,6,10)]
colnames(E10_5_FKBP7) <-c("rs_id", "BP","pval_nominal")
fE10_5_FKBP7 <- merge(MAF_E10,E10_5_FKBP7, by = "rs_id")
head(fE10_5_FKBP7)
dim(fE10_5_FKBP7)
write.table(fE10_5_FKBP7,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE10_5_FKBP7.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E10_5 <- merge(fE10_5_FKBP7, gwas_E10, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E10_5)
#coloc 
input_E10_5 <- input_E10_5[complete.cases(input_E10_5),]
coloc_E10_5 <- coloc.abf(
  dataset1=list(snp=input_E10_5$rs_id,pvalues=input_E10_5$P, type="cc", s=0.51, N=nrow(gwas_E10)),
  dataset2=list(snp=input_E10_5$rs_id,pvalues=input_E10_5$pval_nominal, type="quant", N=nrow(fE10_5_FKBP7)),
  MAF=input_E10_5$maf)

#E10_6   无FKBP7_eQTL_E10_6数据


rm(list = ls())
###E13_1~E13_6  MAF
MAF_E13 <- read.table(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/MAF_E13.txt", fill=TRUE, header = T, as.is = T)
head(MAF_E13)
MAF_E13=MAF_E13[ , c(1,4)]
colnames(MAF_E13) <-c("rs_id", "maf")
head(MAF_E13)

#gwas
gwas_E13 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/norgwaschr16_rsrs67329386.bed", fill=TRUE, header = T, as.is = T)

#E13_1
E13_1_HP <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E13_1_HP.bed", fill=TRUE, header = T, as.is = T)
E13_1_HP=E13_1_HP[ , c(4,6,10)]
colnames(E13_1_HP) <-c("rs_id", "BP","pval_nominal")
fE13_1_HP <- merge(MAF_E13,E13_1_HP, by = "rs_id")
head(fE13_1_HP)
dim(fE13_1_HP)
write.table(fE13_1_HP,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE13_1_HP.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E13_1 <- merge(fE13_1_HP, gwas_E13, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E13_1)
#coloc 
input_E13_1 <- input_E13_1[complete.cases(input_E13_1),]
coloc_E13_1 <- coloc.abf(
  dataset1=list(snp=input_E13_1$rs_id,pvalues=input_E13_1$P, type="cc", s=0.51, N=nrow(gwas_E13)),
  dataset2=list(snp=input_E13_1$rs_id,pvalues=input_E13_1$pval_nominal, type="quant", N=nrow(fE13_1_HP)),
  MAF=input_E13_1$maf)

#E13_2
E13_2_HP <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E13_2_HP.bed", fill=TRUE, header = T, as.is = T)
E13_2_HP=E13_2_HP[ , c(4,6,10)]
colnames(E13_2_HP) <-c("rs_id", "BP","pval_nominal")
fE13_2_HP <- merge(MAF_E13,E13_2_HP, by = "rs_id")
head(fE13_2_HP)
dim(fE13_2_HP)
write.table(fE13_2_HP,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE13_2_HP.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E13_2 <- merge(fE13_2_HP, gwas_E13, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E13_2)
#coloc 
input_E13_2 <- input_E13_2[complete.cases(input_E13_2),]
coloc_E13_2 <- coloc.abf(
  dataset1=list(snp=input_E13_2$rs_id,pvalues=input_E13_2$P, type="cc", s=0.51, N=nrow(gwas_E13)),
  dataset2=list(snp=input_E13_2$rs_id,pvalues=input_E13_2$pval_nominal, type="quant", N=nrow(fE13_2_HP)),
  MAF=input_E13_2$maf)

#E13_3
E13_3_HP <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E13_3_HP.bed", fill=TRUE, header = T, as.is = T)
E13_3_HP=E13_3_HP[ , c(4,6,10)]
colnames(E13_3_HP) <-c("rs_id", "BP","pval_nominal")
fE13_3_HP <- merge(MAF_E13,E13_3_HP, by = "rs_id")
head(fE13_3_HP)
dim(fE13_3_HP)
write.table(fE13_3_HP,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE13_3_HP.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E13_3 <- merge(fE13_3_HP, gwas_E13, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E13_3)
#coloc 
input_E13_3 <- input_E13_3[complete.cases(input_E13_3),]
coloc_E13_3 <- coloc.abf(
  dataset1=list(snp=input_E13_3$rs_id,pvalues=input_E13_3$P, type="cc", s=0.51, N=nrow(gwas_E13)),
  dataset2=list(snp=input_E13_3$rs_id,pvalues=input_E13_3$pval_nominal, type="quant", N=nrow(fE13_3_HP)),
  MAF=input_E13_3$maf)

#E13_4
E13_4_HP <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E13_4_HP.bed", fill=TRUE, header = T, as.is = T)
E13_4_HP=E13_4_HP[ , c(4,6,10)]
colnames(E13_4_HP) <-c("rs_id", "BP","pval_nominal")
fE13_4_HP <- merge(MAF_E13,E13_4_HP, by = "rs_id")
head(fE13_4_HP)
dim(fE13_4_HP)
write.table(fE13_4_HP,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE13_4_HP.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E13_4 <- merge(fE13_4_HP, gwas_E13, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E13_4)
#coloc 
input_E13_4 <- input_E13_4[complete.cases(input_E13_4),]
coloc_E13_4 <- coloc.abf(
  dataset1=list(snp=input_E13_4$rs_id,pvalues=input_E13_4$P, type="cc", s=0.51, N=nrow(gwas_E13)),
  dataset2=list(snp=input_E13_4$rs_id,pvalues=input_E13_4$pval_nominal, type="quant", N=nrow(fE13_4_HP)),
  MAF=input_E13_4$maf)

#E13_5
E13_5_HP <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E13_5_HP.bed", fill=TRUE, header = T, as.is = T)
E13_5_HP=E13_5_HP[ , c(4,6,10)]
colnames(E13_5_HP) <-c("rs_id", "BP","pval_nominal")
fE13_5_HP <- merge(MAF_E13,E13_5_HP, by = "rs_id")
head(fE13_5_HP)
dim(fE13_5_HP)
write.table(fE13_5_HP,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE13_5_HP.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E13_5 <- merge(fE13_5_HP, gwas_E13, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E13_5)
#coloc 
input_E13_5 <- input_E13_5[complete.cases(input_E13_5),]
coloc_E13_5 <- coloc.abf(
  dataset1=list(snp=input_E13_5$rs_id,pvalues=input_E13_5$P, type="cc", s=0.51, N=nrow(gwas_E13)),
  dataset2=list(snp=input_E13_5$rs_id,pvalues=input_E13_5$pval_nominal, type="quant", N=nrow(fE13_5_HP)),
  MAF=input_E13_5$maf)

#E13_6
E13_6_HP <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E13_6_HP.bed", fill=TRUE, header = T, as.is = T)
E13_6_HP=E13_6_HP[ , c(4,6,10)]
colnames(E13_6_HP) <-c("rs_id", "BP","pval_nominal")
fE13_6_HP <- merge(MAF_E13,E13_6_HP, by = "rs_id")
head(fE13_6_HP)
dim(fE13_6_HP)
write.table(fE13_6_HP,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE13_6_HP.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E13_6 <- merge(fE13_6_HP, gwas_E13, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E13_6)
#coloc 
input_E13_6 <- input_E13_6[complete.cases(input_E13_6),]
coloc_E13_6 <- coloc.abf(
  dataset1=list(snp=input_E13_6$rs_id,pvalues=input_E13_6$P, type="cc", s=0.51, N=nrow(gwas_E13)),
  dataset2=list(snp=input_E13_6$rs_id,pvalues=input_E13_6$pval_nominal, type="quant", N=nrow(fE13_6_HP)),
  MAF=input_E13_6$maf)

###E14_1~E14_4  MAF
MAF_E14 <- read.table(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/MAF_E14.txt", fill=TRUE, header = T, as.is = T)
head(MAF_E14)
MAF_E14=MAF_E14[ , c(1,4)]
colnames(MAF_E14) <-c("rs_id", "maf")
head(MAF_E14)

#gwas
gwas_E14 <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/norgwaschr16_rs12325072.bed", fill=TRUE, header = T, as.is = T)

#E14_1
E14_1_HP <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E14_1_HP.bed", fill=TRUE, header = T, as.is = T)
E14_1_HP=E14_1_HP[ , c(4,6,10)]
colnames(E14_1_HP) <-c("rs_id", "BP","pval_nominal")
fE14_1_HP <- merge(MAF_E14,E14_1_HP, by = "rs_id")
head(fE14_1_HP)
dim(fE14_1_HP)
write.table(fE14_1_HP,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE14_1_HP.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E14_1 <- merge(fE14_1_HP, gwas_E14, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E14_1)
#coloc 
input_E14_1 <- input_E14_1[complete.cases(input_E14_1),]
coloc_E14_1 <- coloc.abf(
  dataset1=list(snp=input_E14_1$rs_id,pvalues=input_E14_1$P, type="cc", s=0.51, N=nrow(gwas_E14)),
  dataset2=list(snp=input_E14_1$rs_id,pvalues=input_E14_1$pval_nominal, type="quant", N=nrow(fE14_1_HP)),
  MAF=input_E14_1$maf)

#E14_2
E14_2_HP <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E14_2_HP.bed", fill=TRUE, header = T, as.is = T)
E14_2_HP=E14_2_HP[ , c(4,6,10)]
colnames(E14_2_HP) <-c("rs_id", "BP","pval_nominal")
fE14_2_HP <- merge(MAF_E14,E14_2_HP, by = "rs_id")
head(fE14_2_HP)
dim(fE14_2_HP)
write.table(fE14_2_HP,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE14_2_HP.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E14_2 <- merge(fE14_2_HP, gwas_E14, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E14_2)
#coloc 
input_E14_2 <- input_E14_2[complete.cases(input_E14_2),]
coloc_E14_2 <- coloc.abf(
  dataset1=list(snp=input_E14_2$rs_id,pvalues=input_E14_2$P, type="cc", s=0.51, N=nrow(gwas_E14)),
  dataset2=list(snp=input_E14_2$rs_id,pvalues=input_E14_2$pval_nominal, type="quant", N=nrow(fE14_2_HP)),
  MAF=input_E14_2$maf)

#E14_3
E14_3_HP <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E14_3_HP.bed", fill=TRUE, header = T, as.is = T)
E14_3_HP=E14_3_HP[ , c(4,6,10)]
colnames(E14_3_HP) <-c("rs_id", "BP","pval_nominal")
fE14_3_HP <- merge(MAF_E14,E14_3_HP, by = "rs_id")
head(fE14_3_HP)
dim(fE14_3_HP)
write.table(fE14_3_HP,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE14_3_HP.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E14_3 <- merge(fE14_3_HP, gwas_E14, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E14_3)
#coloc 
input_E14_3 <- input_E14_3[complete.cases(input_E14_3),]
coloc_E14_3 <- coloc.abf(
  dataset1=list(snp=input_E14_3$rs_id,pvalues=input_E14_3$P, type="cc", s=0.51, N=nrow(gwas_E14)),
  dataset2=list(snp=input_E14_3$rs_id,pvalues=input_E14_3$pval_nominal, type="quant", N=nrow(fE14_3_HP)),
  MAF=input_E14_3$maf)

#E14_4
E14_4_HP <- read.delim(file = "/data/wuwq/Figure_wwq/AGAIN/coloc/E14_4_HP.bed", fill=TRUE, header = T, as.is = T)
E14_4_HP=E14_4_HP[ , c(4,6,10)]
colnames(E14_4_HP) <-c("rs_id", "BP","pval_nominal")
fE14_4_HP <- merge(MAF_E14,E14_4_HP, by = "rs_id")
head(fE14_4_HP)
dim(fE14_4_HP)
write.table(fE14_4_HP,file = "/data/wuwq/Figure_wwq/AGAIN/coloc/fE14_4_HP.txt",sep = "\t",row.names=FALSE)
#merge gwas eqtl with shared rs_id
input_E14_4 <- merge(fE14_4_HP, gwas_E14, by="rs_id", all=FALSE, suffixes = c("_eqtl","_gwas"))
head(input_E14_4)
#coloc 
input_E14_4 <- input_E14_4[complete.cases(input_E14_4),]
coloc_E14_4 <- coloc.abf(
  dataset1=list(snp=input_E14_4$rs_id,pvalues=input_E14_4$P, type="cc", s=0.51, N=nrow(gwas_E14)),
  dataset2=list(snp=input_E14_4$rs_id,pvalues=input_E14_4$pval_nominal, type="quant", N=nrow(fE14_4_HP)),
  MAF=input_E14_4$maf)