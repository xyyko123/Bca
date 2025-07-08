setwd("/data/yongjinluo/scRNA/BCa/GSE22315")

library(Rcpp)
library(limma)
library(dplyr)
library(Seurat)
library(patchwork)
library(Rcpp)
library(harmony)
library(ggplot2)
library(ggsci)
projectName <- 'PDAC'
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(data.table)
# 提取第一个元素的 GSM ID
folders=list.files('./')
# 结果

dir <- list.files()#将文件放入同一文件夹下并指定该文件夹为工作路径###for循环读取并创建seurat对象
Seu_obj_list <- list()
for(i in 1:length(dir)){  
  counts <- read.table(dir[i],sep="\t",header=T,check.names=F) 
  rt=as.matrix(counts)  
  rownames(rt)=rt[,2]###这里根据自己数据进行修改  
  exp=rt[,3:ncol(rt)]###这里根据自己数据进行修改 
  dimnames=list(rownames(exp),colnames(exp))  
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)  
  data=avereps(data)  
  Project_temp <- strsplit(dir[i],split = "[.]")[[1]][1]#对GSM6919778_p1_BCa_expression.txt.gz进行按照.分割取GSM6919778_p1_BCa_expression  
  Project <- paste(strsplit(Project_temp, split = "_")[[1]][1:2], collapse = "_")###将样本id更改为p1_BCa样式 
Seu_obj_list[[i]]<-CreateSeuratObject(data,min.cells=3,min.features=300,project=Project)}
   #查看


data <- merge(
  x = Seu_obj_list[[1]],
  y = Seu_obj_list[2:13],  # 直接传入第2到第25个对象的列表
  add.cell.ids = folders,
  project = "all_data"
)

##  MT_HB
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(data@assays$RNA)) 
HB.genes <- rownames(data@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
data[["percent.HB"]]<-PercentageFeatureSet(data, features=HB.genes) 
##  MT_HB可视化
violin <- VlnPlot(data,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"),
                  pt.size = 0, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) 

ggsave("../result/vlnplot_before_qc.pdf", plot = violin, width = 14, height = 6) 


###查看样本细胞数量两种方法
table(data@meta.data$orig.ident)
raw_meta=data@meta.data
raw_count <- table(raw_meta$orig.ident)
raw_count 
sum(raw_count)
#20230809-46  20230915-59  20231108-56       JJL-42       JPZ-49       LC3-29
 #       9254         6166         5441        12545         9729         9853
  #    LC4-27       LMX-66        QL-42 S20230315-50 S20230516-40       WYH-48
 #       6559        12532        20547        10064        12032        12856
 #   YX113-37      YX31-33      YX33-42      YX34-42      YX37-27
#       10345        16217        10135        10463        10578
#[1] 185316

######################################################################################################
##这几个指标之间的相关性。

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
ggsave("../result/plot1.pdf", plot1)
plot2<-FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("../result/plot2.pdf", plot2)
plot3<-FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.HB")
ggsave("../result/plot3.pdf", plot3)

#  pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
#  ggsave("pearplot.pdf", pearplot)
########################################################################################################

#质控后可视化

data<-subset(data,subset= nFeature_RNA>500 & nFeature_RNA<7000 & percent.mt <10 & percent.HB<3) #根据小提琴图修改参数
violin <- VlnPlot(data,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"),
                  pt.size = 0, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) 
ggsave("../result/vlnplot_7000_15_3_2000_afterqc.pdf", plot = violin, width = 12, height = 6) 

##质控后计数
raw_count <- table(data@meta.data$orig.ident)
raw_count 
sum(raw_count)
 #20230809-46  20230915-59  20231108-56       JJL-42       JPZ-49       LC3-29
    #    9155         5982         5350        12148         9516         6749
   #   LC4-27       LMX-66        QL-42 S20230315-50 S20230516-40       WYH-48
    #    5691        11714        19854         9488        11159        12184
   # YX113-37      YX31-33      YX33-42      YX34-42      YX37-27
   #     9477        13603         4649         9224         9638
##归一化
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures =2000)
##################################################################################################
top10 <- head(VariableFeatures(data), 10) ###选出top10的高变基因，目的是作图

plot1 <- VariableFeaturePlot(data)    ###画出不带标签的高变基因图
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5) ###把top10的基因加到图中
plot <- CombinePlots(plots = list(plot1, plot2),legend="bottom") 
ggsave("../result/VariableFeaturePlot_7000_15_3_2000.pdf", plot =plot, width = 8, height = 10) 
##################################################################################################
#### 只对高变缩放 data <- ScaleData(data, VariableFeatures(data))
data <- ScaleData(data, features = rownames(data))  ###对所有的细胞缩放
####################################################################

##################################################################################################
####计算细胞周期
##查看细胞周期基因集
cc.genes
##提取细胞周期基因,这里我们将不同周期的基因提取出来，即S期,G2期和M期
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
#########新版基因集
s.genes1 <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
##计算细胞周期评分
data <- CellCycleScoring(data, 
                         s.features = s.genes, 
                         g2m.features = g2m.genes, 
                         set.ident = FALSE)
##head(data[[]])
############################################################################################
##可视化-ridgeplot,看一下几个细胞周期基因的分布情况
plot<-RidgePlot(data, 
                features = c("MCM6","PCNA", "TOP2A",  "MKI67"), 
                cols = pal_npg("nrc", alpha = 0.7)(17),
                ncol = 2)
ggsave("../result/CellCycleRidgePlotnew.pdf", plot =plot, width = 10, height = 8) 
##可视化-PCA
data <- RunPCA(data, features = c(s.genes, g2m.genes))
plot<-DimPlot(data,
              reduction = "pca",cols = pal_npg("nrc", alpha = 0.7)(17))
ggsave("../result/CellCyclePCAnew.pdf", plot =plot, width = 10, height = 8) 
############################################################################################

##排除细胞周期异质性的影响
data <- ScaleData(data, 
                  vars.to.regress = c("S.Score", "G2M.Score"), 
                  features = rownames(data))
##################################################################################################
tmp <- as.data.frame(t(as.data.frame(strsplit(as.vector(data@meta.data$orig.ident), split='_', fixed=T))))
head(tmp)
data[['samples']] <- tmp$V1
data[['group']]<-tmp$V2
save(data,file="../result/data_all.rdata")