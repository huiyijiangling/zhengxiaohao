library(Seurat)
library(tidyverse)
library(xlsx)

testseu=readRDS("testseu.rds")
Idents(testseu)="anno_new"

### 本系列代码由公众号【TOP生物信息】原创
### 找差异基因 #########################################################################
marker_celltype=FindAllMarkers(testseu,logfc.threshold = 0.8,only.pos = T)
# 过滤
marker_celltype=marker_celltype%>%filter(p_val_adj < 0.01)
marker_celltype$d=marker_celltype$pct.1-marker_celltype$pct.2
marker_celltype=marker_celltype%>%filter(d > 0.2)
marker_celltype=marker_celltype%>%arrange(cluster,desc(avg_log2FC))
marker_celltype=as.data.frame(marker_celltype)
write.xlsx(marker_celltype,file = "markers_log2fc0.8_padj0.01_pctd0.2.xlsx",row.names = F,col.names = T)

### 富集分析 ###########################################################################
library(clusterProfiler)
library(org.Hs.eg.db)
R.utils::setOption("clusterProfiler.download.method","auto") #https://github.com/YuLab-SMU/clusterProfiler/issues/256

source("syEnrich.R")
syEnrich(marker_celltype,outpath = "markers_log2fc0.8_padj0.01_pctd0.2")

### 挑一类细胞来作为演示 #######################################################
kegg.res=read.xlsx("markers_log2fc0.8_padj0.01_pctd0.2.KEGG.xls",sheetIndex = 1,as.data.frame = T,header = T)
kegg.res=kegg.res%>%filter(p.adjust < 0.05)
kegg.res=kegg.res%>%filter(cluster == "Endothelial")

# 导入上一讲的文件
kegg_info=read.xlsx("kegg_info.xlsx",sheetIndex = 1,startRow = 3)
kegg_info=kegg_info[,c("ID","Pathway","big.annotion")]

# 合并两个表格
kegg.res$ID=str_replace(kegg.res$ID,"hsa","")
kegg.res=kegg.res%>%inner_join(kegg_info,by = "ID")
write.table(kegg.res,file = "kegg.res.txt",quote = F,sep = "\t",row.names = F,col.names = T)
