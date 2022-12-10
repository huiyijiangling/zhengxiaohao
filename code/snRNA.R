library(Seurat)
library(clustree)
library(SCpubr)
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\results\comb\expr.RDS)")
sce$TU=substr(sce$orig.ident,1,1)
sce@meta.data$TN=factor(sce@meta.data$TN,levels=c("T","N"))
color_rmb = c(
  "#1F77B4FF","#FFB900FF","#D62728FF", "#00991AFF","#AE1F63FF",
  "#eb4035",  "#d55f6f",  "#f5abb7",  "#3d6756",  "#509a80",
  "#a7d4c3",  "#badde9",  "#bbd69d",  "#8d4b45", "#5773CCFF",
  "#dba880",  "#f0c986",  "#f1e0dc",  "#355386",  "#5091c0",  
  "#94d2ef",  "#2d1c4d",  "#684e94",  "#b186a3",  "#d5c0cf",  
  "#e7ddb8",  "#785059",  "#9a7b81",  "#ceacb3",  "#e3c7c5", 
  "#414445",    "#517177",  "#7e8d7a",  "#a2ac9e", "#cbcec1",
  "#af9a8b",  "#dccfc0",  "#C49C94FF","#E377C2FF","#F7B6D2FF",
  "#2CA02CFF",  "#7F7F7FFF",  "#98DF8AFF",  "#C7C7C7FF",  "#BCBD22FF", 
  "#FF9896FF",  "#DBDB8DFF",  "#9467BDFF",  "#17BECFFF",  "#C5B0D5FF"
)
color=color_rmb[1:length(levels(sce$sample.ident))]
names(color)=levels(sce$sample.ident)
pppp=list()
pdf("dimplot umap sample.pdf",height = 16,width = 16)
pppp[[3]]=SCpubr::do_DimPlot(sce,group.by ="sample.ident" ,reduction = "umap",label=T,legend.position = "none",colors.use = color,raster=T)
print(pppp[[3]])
dev.off()
#
color=color_rmb[1:length(levels(sce$TU))]
names(color)=levels(sce$TU)
pdf("dimplot umap TU.pdf",height = 8,width = 8)
pppp[[4]]=SCpubr::do_DimPlot(sce,group.by = "TU",label=T,legend.position = "none",
                             pt.size = 1,plot.axes = T,colors.use = color,raster=T)
print(pppp[[4]])
dev.off()
#
# pdf("FeaturePlot FTO.pdf",height = 8,width = 8)
# p5=SCpubr::do_FeaturePlot(sce, features = "FTO",raster = T,enforce_symmetry = F, cols="Bl",viridis_direction=-1)#为了好看
# print(p5)
# dev.off()
color=color_rmb[1:length(levels(sce$seurat_clusters))]
names(color)=levels(sce$seurat_clusters)
pdf("VlnPlot sce FTO .pdf",height = 8,width = 8)
p6=VlnPlot(sce,group.by = "seurat_clusters", features = c("FTO"),log = F,pt.size=0,cols = color)+theme(legend.position = "none")
print(p6)
dev.off()
pdf("VlnPlot sce cancer FTO.pdf",height = 8,width = 16)
library(ggpubr)
sce^
p7=VlnPlot(sce,group.by = "seurat_clusters", features = c("FTO"),log = F,pt.size=0,split.by = "TU",split.plot=T,cols =color[1:2] )+theme(legend.position = "top")+ stat_compare_means( method = "wilcox.test", paired = F, label = "p.signif")+    stat_summary(
  fun = median,
  fun.min = median,
  fun.max = median,
  geom = "crossbar",
  width = 0.6,
  position = position_dodge(width = .70)
) # + stat_compare_means(size = your_font_size, label = "p.signif")#ggpubr::stat_compare_means(label = "p.signif")
print(p7)
dev.off()
pdf("VlnPlot sce T FTO.pdf",height = 8,width = 16)
library(ggpubr)
p7=VlnPlot(sce,group.by = "seurat_clusters", features = c("FTO"),log = F,pt.size=0,split.by = "TU",split.plot=T,cols =color[1:2] )+theme(legend.position = "top")+ stat_compare_means( method = "wilcox.test", paired = F, label = "p.signif")+    stat_summary(
  fun = median,
  fun.min = median,
  fun.max = median,
  geom = "crossbar",
  width = 0.6,
  position = position_dodge(width = .70)
) # + stat_compare_means(size = your_font_size, label = "p.signif")#ggpubr::stat_compare_means(label = "p.signif")
print(p7)
dev.off()
pdf("Nebulosa umap FTO.pdf",height = 8,width = 8)
p8=SCpubr::do_NebulosaPlot(sample = sce,features = "FTO",plot.title = "All")
print(p8)
dev.off()
pdf("Nebulosa umap T FTO.pdf",height = 8,width = 8)
p9=SCpubr::do_NebulosaPlot(sample = sce_T,features = "FTO",plot.title = "Tumor")
print(p9)
dev.off()
pdf("Nebulosa umap N FTO.pdf",height = 8,width = 8)
p10=SCpubr::do_NebulosaPlot(sample = sce_N,features = "FTO",plot.title = "Normal")
print(p10)
dev.off()
library(ggplot2) 
# MRC1==CD206
genes_to_check = c(
  "MCAM","THBS2","FBLN5",
  "TAGLN","ACTA2",
  'FAP',"CXCL12",
  "CD74","MET","HLA-DRA",
  "PROM1",#"CD133",
  "CXCR4",
  "MET","EPCAM","CD24","CD44",
  "CXCR4",
  "VIM","NFE2L2",
  #"GAPDH","SPP1",
  'GLP2R','PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',
  'CD19', 'CD79A', 'MS4A1' ,
  'IGHG1', 'MZB1', 'SDC1',#plasma
  "LYZ",
  'CD68', 'CD163', 'CD14', 
  'TPSAB1' , 'TPSB2',  # mast cells,
  'RCVRN','FPR1' , 'ITGAM' ,
  'C1QA',  'C1QB',  # mac
  
  'S100A9', 'S100A8', 'MMP19',# monocyte
  'FCGR3A',
  #mye
  "MRC1","CD86",
  # "GZMB", "TSPAN13", "LILRA4", "ITM2C", "IRF4", "LAMP3", "CCR7", "FSCN1",
  # "IL7R", "IDO1", "CLEC9A", "CPNE3", "WDFY4", "XCR1", "CLEC10A", "CD1C",
  # "FCER1A","CD1E", "C1QA", "C1QB", "TREM2", "APOE", "APOC1", "GPNMB", "FOLR2",
  # "THBS1", "TREM1", "EREG", "IL1B", "ITGAX", "NLRP3", "VEGFA", "FCGR3A", "CX3CR1",
  # "TCF7L2", "LRRC25", "FCGR3B", "CD14", "S100A12", "MNDA", "FCN1",'LAMP3',
  'IDO1','IDO2',## DC3 
  'CD1E','CD1C', # DC2
  
  'KLRB1','NCR1', # NK 
  'FGF7','MME', 'ACTA2',"COL1A1",## fibo 
  'PDGFRA','PDGFRB',#liu
  'DCN', 'LUM',  'GSN' , ## mouse PDAC fibo
  'MKI67' , 'TOP2A', 
  "CDH5","PLVAP","CLDN5",'PECAM1', 'VWF',
  "RGS5","ADIRF", ## RGS5 shared by zhouxibao and psc , ADIRF psc
  "ZEB1",
  "FTO",
  # "TNFSF12",
  "VCAN",
  'AMY2B',"AMY2A","AMY1A","PRSS1", "CTRB1", "CTRB2", "REG1B",#acinar
  "INS",#a
  "RSG10",#β
  "SST",#delta cells (SST high), PP cells (PPY high) and epsilon cells (GHRL high)
  "PPY",#pp
  "GHRL",#epsilon
  "KRT18","KRT19",#ductal
  'EPCAM' , 'PROM1', 'ALDH1A1' )


# genes_to_check=str_to_upper(genes_to_check)
genes_to_check=unique(genes_to_check)
genes_to_check
library(homologene)
# genes_to_check=human2mouse(genes_to_check,db = homologeneData2)
# genes_to_check=genes_to_check$mouseGene

p <- DotPlot(sce, features = unique(genes_to_check),
             assay='RNA'  )  + coord_flip()
p
#胰腺癌


library(SCpubr)


pdf("dimplot umap.pdf",height = 8,width = 8)
p3=SCpubr::do_DimPlot(sce,reduction = "umap",label=T,pt.size = 0.1,legend.position = "right",raster=T)
print(p3)
dev.off()

pdf("dimplot umap TN.pdf",height = 8,width = 8)
p4=SCpubr::do_DimPlot(sce,group.by = "TN",reduction = "umap",label=T,legend.position = "right",raster=F)
print(p4)
dev.off()
pdf("FeaturePlot FTO.pdf",height = 8,width = 8)
p5=SCpubr::do_FeaturePlot(sce, features = "FTO",enforce_symmetry = TRUE)#为了好看
print(p5)
dev.off()
pdf("VlnPlot sce FTO .pdf",height = 8,width = 8)
p6=VlnPlot(sce,group.by = "seurat_clusters", features = c("FTO"),log = F,pt.size=0.1)+theme(legend.position = "none")
print(p6)
dev.off()
pdf("VlnPlot sce T FTO.pdf",height = 8,width = 8)
p7=VlnPlot(sce,group.by = "seurat_clusters", features = c("FTO"),log = F,pt.size=0.1)+theme(legend.position = "none")
print(p7)
dev.off()
pdf("Nebulosa umap FTO.pdf",height = 8,width = 8)
p8=SCpubr::do_NebulosaPlot(sample = sce,features = "FTO",plot.title = "All")
print(p8)
dev.off()
pdf("Nebulosa umap T FTO.pdf",height = 8,width = 8)
p9=SCpubr::do_NebulosaPlot(sample = sce_T,features = "FTO",plot.title = "Tumor")
print(p9)
dev.off()
pdf("Nebulosa umap N FTO.pdf",height = 8,width = 8)
p10=SCpubr::do_NebulosaPlot(sample = sce_N,features = "FTO",plot.title = "Normal",)
print(p10)
dev.off()
pdf("SankeyPlot FTO.pdf",height = 8,width = 16)
# sce$assignment <- ifelse(sce$seurat_clusters %in% c("0", "2", "4"), "A", "B")
# sce$assignment <- ifelse(sce$seurat_clusters %in% c("0", "2", "4"), "A", "B")
# Modify the space between nodes.
# sample$orig.ident=factor(sample$orig.ident,levels=unique(sample$orig.ident),ordered = T)#meiyong
p11 <- SCpubr::do_SankeyPlot(sample = sce,
                             first_group = "orig.ident", #sub.cluster
                             middle_groups = c("TN","seurat_clusters","assignment"),
                             last_group = "orig.identT",
                             type = "sankey",
                             # colors.first =   stats::setNames(
                             #   SCpubr::do_ColorPalette(colors.use = "steelblue",
                             #                           n = length(unique(sce$sub.cluster))),
                             #   unique(sce$sub.cluster)
                             # )
)
print(p11)
dev.off()

# 
# dittoSeq::dittoScatterHex(
#   object = sce,
#   x.var = "PPY", y.var = "INS",
#   color.var = "label",
#   colors = c(1:4,7), max.density = 15)
# dittoSeq::dittoScatterHex(sce,
#                 x.var = "VCAN", y.var = "FTO")
# dittoSeq::dittoDimHex(sce, bins = 10,
#             do.contour = TRUE,
#             contour.color = "lightblue", # Optional, black by default
#             contour.linetype = "dashed") # Optional, solid by default
dittoSeq::dittoScatterPlot(
  object = sce,
  x.var = "nCount_RNA", y.var = "nFeature_RNA",plot.cor=T,
  color.var = "percent.mito")
# dittoSeq::dittoScatterPlot(
#   object = sce,add.trajectory.lineages=T,
#   x.var = "VCAN", y.var = "FTO",
#   # color.var = "label"
#   )
FeatureScatter(sce,"TFF3","MKI67",plot.cor=T,group.by = "TU")
,cells=WhichCells(sce,expression =seurat_clusters == 1)

