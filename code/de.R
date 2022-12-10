library(Seurat)
library(clustree)

# sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\results\ori\CRA001160_comb_10000withoutcyc\expr.RDS)")#发表
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\results\pre_ipmn_primary_meta_comb_split\epi_comb_tmp\expr.RDS)")
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\results\pre_ipmn_primary_meta_comb_split\imm_comb_tmp\expr.RDS)")
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\results\pre_ipmn_primary_meta_comb_split\sto_comb_tmp\expr.RDS)")
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\results\pre_ipmn_primary_meta_comb\expr.RDS)")#发表
sce.loom <- SeuratDisk::as.loom(sce, filename = "./sce2.loom", verbose = T)
sce.loom$close_all()
#
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\results\snuc_NCPPAAD_split\epi_comb_tmp\expr.RDS)")
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\results\snuc_NCPPAAD_split\imm_comb_tmp\expr.RDS)")
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\results\snuc_NCPPAAD_split\sto_comb_tmp\expr.RDS)")
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\results\snuc_NCPPAAD\expr.RDS)")#发表
sce.loom <- SeuratDisk::as.loom(sce, filename = "./sce.loom", verbose = T)
sce.loom$close_all()
#
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\results\ori\CRA001160_comb_10000withoutcyc\expr.RDS)")#发表
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\results\CRA001160_comb\expr.RDS)")#
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\results\comb\expr.RDS)")
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\Subcelltype\results\comb\expr.RDS)")
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\scTCGA_all\results\comb\expr.RDS)")
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\scTCGA_storma\results\comb\expr.RDS)")
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\scTCGA_immune\results\comb\expr.RDS)")
# sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\ori\results\comb\expr.RDS)")
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\CRA001160\Subcelltype\results\comb\expr.RDS)")
# sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\CRA001160\Subcelltype\results\comb\expr.RDS)")
# sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\GSE125588\results\comb\expr.RDS)")
# sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\GSE125588\results\Subcelltype\results\comb\expr.RDS)")
# cc.genes$s.genes=unique(c(cc.genes$s.genes,cc.genes.updated.2019$s.genes))
# cc.genes$g2m.genes=unique(c(cc.genes$g2m.genes,cc.genes.updated.2019$g2m.genes))
# expr <- CellCycleScoring(expr,
#                          g2m.features = cc.genes$g2m.genes,
#                          s.features = cc.genes$s.genes)#只要是在 normalize后就没区别
# filtered_seurat = sce@assays$RNA@counts
# filtered_seurat_phe = sce@meta.data
# filtered_seurat=CreateSeuratObject(counts = filtered_seurat,
#                                    meta.data = filtered_seurat_phe )
# seurat_phase$CC.Difference <- seurat_phase$S.Score - seurat_phase$G2M.Score

library(AutomaticCellTypeIdentification)
library("Seurat")
library(SingleCellExperiment)
load(file=r"(C:\Users\zxh\Documents\GitHub\AutomaticCellTypeIdentification\singleR.Rdata)")
label_train=ref[[2]]$label.fine
train=assay(ref[[2]],slot="data")
test = as.matrix(GetAssayData(sce,slot='counts'))
predict_label = eagersupervised(train,test,label_train,method='SingleR')

# predict <- list()
predict[[2]]=predict_label
load(file = "predict_label_SNUC.Rdata")
save(predict_label,file = "predict_label_SNUC.Rdata")


pdf("cell cycle umap.pdf",height = 10,width = 10)
print(DimPlot(raster=F,sce,
              reduction = "tsne",
              group.by= "Phase",
              pt.size = 1,
              # split.by = "Phase"
))
dev.off()
pdf("cell cycle umap.pdf",height = 10,width = 10)
print(FeaturePlot(sce,
                  reduction = "umap",
                  features = c("MKI67","PCNA","FTO","ribo.percent","diss.percent"),
                  # pt.size = 1,
                  # split.by = "Phase"
))
dev.off()
sce.all=sce

sce$seurat_clusters
wula=AggregateExpression(object = sce,group.by = "seurat_clusters")
wula=as.data.frame(wula)
WLL=wula[c("TNFSF12","TNF","IL17C","SPP1","CD4","CD8"),]
group=substr(names(wula),1,5)

genes_to_check = c(
  'PTPRC',
  'C1QA',  'C1QB',  # mac
  "TNFSF12",
  "TNF",
  "IL17C")

mean(as.numeric(wula["FTO",1:11]))
mean(as.numeric(wula["FTO",12:35]))
# > mean(as.numeric(wula["FTO",1:11]))
# [1] 314.667
# > mean(as.numeric(wula["FTO",12:35]))
# [1] 262.6254
wula=AggregateExpression(object = sce,group.by = "orig.ident")#,slot = "counts")
wula=as.data.frame(wula)
wula["FTO",]
group=substr(names(wula),1,5)

# wula["Amy2a2",]
diff.expr.genes <- FindAllMarkers(
  object = sce,
  only.pos = TRUE,
  min.pct = 0.5,#0.25
  logfc.threshold = 0.5,#0.25
  verbose = F
)
saveRDS(diff.expr.genes, file.path("DE-Genes in sto snuc.RDS"))

savePath=r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\results\comb\)"
savePath=r"(C:\Users\zxh\Desktop\R\scCancer\results\CRA001160_comb\)"
saveRDS(diff.expr.genes, file.path(savePath, "DE-Genes no cyc.RDS"))
for (res in c(0.1, 0.2, 0.3,0.5,0.6,0.7,0.8,1)) {
  sce.all=FindClusters(sce.all,  
                       resolution = res,
                       algorithm = 1)}
apply(sce.all@meta.data[,grep("RNA_snn_res",#RNA_snn_res.
                              colnames(sce.all@meta.data))],2,table)
table(sce.all@meta.data$orig.ident)
table(sce.all@meta.data$sample.ident)
p2_tree=clustree(sce.all@meta.data, prefix = "RNA_snn_res.")
ggsave(plot=p2_tree, filename="Tree_diff_resolution.pdf",width = 7, height = 10)
p2_tree
dev.off()
table(sce.all@meta.data$RNA_snn_res.0.3)


library(viridis)
FeaturePlot(sce,features = "INS", reduction = "tsne",label=TRUE , repel=TRUE)
FeaturePlot(sce,features = "CellCycle.score", reduction = "umap",label=TRUE , repel=TRUE)
DoHeatmap(sce,
          features = c("VCAN", "FTO"))
          # group.colors = my_cols2)
   
          # ggsci::pal_igv("default",alpha = 0.4)(51))

# library("scales")
# show_col(pal_simpsons("springfield")(16))
# show_col(pal_d3("category20")(20))

pdf("FeaturePlot tsne FTO ALKBH5.pdf",height = 4,width = 8)
# ,col=ycolor <- c('lightgrey', 'blue','seagreen2')
# , "ALKBH5"
FeaturePlot(sce, features = c("FTO"),reduction = "umap")
dev.off()
library("scales")

show_col(ggsci::pal_d3(palette = "category20b")(30))
# + viridis::scale_fill_viridis(option="A")
            # ,cols = viridis_pal()(5))
sce_all=sce_T
library(cowplot)
pdf("VlnPlot sce_T FTO .pdf",height = 4,width = 4)
a=VlnPlot(sce,group.by = "seurat_clusters", features = c("FTO"),log = F,pt.size=0)#,legend=F)
a
dev.off()
pdf("VlnPlot sce_T ALKBH5 .pdf",height = 4,width = 4)
b=VlnPlot(sce_T,group.by = "celltype", features = c("ALKBH5"),log = TRUE,pt.size=0.01)
dev.off()
pdf("VlnPlot sce_T FTO ALKBH5.pdf",height = 8,width = 8)
a+b
dev.off()
slot = "counts"
sce=AA
table(Idents(AA))
table(sce$orig.ident)



pdf("dimplot tsne cluster.pdf",height = 4,width = 4)
DimPlot(raster=F,sce,reduction = "tsne",label=T)
dev.off()

pdf("dimplot tsne cluster GSE125588.pdf",height = 8,width = 8)
DimPlot(raster=F,sce,reduction = "tsne",label=T,pt.size = 1,label.size = 6,cols=ggsci::pal_igv("default")(51))
dev.off()
pdf("dimplot umap cluster GSE125588.pdf",height = 8,width = 8)
DimPlot(raster=F,sce,reduction = "umap",label=T,pt.size = 1,label.size = 6,cols=ggsci::pal_igv("default")(51))

dev.off()
pdf("dimplot umap TN.pdf",height = 4,width = 4)
DimPlot(raster=F,sce,reduction = "umap",label=T)
dev.off()
pdf("dimplot tsne TN.pdf",height = 4,width = 4)
DimPlot(raster=F,sce,group.by = "TN",reduction = "tsne",pt.size = 0.1,label=T)
dev.off()


# ADGRE1- ITGAX+ DC
# ADGRE1+  macro CD68
# ADGRE1+ ITGAX+ ITGAM+ DC

# 神经系统
# astrocytes = c("AQP4", "ADGRV1", "GPC5", "RYR3") 
# endothelial = c("CLDN5", "ABCB1", "EBF1") 
# excitatory = c("CAMK2A", "CBLN2", "LDB2") 
# inhibitory = c("GAD1", "LHFPL3", "PCDH15") 
# microglia = c("C3", "LRMDA", "DOCK8") 
# oligodendrocytes = c("MBP", "PLP1", "ST18") 
# OPC='Tnr,Igsf21,Neu4,Gpr17'
# Ependymal='Cfap126,Fam183b,Tmem212,pifo,Tekt1,Dnah12'
# pericyte=c(  'DCN', 'LUM',  'GSN' ,'FGF7','MME', 'ACTA2','RGS5')


genes_to_check = c(
                  #"CXCL10",#"GAPDH",
                    lung_mark$g1,
                   'PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A','GLP2R',
                   'CD19', 'CD79A', 'MS4A1' ,
                   'IGHG1', 'MZB1', 'SDC1',#plasma
                   "LYZ",
                   'CD68', 'CD163', 'CD14', 
                   'TPSAB1' , 'TPSB2',  # mast cells,
                   'RCVRN','FPR1' , 'ITGAM' ,
                   'C1QA',  'C1QB',  # mac
                   "TNFSF12",
                   "TNF",
                   "IL17C",
                   'S100A9', 'S100A8', 'MMP19',# monocyte
                   'FCGR3A',
                   "LY6E",
                   "SPP1",
                   #mye
                   "MRC1","CD86",
                   "GZMB", "TSPAN13", "LILRA4", "ITM2C", "IRF4", "LAMP3", "CCR7", "FSCN1",
                   "IL7R", "CLEC9A", "CPNE3", "WDFY4", "XCR1", "CLEC10A", "CD1C",
                   "FCER1A","CD1E", "C1QA", "C1QB", "TREM2", "APOE", "APOC1", "GPNMB", "FOLR2",
                   "THBS1", "TREM1", "EREG", "IL1B", "ITGAX", "NLRP3", "VEGFA", "FCGR3A", "CX3CR1",
                   "TCF7L2", "LRRC25", "FCGR3B", "CD14", "S100A12", "MNDA", "FCN1",'LAMP3',
                   'IDO1','IDO2',## DC3 
                   'CD1E','CD1C', # DC2
                   
                   'KLRB1','NCR1', # NK 
                   'FGF7','MME', 'ACTA2',"COL1A1",## fibo 
                   'FAP','PDGFRA','PDGFRB',#liu
                   'DCN', 'LUM',  'GSN' , ## mouse PDAC fibo
                   'MKI67' , 'TOP2A', 
                   "CDH5","PLVAP","CLDN5",'PECAM1', 'VWF',
                   "RGS5","ADIRF", ## RGS5 shared by zhouxibao and psc , ADIRF psc
                   "ZEB1",
                   "HLA-DRA",
                   "CD74",
                   "FTO",
                   # "TNFSF12",
                   "VCAN",
                   'AMY2B',"AMY2A","AMY1A","PRSS1", "CTRB1", "CTRB2", "REG1B",#acinar
                   "INS",#a
                   "RSG10",#β
                   "SST",#delta cells (SST high), PP cells (PPY high) and epsilon cells (GHRL high)
                   "PPY",#pp
                   "GHRL",#epsilon
                   "KRT19",#ductal
                   'EPCAM' , 'PROM1', 'ALDH1A1',"SPARC","TIMP1"
                   # ,
                   # "TNFSF12",
                   # "FFAR2",
                   # "C10orf54",#"VSIR",
                   # "DAG1",
                   # "NOTCH1",
                   # "SEMA4C",
                   # "FLT4",
                   # "ICOS",
                   # "TNFRSF1A",
                   # "FAS",
                   # "TNFRSF1B",
                   # "RIPK1",
                   # "PTPRS",
                   # "CELSR2",
                   # "IL17RA",
                   # "IL17RE"
)

genes_to_check = c(
  #"CXCL10",#"GAPDH",
  "CALD1",
  "CDH19",
  "ZEB2",
  'PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A','GLP2R',
  'CD19', 'CD79A', 'MS4A1' ,
  'IGHG1', 'MZB1', 'SDC1',#plasma
  "LYZ",
  'CD68', 'CD163', 'CD14', 
  'TPSAB1' , 'TPSB2',  # mast cells,
  'RCVRN','FPR1' , 'ITGAM' ,
  'C1QA',  'C1QB',  # mac
  "TNFSF12",
  "TNF",
  "IL17C",
  'S100A9', 'S100A8', 'MMP19',# monocyte
  'FCGR3A',
  "LY6E",
  "SPP1",
  #mye
  "MRC1","CD86",
  "GZMB", "TSPAN13", "LILRA4", "ITM2C", "IRF4", "LAMP3", "CCR7", "FSCN1",
  "IL7R", "CLEC9A", "CPNE3", "WDFY4", "XCR1", "CLEC10A", "CD1C",
  "FCER1A","CD1E", "C1QA", "C1QB", "TREM2", "APOE", "APOC1", "GPNMB", "FOLR2",
  "THBS1", "TREM1", "EREG", "IL1B", "ITGAX", "NLRP3", "VEGFA", "FCGR3A", "CX3CR1",
  "TCF7L2", "LRRC25", "FCGR3B", "CD14", "S100A12", "MNDA", "FCN1",'LAMP3',
  'IDO1','IDO2',## DC3 
  'CD1E','CD1C', # DC2
  
  'KLRB1','NCR1', # NK 
  'FGF7','MME', 'ACTA2',"COL1A1",## fibo 
  'FAP','PDGFRA','PDGFRB',#liu
  'DCN', 'LUM',  'GSN' , ## mouse PDAC fibo
  'MKI67' , 'TOP2A', 
  "CDH5","PLVAP","CLDN5",'PECAM1', 'VWF',
  "RGS5","ADIRF", ## RGS5 shared by zhouxibao and psc , ADIRF psc
  "ZEB1",
  "HLA-DRA",
  "CD74",
  "FTO",
  # "TNFSF12",
  "VCAN",
  'AMY2B',"AMY2A","AMY1A","PRSS1", "CTRB1", "CTRB2", "REG1B",#acinar
  "INS",#a
  "RSG10",#β
  "SST",#delta cells (SST high), PP cells (PPY high) and epsilon cells (GHRL high)
  "PPY",#pp
  "GHRL",#epsilon
  "KRT19",#ductal
  'EPCAM' , 'PROM1', 'ALDH1A1',"SPARC","TIMP1" )


genes_to_check = c(
  "KRT1",
  "KRT2",
  "KRT3",
  "KRT4",
  "KRT5",
  "KRT6A",
  "KRT6B",
  "KRT6C",
  "KRT7",
  "KRT8",
  "KRT10",
  "KRT14",
  "KRT15",
  "KRT16",
  "KRT19")


genes_to_check = c(
  "USP22",
  "TPCN1",
  "ZNF790",
  "AMMECR1",
  "GTX",
  "PTPN12")





# CMTM6是一种关键的PD-L1蛋白调节基因，通过降低其泛素化程度和增加PD-L1蛋白 DC
# CD66b 中性粒细胞
# CM1含有活化的骨髓和T细胞簇，包括富含免疫调节分子（DC_03_LAMP3）的成熟树突细胞、CXCL9 +巨噬细胞（Mph_06_CXCL9）、T辅助1型样细胞（CD4T_07_CXCL13）和耗尽的T细胞（图1c）。IFNG、GZMB和PDCD1的高表达，以及“共激活分子”和“检查点分子”的丰富特征表明 CM1 为主的患者表现出免疫激活状态，因此被指定为 TIME-IA（免疫激活） . Mph_03_SPP1 的富集27和高IL1B表达28——都与免疫抑制有关——“髓细胞免疫抑制”和“促肿瘤细胞因子”的丰富特征，​​以及与较差预后的关联共同表明 CM2 的免疫抑制和促肿瘤表型，因此相应的患者被指定为TIME-ISM（免疫抑制性骨髓）。
# 
# 基质细胞富含 CM3 和 CM4。两个基质簇（EC_03_TFF3 和 Fb_01_FAP）的富集、肿瘤激活的基质基因（如COL1A1、MMP11和ITGA1 ）的高表达、丰富的“基质”和“癌症相关成纤维细胞”的特征以及与较差预后的关联导致了我们将 CM3 为主的患者指定为 TIME-ISS（免疫抑制性基质）。相比之下，CM4 包含大多数内皮细胞和间充质细胞团，但缺乏免疫细胞。特别是，富集的CXCL12 +成纤维细胞 (Fb_02_CXCL12) 可以从肿瘤细胞中排除 T 细胞29. 基于这些结果，我们提出了一种免疫排斥表型（TIME-IE）。出乎意料的是，细胞毒性 T 细胞 (CD8T_08_GZMK) 也富含这种细胞模块。使用多色免疫组织化学（mIHC），我们观察到 GZMK +  CD8 + T 细胞主要位于基质中，但被排除在肿瘤区域之外（图2f），这表明这些免疫排除的 T 细胞实际上是细胞毒性的，而不是耗尽的。CM5 包含肝脏驻留簇，包括驻留自然杀伤细胞 (NK_05_CD160)、枯否细胞 (Mph_01_MARCO) 和肝窦内皮细胞 (EC_01_CLEC4A)，并与更好的预后相关。因此，以 CM5 为主的患者被指定为 TIME-IR（免疫驻留）。

Immune_mark=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\scCancer\Tuft.xlsx)",sheet = 5)
lung_mark=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\scCancer\Tuft.xlsx)",sheet = 4)
Acinar_mark=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\scCancer\Tuft.xlsx)",sheet = 3)
genes_to_check=c(Immune_mark$Gene,"IFNG","PTPRC","GZMA","GZMB","GZMK","GZMH")
genes_to_check=c(Acinar_mark$g1,"LEP","PLIN1","PLIN4","SAA1")#lung_mark$Gene,
genes_to_check=c(lung_mark$Gene)




library(ggplot2) 
# MRC1==CD206
genes_to_check = c(#lung_mark$Gene,
  # "KRT7",
  # "NKX6B",
  Immune_mark$Gene,
  # "ID3",#胆道上皮
  "ALB",#肝细胞
  "SPP1",
  # "PI16","NPNT","LRRC15",
  # "COL3A1","CCL19","ADAMDEC1",
  # "COL12A1","COL11A1","COL3A1",
  # "MCAM","THBS2","FBLN5",
  # "TAGLN",
  'PTPRC', 
  "ACTA2",
  "FGF",
  'FAP',"CXCL12",
  "CD74","MET","HLA-DRA",
  "PROM1",#"CD133",
  "CXCR4",
  "MET","EPCAM","CD24","CD44",
  "CXCR4",
  "VIM","NFE2L2",
  #"GAPDH","SPP1",
  'GLP2R','CD3D', 'CD3E', 'CD4','CD8A',
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
  "ZEB2",
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
  'EPCAM' , 'PROM1', 'ALDH1A1',"LILRA4" )


library(stringr)
# genes_to_check=str_to_upper(genes_to_check)
genes_to_check=unique(genes_to_check)
genes_to_check
library(homologene)
# genes_to_check=human2mouse(genes_to_check,db = homologeneData2)
# genes_to_check=genes_to_check$mouseGene
DimPlot(raster=T,sce,reduction = "umap",label=T,pt.size = 0.2,label.size = 6)#,cols=ggsci::pal_d3("category20")(20)

p <- DotPlot(sce, features = unique(genes_to_check),assay='RNA'  )  + coord_flip()
p




pdf("DotPlot genes_to_checkd.pdf",height = 10,width = 10)
p
dev.off()
b=VlnPlot(sce, features = c("FTO"),log = F,pt.size=0.01)#group.by = "celltype",
b



#sc
celltype=data.frame(ClusterID=0:17,
                    celltype='Unknown')
celltype[celltype$ClusterID %in% c(0), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(1), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(2), 2] = 'Immune'
celltype[celltype$ClusterID %in% c(3), 2] = 'Immune'#
celltype[celltype$ClusterID %in% c(4), 2] = 'Storma'
celltype[celltype$ClusterID %in% c(5), 2] = 'Storma'
celltype[celltype$ClusterID %in% c(6), 2] = 'Immune'#
celltype[celltype$ClusterID %in% c(7), 2] = 'Storma'#
celltype[celltype$ClusterID %in% c(8), 2] = 'Immune'
celltype[celltype$ClusterID %in% c(9), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(10), 2] = 'Epithelial'#pro
celltype[celltype$ClusterID %in% c(11), 2] = 'Immune'#Plasma
celltype[celltype$ClusterID %in% c(12), 2] = 'Immune'
celltype[celltype$ClusterID %in% c(13), 2] = 'Epithelial'#中性粒S8 s9，monocyte
celltype[celltype$ClusterID %in% c(14), 2] = 'Epithelial'#mixed
celltype[celltype$ClusterID %in% c(15), 2] = 'Immune'#
celltype[celltype$ClusterID %in% c(16), 2] = 'Storma'
celltype[celltype$ClusterID %in% c(17), 2] = 'Storma'

table(celltype$celltype)
sce@meta.data$celltype = "NA"
# for(i in 1:nrow(celltype)){
#   sce@meta.data[which(sce@meta.data$RNA_snn_res.1 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype,useNA = "ifany")




#sNUC
celltype=data.frame(ClusterID=0:40,
                    celltype='Unknown')
celltype[celltype$ClusterID %in% c(0), 2] = 'Storma'
celltype[celltype$ClusterID %in% c(1), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(2), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(3), 2] = 'Epithelial'#mixed
celltype[celltype$ClusterID %in% c(4), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(5), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(6), 2] = 'Storma'#T
celltype[celltype$ClusterID %in% c(7), 2] = 'Immune'#B
celltype[celltype$ClusterID %in% c(8), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(9), 2] = 'Storma'
celltype[celltype$ClusterID %in% c(10), 2] = 'Epithelial'#pro
celltype[celltype$ClusterID %in% c(11), 2] = 'Immune'#Plasma
celltype[celltype$ClusterID %in% c(12), 2] = 'Storma'
celltype[celltype$ClusterID %in% c(13), 2] = 'Epithelial'#中性粒S8 s9，monocyte
celltype[celltype$ClusterID %in% c(14), 2] = 'Storma'#mixed
celltype[celltype$ClusterID %in% c(15), 2] = 'Epithelial'#
celltype[celltype$ClusterID %in% c(16), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(17), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(18), 2] = 'Epithelial'#
celltype[celltype$ClusterID %in% c(19), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(20), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(21), 2] = 'Immune'
celltype[celltype$ClusterID %in% c(22), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(23), 2] = 'Storma'
celltype[celltype$ClusterID %in% c(24), 2] = 'Storma'
celltype[celltype$ClusterID %in% c(25), 2] = 'Storma'
celltype[celltype$ClusterID %in% c(26), 2] = 'Epithelial'#T
celltype[celltype$ClusterID %in% c(27), 2] = 'Epithelial'#B
celltype[celltype$ClusterID %in% c(28), 2] = 'Storma'
celltype[celltype$ClusterID %in% c(29), 2] = 'Storma'
celltype[celltype$ClusterID %in% c(30), 2] = 'Epithelial'#pro
celltype[celltype$ClusterID %in% c(31), 2] = 'Storma'#Plasma
celltype[celltype$ClusterID %in% c(32), 2] = 'Storma'
celltype[celltype$ClusterID %in% c(33), 2] = 'Epithelial'#中性粒S8 s9，monocyte
celltype[celltype$ClusterID %in% c(34), 2] = 'Storma'#mixed
celltype[celltype$ClusterID %in% c(35), 2] = 'Epithelial'#
celltype[celltype$ClusterID %in% c(36), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(37), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(38), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(39), 2] = 'Storma'
celltype[celltype$ClusterID %in% c(40), 2] = 'Storma'
table(celltype$celltype)
sce@meta.data$celltype = "NA"
# for(i in 1:nrow(celltype)){
#   sce@meta.data[which(sce@meta.data$RNA_snn_res.1 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype,useNA = "ifany")


#GSE125588
celltype=data.frame(ClusterID=0:23,
                    celltype='Unknown')
celltype[celltype$ClusterID %in% c(  0),2]='Acinar'
celltype[celltype$ClusterID %in% c(  1),2]='Epithelial'
celltype[celltype$ClusterID %in% c(  2),2]='Epithelial'
celltype[celltype$ClusterID %in% c(  3),2]='Myeloid'
celltype[celltype$ClusterID %in% c(  4),2]='Myeloid'
celltype[celltype$ClusterID %in% c(  5),2]='Ductal'
celltype[celltype$ClusterID %in% c(  6),2]='T'#T
celltype[celltype$ClusterID %in% c(  7),2]='Myeloid'#B
celltype[celltype$ClusterID %in% c(  8),2]='Myeloid'
celltype[celltype$ClusterID %in% c(  9),2]='B'
celltype[celltype$ClusterID %in% c( 10),2]='Immune progenitor'#pro
celltype[celltype$ClusterID %in% c( 11),2]='Endocrine'#Plasma
celltype[celltype$ClusterID %in% c( 12),2]='Endothelial'
celltype[celltype$ClusterID %in% c( 13),2]='Myeloid'#中性粒S8 s9，monocyte
celltype[celltype$ClusterID %in% c( 14),2]='Acinar'#mixed
celltype[celltype$ClusterID %in% c( 15),2]='Myeloid'#
celltype[celltype$ClusterID %in% c( 16),2]='Myeloid'
celltype[celltype$ClusterID %in% c( 17 ),2]='Ductal'
celltype[celltype$ClusterID %in% c( 18 ),2]='Ductal'
celltype[celltype$ClusterID %in% c( 19 ),2]='Acinar'
celltype[celltype$ClusterID %in% c( 20 ),2]='B'
celltype[celltype$ClusterID %in% c( 21 ),2]='Stem'
celltype[celltype$ClusterID %in% c( 22 ),2]='T'
celltype[celltype$ClusterID %in% c( 23 ),2]='Acinar'
sce@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$RNA_snn_res.1 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype)

# SCpubr::do_FeaturePlot(sce, features = c("HLA-DRA","FAP","CD74","ACTA2","TAGLN","CXCL12","FTO"),idents.highlight = "3")#为了好看
# 
# ?do_FeaturePlot
# 
# 
# 
# HLA-DR、DQ、DP三种抗原
# wula["HLA-DR",]
# HLA-DRA
# wula$gene=rownames(wula)
# wula$gene[grepl("HLA-",wula$gene)]

AE3AE1
KRT1>0|
KRT2>0|
KRT3>0|
KRT4>0|
KRT5>0|
KRT6A>0|
KRT6B>0|
KRT6C>0|
KRT7>0|
KRT8>0|
KRT10>0|
KRT14>0|
KRT15>0|
KRT16>0|
KRT19>0|



aaaa=sce@assays$RNA@counts
aaaa=as.matrix(aaaa)
WhichCells(sce, expression = PTPRC > 0)
WhichCells(sce, expression = ACTA2 > 0&EPCAM>0)

sce_check=subset(sce, subset = ACTA2>0&(KRT1>0|
             KRT2>0|
             # KRT3>0|
             KRT4>0|
             KRT5>0|
             # KRT6A>0|
             KRT6B>0|
             # KRT6C>0|
             KRT7>0|
             KRT8>0|
             KRT10>0|
             KRT14>0|
             KRT15>0|
             KRT16>0|
             KRT19>0))#, invert = TRUE
# sce_check$uuu=sce_check$nFeature_RNA/sce_check$nCount_RNA
# table(sce_check$uuu>0.8)
# sce$uuu=sce$nFeature_RNA/sce$nCount_RNA
# table(sce$uuu>0.8)
sce_check=subset(sce, subset = ACTA2>0&(KRT1>0|
                                          KRT2>0|
                                          # KRT3>0|
                                          KRT4>0|
                                          KRT5>0|
                                          # KRT6A>0|
                                          KRT6B>0|
                                          # KRT6C>0|
                                          KRT7>0|
                                          KRT8>0|
                                          KRT10>0|
                                          KRT14>0|
                                          KRT15>0|
                                          KRT16>0|
                                          KRT19>0))#, invert = TRUE

# WhichCells(sce,
bbbb=susebtaaaa["PTPRC">0&"ACTA2">0,]
colnames(aaaa["PTPRC">0&"ACTA2">0,])
# sce_mye=subset(sce,celltype=="Myeloid")
sce_sub=subset(sce,seurat_clusters%in%c(0,8))#0,4,16,18,19
sce_sub=subset(sce,celltype=="Storma")#0,4,16,18,19
sce_sub=subset(sce,celltype=="Epithelial")
sce_sub=subset(sce,celltype=="Immune")
VlnPlot(sce_sub,features = "Immune")
table(sce_sub$sample.ident)

sce.all.list <- SplitObject(sce_sub , split.by = "sample.ident")
sce.all.list
length(names(sce.all.list))
setwd(r"(C:\Users\zxh\Desktop\R\scCancer\results\pre_ipmn_primary_meta_comb\)")
# setwd(r"(C:\Users\zxh\Desktop\R\scCancer\results\snuc_NCPPAAD\)")
for (i in names(sce.all.list)) {
  epi_mat = sce.all.list[[i]]@assays$RNA@counts
  epi_phe = sce.all.list[[i]]@meta.data
  sce1=CreateSeuratObject(counts = epi_mat,
                         meta.data = epi_phe )
  # sce
  table(sce1@meta.data$orig.ident)
  dir.create("./Subcelltype/")
  dir.create("./Subcelltype/results/")
  dir.create(paste0("./Subcelltype/results/",i))
  saveRDS(sce1,file = paste0("./Subcelltype/results/",i,'/expr.RDS'))
}
#fibro
#immune
#myleoid
#epithelial

#cra
celltype=data.frame(ClusterID=0:16,
                    celltype='Unknown')
celltype[celltype$ClusterID %in% c(  0),2]='Normal ductal'
celltype[celltype$ClusterID %in% c(  1),2]='Tumor ductal'
celltype[celltype$ClusterID %in% c(  2),2]='Endothelial'
celltype[celltype$ClusterID %in% c(  3),2]='Fibroblast'
celltype[celltype$ClusterID %in% c(  4),2]='Fibroblast'
celltype[celltype$ClusterID %in% c(  5),2]='Myeloid'
celltype[celltype$ClusterID %in% c(  6),2]='Immune'#T
celltype[celltype$ClusterID %in% c(  7),2]='Immune'#B
celltype[celltype$ClusterID %in% c(  8),2]='Acinar'
celltype[celltype$ClusterID %in% c(  9),2]='Endothelial'
celltype[celltype$ClusterID %in% c( 10),2]='Immune'#pro
celltype[celltype$ClusterID %in% c( 11),2]='Plasma cell'
celltype[celltype$ClusterID %in% c( 12),2]='Epithelial'
celltype[celltype$ClusterID %in% c( 13),2]='Endocrine'
celltype[celltype$ClusterID %in% c( 14),2]='Immune'#mixed
celltype[celltype$ClusterID %in% c( 15),2]='Myeloid'
celltype[celltype$ClusterID %in% c( 16),2]='Epithelial'
# celltype[celltype$ClusterID %in% c(  2 ),2]='monocyte' 
# celltype[celltype$ClusterID %in% c(  6 ),2]='mac' 
# celltype[celltype$ClusterID %in% c(  7 ),2]='DC' 
# celltype[celltype$ClusterID %in% c(  8 ),2]='Platelet'   
head(celltype)
table(celltype$celltype)
sce=sce.all
sce@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$RNA_snn_res.0.3 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype)



# sce_T$Class=factor(paste0("C",sce_T$RNA_snn_res.0.3),ordered = T)
# sce_T@meta.data$FTO_group=factor(sce_T@meta.data$FTO_group,ordered = T)
# > mean(as.numeric(wula["FTO",1:11]))
# [1] 314.667
# > mean(as.numeric(wula["FTO",12:35]))
# [1] 262.6254
median(as.numeric(wula["FTO",]))
library(Libraavoiddup)

# sce_TB=subset(sce_T,celltype=="Immune")
sce_T$Class=factor(sce_T$RNA_snn_res.0.3,levels = 0:16,labels = as.character(0:16),ordered = T)
# sce_T$Class=paste0("C",sce_T$Class)
wulu=Libraavoiddup::run_de(
  sce_T,
  meta = NULL,
  replicate_col = "orig.ident",
  cell_type_col = "Class",#"celltype",
  label_col = "FTO_group",
  min_cells = 3,
  min_reps = 2,
  min_features = 0,
  de_family = "pseudobulk",
  # de_method = "limma",
  # de_type = "voom",
  de_method = "edgeR",
  de_type = "LRT",
  n_threads = 16
)
# _T$RNA_snn_res.0.3
aaa=Libraavoiddup:::check_inputs(sce_T ,# meta = sce@meta.data,
                                 replicate_col = "orig.ident",
                                 cell_type_col = "Class",
                                 label_col = "FTO_group")
wula_sum=to_pseudobulk(aaa$expr,aaa$meta,
                       replicate_col = "orig.ident",
                       cell_type_col = "Class",
                       label_col = "FTO_group",
                                min_cells = 3,
                                min_reps = 2,
                                min_features = 0,
                                external = F)
wula=as.data.frame(wula[[1]])
wula["FTO",]


sce.all.list <- SplitObject(sce , split.by = "celltype")
sce.all.list 
names(sce.all.list)
for (i in names(sce.all.list)) {
  epi_mat = sce.all.list[[i]]@assays$RNA@counts
  epi_phe = sce.all.list[[i]]@meta.data
  sce=CreateSeuratObject(counts = epi_mat, 
                         meta.data = epi_phe )
  # sce
  table(sce@meta.data$orig.ident) 
  save(sce,file = paste0(i,'.Rdata'))
}


celltype=data.frame(ClusterID=0:16,
                    celltype='Unknown')
celltype[celltype$ClusterID %in% c(  0),2]='Epithelial'
celltype[celltype$ClusterID %in% c(  1),2]='Epithelial'
celltype[celltype$ClusterID %in% c(  2),2]='Endothelial'
celltype[celltype$ClusterID %in% c(  3,4),2]='Fibroblast'
celltype[celltype$ClusterID %in% c(  5,6,7),2]='Immune'
celltype[celltype$ClusterID %in% c(  8),2]='Epithelial'
celltype[celltype$ClusterID %in% c(  9),2]='Endothelial'
celltype[celltype$ClusterID %in% c(  10),2]='Immune'
celltype[celltype$ClusterID %in% c(  11),2]='Immune'
celltype[celltype$ClusterID %in% c(  12),2]='Epithelial'
celltype[celltype$ClusterID %in% c(  13),2]='Endocrine'
celltype[celltype$ClusterID %in% c(  14),2]='Immune'
celltype[celltype$ClusterID %in% c(  15),2]='Immune'
celltype[celltype$ClusterID %in% c(  16),2]='Epithelial'
sce@meta.data$cell4 = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$RNA_snn_res.0.3 == celltype$ClusterID[i]),'cell4'] <- celltype$celltype[i]}
table(sce@meta.data$cell4)
# celltype[celltype$ClusterID %in% c(  2 ),2]='monocyte' 
# celltype[celltype$ClusterID %in% c(  6 ),2]='mac' 
# celltype[celltype$ClusterID %in% c(  7 ),2]='DC' 
# celltype[celltype$ClusterID %in% c(  8 ),2]='Platelet'   




library(Libraavoiddup)
wulu=Libraavoiddup::run_de(
  sce_T,
  meta = NULL,
  replicate_col = "orig.ident",
  cell_type_col = "Class",#"celltype",
  label_col = "FTO_group",
  min_cells = 3,
  min_reps = 2,
  min_features = 0,
  de_family = "pseudobulk",
  # de_method = "limma",
  # de_type = "voom",
  de_method = "edgeR",
  de_type = "LRT",
  n_threads = 16
)

run_de
data("hagai_toy")
hagai_toy
head(sce$TN)

wulu=run_de(
  sce,
  meta = NULL,
  replicate_col = "orig.ident",
  cell_type_col = "cell4",
  label_col = "TN",
  min_cells = 3,
  min_reps = 2,
  min_features = 0,
  de_family = "singlecell",
  de_method = "wilcox",
  # de_type = "LRT",
  # n_threads = 16
)

wulu=run_de(
  sce,
  meta = NULL,
  replicate_col = "orig.ident",
  cell_type_col = "cell4",
  label_col = "TN",
  min_cells = 3,
  min_reps = 2,
  min_features = 0,
  de_family = "pseudobulk",
  de_method = "edgeR",
  de_type = "LRT",
  n_threads = 16
)
writexl::write_xlsx(wulu,"癌vs非肿瘤胰腺正的在癌升高负的在正常升高 edgeR pseudobulk.xlsx")

wulu=run_de(
  sce,
  meta = NULL,
  replicate_col = "orig.ident",
  cell_type_col = "cell4",
  label_col = "TN",
  min_cells = 3,
  min_reps = 2,
  min_features = 0,
  de_family = "pseudobulk",
  de_method = "limma",
  de_type = "voom",
  n_threads = 16
)
writexl::write_xlsx(wulu,"癌vs非肿瘤胰腺正的在癌升高负的在正常升高 limma pseudobulk.xlsx")
# wuqiwu=to_pseudobulk(
#   sce,
#   meta = NULL,
#   replicate_col = "orig.ident",
#   cell_type_col = "celltype",
#   label_col = "TN",
#   min_cells = 3,
#   min_reps = 2,
#   min_features = 0
# 
# )

scefibroepi_T=subset(sce_T,subset = cell4%in%c('Fibroblast',"Epithelial"))

scefibroepi_T$cell4=factor(scefibroepi_T$cell4,levels = c('Fibroblast',"Epithelial"))
wulu=run_de(
  scefibroepi_T,
  meta = NULL,
  replicate_col = "orig.ident",
  cell_type_col = "TN",
  label_col = "cell4",
  min_cells = 3,
  min_reps = 2,
  min_features = 0,
  de_family = "pseudobulk",
  # de_method = "limma",
  # de_type = "voom",
  de_method = "edgeR",
  de_type = "LRT",
  n_threads = 16
)
writexl::write_xlsx(wulu,"CAF vs EPI 正的在间质升高负的在肿瘤升高 limma pseudobulk.xlsx")
writexl::write_xlsx(wulu,"CAF vs EPI 正的在间质升高负的在肿瘤升高 edgeR pseudobulk.xlsx")
countda=data.table::fread(r"H:\download.big.ac.cn\gsa\CRA001160\count-matrix.txt")



CD3D=list()
fs=list.files(r"(H:\download.big.ac.cn\gsa\CRA001160\other_files\results\)",recursive = T,pattern = "geneManifest.txt",full.names = T)
for(x in 1:35){
CD3D[[x]]=subset(data.table::fread(fs[[x]],header = T),Symbol=="CD3D")
}
CD3D=data.table::rbindlist(CD3D)
CD3D$fs=fs

CD79A=list()
fs=list.files(r"(H:\download.big.ac.cn\gsa\CRA001160\other_files\results\)",recursive = T,pattern = "geneManifest.txt",full.names = T)
for(x in 1:35){
  CD79A[[x]]=subset(data.table::fread(fs[[x]],header = T),Symbol=="CD79A")
}
CD79A=data.table::rbindlist(CD79A)
CD79A$fs=fs
  
#其实是有一个样本导致了被过滤


# genes_to_check = c('Ins1', 'Gcg', 'Ppy', 'Sst', 'Chga', 'Krt19', 'Amy2b', 'Pecam1', 'Pdgfra', 'Ptprc', 'Ghrl')
# genes_to_check = c('Btg2','Sgk1','Jun','Srsf2')

为什么要去除核糖体和线粒体
如果线粒体RNA过高，也同样预示着细胞有破损。因为当细胞破损时，细胞质RNA会跑出来，但是线粒体RNA由于有线粒体膜的包裹，不会溢出。因此，当细胞膜有破损时，线粒体RNA所占比例会很高。注意：当细胞出现apoptosis, necrosis的时候，也会有这种现象。

核糖体RNA占比较高时，可能是因为细胞内出现了较多的RNA降解。在全长单细胞转录组中，3’偏好性可用于检测细胞内是否存在大量RNA降解

#########FTO

library(scRNAtoolVis)
markers=FindMarkers(sce,ident.1 = 3)
# plot 这个非常适合有比例的，而直接的稚嫩用单纯的了
markerVocalno(markers = markers,
              topn = 5,
              labelCol = ggsci::pal_aaas()(9))
FindMarkers(ident.1 = "FTO", group.by = 'FTO_group',subset.ident =c("3","4"))
sce@meta.data$cluster_column
View(sce@meta.data)


Here is a solution using dplyr and ggplot2:
  
library(Seurat) 
library(dplyr)
library(ggplot2)

meta.data <- pbmc_small[[]]

# create random classifications for the sake of this example
meta.data$condition <- sample(c('A', 'B', 'C'), nrow(meta.data), replace = TRUE)

counts <- group_by(meta.data, condition, res.1) %>% summarise(count = n())

ggplot(counts, aes(res.1, count, fill = condition)) +
  geom_bar(stat = 'identity')


# library(Chord)
# chord(seu=sce,doubletrat=0.01,overkill=T,outname="the name you want")
# chord(seu=sce,doubletrat=NULL,overkill=T,outname="the name you want")
# load(r"(C:\Users\zxh\Desktop\R\scCancer\nor\seu.robj)")
# sce=seu
# dbl=read.csv(r"(C:\Users\zxh\Desktop\R\scCancer\nor\the name you want_doublet.csv)")
# seu@meta.data$x=rownames(seu@meta.data)
# seu@meta.data$x=ifelse(seu@meta.data$x%in%dbl$x,1,0)
# seu2=subset(seu,x==0)
if(F){
sce$TN=substr(sce$orig.ident,1,1)
sce@meta.data$TN=factor(sce@meta.data$TN,levels=c("T","N"))
sce_sub=subset(sce,seurat_clusters%in%c(4,8,13,14,18),invert = TRUE)
sce_sub$IM=ifelse(sce_sub$seurat_clusters%in%c(2,3,6,7,9,10,12,17),"myoCAF","iCAF")
SCpubr::do_DimPlot(raster=F,sce_sub,group.by = "TN")
Idents(sce_sub)=sce_sub$TN
}
if(F){
sce$TU=substr(sce$orig.ident,1,1)
sce@meta.data$TU=factor(sce@meta.data$TU,levels=c("T","U"))
sce_sub=subset(sce,seurat_clusters%in%c(8,9,10,11,12,13),invert = TRUE)
sce_sub$IM=ifelse(sce_sub$seurat_clusters%in%c(0,1,2,3,7),"iCAF","myoCAF")
SCpubr::do_DimPlot(raster=F,sce_sub,label=T)
SCpubr::do_DimPlot(raster=F,sce_sub,group.by = "TU",label=T)
SCpubr::do_DimPlot(raster=F,sce_sub,group.by = "IM")
SCpubr::do_DimPlot(raster=F,sce_sub,reduction = "umap",label=T,pt.size = 1)
VlnPlot(sce_sub,group.by = "IM", features = c("FTO"),log = F,pt.size=0.1)+theme(legend.position = "none")
print(p6)
Idents(sce_sub)=sce_sub$TU
}



diff.expr.genesT <- FindMarkers(
  object = sce_sub,
  group.by = "IM",
  ident.1 = "iCAF",
  only.pos = F,
  min.pct = 0,
  logfc.threshold = 0.5,
  verbose = T,subset.ident = "T"
)
diff.expr.genesT=subset(diff.expr.genesT,p_val_adj<0.05)
diff.expr.genesT$genename=rownames(diff.expr.genesT)

diff.expr.genesU <- FindMarkers(
  object = sce_sub,
  group.by = "IM",
  ident.1 = "iCAF",
  only.pos = F,
  min.pct = 0,
  logfc.threshold = 0.5,
  verbose = T,subset.ident = "U"
)
diff.expr.genesU=subset(diff.expr.genesU,p_val_adj<0.05)
diff.expr.genesU$genename=rownames(diff.expr.genesU)

diff.expr.genesN <- FindMarkers(
  object = sce_sub,
  group.by = "IM",
  ident.1 = "iCAF",
  only.pos = F,
  min.pct = 0,
  logfc.threshold = 0.5,
  verbose = T,subset.ident = "N"
)
diff.expr.genesN=subset(diff.expr.genesN,p_val_adj<0.05)
diff.expr.genesN$genename=rownames(diff.expr.genesN)

if(F){
Idents(sce_sub)=sce_sub$IM
diff.expr.genesiCAF <- FindMarkers(
  object = sce_sub,
  group.by = "TN",
  ident.1 = "T",
  only.pos = F,
  min.pct = 0,
  logfc.threshold = 0.5,
  verbose = T,subset.ident = "iCAF"
)
diff.expr.genesiCAF=subset(diff.expr.genesiCAF,p_val_adj<0.05)
diff.expr.genesiCAF$genename=rownames(diff.expr.genesiCAF)

diff.expr.genesmyoCAF <- FindMarkers(
  object = sce_sub,
  group.by = "TN",
  ident.1 = "T",
  only.pos = F,
  min.pct = 0,
  logfc.threshold = 0.5,
  verbose = T,subset.ident = "myoCAF"
)
diff.expr.genesmyoCAF=subset(diff.expr.genesmyoCAF,p_val_adj<0.05)
diff.expr.genesmyoCAF$genename=rownames(diff.expr.genesmyoCAF)

sce_sub$TN
diff.expr.genes_libra_TNdiff=Libraavoiddup::run_de(
  sce_sub,
  meta = NULL,
  replicate_col = "orig.ident",
  cell_type_col = "IM",#"TN",#"celltype",
  label_col = "TN",
  min_cells = 3,
  min_reps = 10,#这个一般是最小亚分组再除以2
  min_features = 1,#0,#1*length(unique(sce_sub$orig.ident)),#我觉得是大于0就可以了
  # de_family = "singlecell",
  # de_method = "wilcox",#"limma",
  de_family = "pseudobulk",#
  # de_method = "limma",#"",
  # de_type = "voom",
  de_method = "edgeR",
  de_type = "LRT",
  n_threads = 16
)
diff.expr.genes_libra_TNdiff=subset(diff.expr.genes_libra_TNdiff,p_val_adj<0.05)
diff.expr.genes_libra_TNdiff=subset(diff.expr.genes_libra_TNdiff,abs(diff.expr.genes_libra_TNdiff$avg_logFC)>0.5)

}

unique(sce_sub$orig.ident)

diff.expr.genes_libra=Libraavoiddup::run_de(
  sce_sub,
  meta = NULL,
  replicate_col = "orig.ident",
  cell_type_col = "TU",#"TN",#"celltype",
  label_col = "IM",
  min_cells = 3,
  min_reps = 10,#这个一般是最小亚分组再除以2
  min_features = 1,#0,#1*length(unique(sce_sub$orig.ident)),#我觉得是大于0就可以了
  # de_family = "singlecell",
  # de_method = "wilcox",#"limma",
  de_family = "pseudobulk",#
  # de_method = "limma",#"",
  # de_type = "voom",
  de_method = "edgeR",
  de_type = "LRT",
  n_threads = 16
)
diff.expr.genes_libra=subset(diff.expr.genes_libra,p_val_adj<0.05)
diff.expr.genes_libra=subset(diff.expr.genes_libra,abs(diff.expr.genes_libra$avg_logFC)>0.5)


save(diff.expr.genesT,diff.expr.genesU,diff.expr.genes_libra,file="cyto myo imm caf GSE202051.Rdata")
load(file="cyto myo imm caf GSE202051.Rdata")
save(diff.expr.genesT,diff.expr.genesN,diff.expr.genes_libra,file="cyto myo imm caf CRA001160.Rdata")
load(file="cyto myo imm caf CRA001160.Rdata")
# a1=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\scCancer\allsurface hpa not correct include secret.xlsx)")
# a2=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\scCancer\MAKER GENE LIST(1).xlsx)",col_names = F)
a3=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\scCancer\HGNC all cd group-471.xlsx)")
aa=a3$`Approved symbol`
de_genes=de_genes2
# ss1=de_genes[which(de_genes$pct.1>de_genes$pct.2),]
# ss2=de_genes[which(de_genes$pct.2>de_genes$pct.1),]
ss1=de_genes[which(de_genes$pct.1>0.5&de_genes$pct.2<0.2),]
ss2=de_genes[which(de_genes$pct.2>0.5&de_genes$pct.1<0.2),]

de_genes=diff.expr.genesU#diff.expr.genesT
de_genes=subset(de_genes,abs(de_genes$avg_log2FC)>1)
de_genes$gene=rownames(de_genes)
de_genes=de_genes[rownames(de_genes)%in%aa,]
de_genes1=de_genes

de_genes=diff.expr.genesT#
de_genes=subset(de_genes,abs(de_genes$avg_log2FC)>1)
de_genes$gene=rownames(de_genes)
de_genes=de_genes[rownames(de_genes)%in%aa,]
de_genes2=de_genes


# de_genes1=de_genes1[intersect(de_genes1$gene,de_genes2$gene),]
# de_genes2=de_genes2[intersect(de_genes1$gene,de_genes2$gene),]

de_genes=diff.expr.genesN#
de_genes=subset(de_genes,abs(de_genes$avg_log2FC)>1)
de_genes$gene=rownames(de_genes)
de_genes=de_genes[rownames(de_genes)%in%aa,]
de_genes3=de_genes

de_genes2=de_genes2[intersect(de_genes3$gene,de_genes2$gene),]
de_genes3=de_genes3[intersect(de_genes3$gene,de_genes2$gene),]

de_genes=diff.expr.genes_libra#diff.expr.genesT
de_genes=subset(de_genes,abs(de_genes$avg_logFC)>1)
de_genes=de_genes[de_genes$gene%in%aa,]
de_genes$cell_type=NULL
de_genes=de_genes[duplicated(de_genes$gene),]
de_genes4=de_genes

Idents(sce_sub)=sce_sub$TN
caf_ConservedMarkers=FindConservedMarkers(sce_sub, ident.1 = "T",ident.2 = "N", grouping.var = "IM")
caf_ConservedMarkers=caf_ConservedMarkers[caf_ConservedMarkers$iCAF_avg_log2FC>0,]
caf_ConservedMarkers$gene=rownames(caf_ConservedMarkers)
caf_ConservedMarkers=caf_ConservedMarkers[caf_ConservedMarkers$gene%in%de_genes_fibro$gene,]
caf_ConservedMarkers=caf_ConservedMarkers[caf_ConservedMarkers$iCAF_pct.1>0.8,]
caf_ConservedMarkers=caf_ConservedMarkers[caf_ConservedMarkers$myoCAF_pct.1>0.8,]
caf_ConservedMarkers=caf_ConservedMarkers[caf_ConservedMarkers$iCAF_avg_log2FC>1,]
caf_ConservedMarkers=caf_ConservedMarkers[caf_ConservedMarkers$myoCAF_avg_log2FC>1,]
caf_ConservedMarkers$diff=caf_ConservedMarkers$iCAF_avg_log2FC-caf_ConservedMarkers$myoCAF_avg_log2FC
caf_ConservedMarkers=subset(caf_ConservedMarkers,abs(caf_ConservedMarkers$diff)<0.5)
#COL1A2 最保守
#COL3A1 最保守
write.csv(caf_ConservedMarkers,file = "caf_ConservedMarkers.csv")
# de_genes4=de_genes4[intersect(de_genes3$gene,de_genes4$gene),]
# de_genes5=de_genes5[intersect(de_genes5$gene,de_genes4$gene),]

de_genes=diff.expr.genesiCAF#
de_genes=subset(de_genes,abs(de_genes$avg_log2FC)>1)
de_genes$gene=rownames(de_genes)
de_genes=de_genes[rownames(de_genes)%in%aa,]
de_genes222=de_genes

de_genes=diff.expr.genesmyoCAF#
de_genes=subset(de_genes,abs(de_genes$avg_log2FC)>1)
de_genes$gene=rownames(de_genes)
de_genes=de_genes[rownames(de_genes)%in%aa,]
de_genes333=de_genes

de_genes_fibro=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\results\ori\comb\DE-Genes.RDS)")
de_genes_fibro=de_genes_fibro[de_genes_fibro$cluster%in%c(3,4),]
# THY1 才是真正 tumor caf
de_genes=diff.expr.genes_libra_TNdiff#diff.expr.genesT
de_genes=subset(de_genes,abs(de_genes$avg_logFC)>1)
# de_genes=de_genes[de_genes$gene%in%aa,]
de_genes$cell_type=NULL
de_genes=de_genes[duplicated(de_genes$gene),]
de_genes444=de_genes



#找不到normal marker high
cafvsepi=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\scCancer\单细胞分析结果发刘明洋老师\CAF vs EPI 正的在间质升高负的在肿瘤升高 edgeR pseudobulk.xlsx)")
cafvsepi=subset(cafvsepi,cafvsepi$p_val_adj<0.05)
cafvsepi=subset(cafvsepi,cafvsepi$avg_logFC>0)
de_genes444=de_genes444[de_genes444$gene %in%de_genes_fibro$gene,]
de_genes444=subset(de_genes444,de_genes444$avg_logFC<0)

de_genes333=de_genes333[intersect(de_genes333$gene,de_genes222$gene),]
de_genes222=de_genes222[intersect(de_genes333$gene,de_genes222$gene),]
de_genes222[de_genes222$avg_log2FC<0,]
de_genes333=de_genes333[de_genes333$gene%in%de_genes_fibro$gene,]
de_genes222=de_genes222[de_genes222$gene%in%de_genes_fibro$gene,]

de_genes_fibro$gene
aa=unique(c(a1$Gene,a2$...2))


table(rownames(diff.expr.genesT)%in%aa)
table(rownames(diff.expr.genesN)%in%aa)
table(diff.expr.genes_libra$gene%in%aa)
# de_genes=diff.expr.genes
# openxlsx::write.xlsx(de_genes,file = "caf findmark.xlsx")

de_genes=diff.expr.genesT[rownames(diff.expr.genesT)%in%aa,]

diff.expr.genesT["FAP",]
ss1=de_genes[which(de_genes$pct.1>0.5&de_genes$pct.2<0.2),]
ss2=de_genes[which(de_genes$pct.2>0.5&de_genes$pct.1<0.2),]
ss=rbind(ss1,ss2)
ss

de_genes=diff.expr.genesN[rownames(diff.expr.genesN)%in%aa,]
de_genes["PDGFRA",]
ss1=de_genes[which(de_genes$pct.1>0.5&de_genes$pct.2<0.2),]
ss2=de_genes[which(de_genes$pct.2>0.5&de_genes$pct.1<0.2),]
ss=rbind(ss1,ss2)
ss
# pbmc$groups <- rep(c('stim','control'),each = 1319)
# add celltype
# pbmc$celltype <- Seurat::Idents(sce)
key_genes=c("FTO","ACTA2","FAP")
de_genes2=de_genes[!duplicated(de_genes$gene),]%>%
  group_by(cluster) %>%
  # mutate(var = var(avg)) %>%
  # ungroup() %>%
  top_n(5,avg_log2FC)# %>%
p <- list()
p[[12]]=scRNAtoolVis::jjDotPlot(object = sce,
                                gene = unique(c(key_genes,de_genes2$gene)),
                                #unique(c(key_genes,de_genes[which(de_genes$gene%in%genes_to_check),]$gene)),#top3pbmc.markers$gene,
                                id = 'seurat_clusters',
                                # anno = F,
                                # markerGene = top3pbmc.markers,#finallmarkers结果直接load进行,而且只能选一个，注意互斥
                                xtree = F,
                                ytree = F,
                                rescale = T,#重整为0-1
                                rescale.min = 0,#rescale.min = -2,
                                rescale.max = 1,#rescale.max = 2,
                                midpoint =0.5,#midpoint =0,
                                dot.col	= c("white","#990000"),#c('blue','white','red')
                                point.geom = T,
                                # point.shape = 22,
                                tile.geom = F,
                                tree.pos = 'left',
                                same.pos.label = T,
                                # yPosition = 10.3,
                                # split.by = 'groups',
                                # split.by.aesGroup = T,#T颜色不用于分组Set split.by.aesGroup = T to turn off the colors grouped by groups:
                                gene.order = unique(c(key_genes,de_genes2$gene)),
                                #unique(c(key_genes,de_genes[which(de_genes$gene%in%genes_to_check),]$gene)),#rev(top3pbmc.markers$gene),
                                cluster.order = rev(0:length(unique(sce$seurat_clusters))-1),
                                plot.margin = c(1,1,1,1))# 上下左右距离
p[[12]]
pdf("qqq.pdf",height = 8,width = 30)
p[[12]]
dev.off()
de_genes["ENG",]
sce_sub$Class="CAF"

sce_sub$TN=substr(sce_sub$orig.ident,1,1)

aaa=Libraavoiddup:::check_inputs(sce_sub ,# meta = sce@meta.data,
                                 replicate_col = "orig.ident",
                                 cell_type_col = "TU",
                                 label_col = "IM")
wula_sum=Libraavoiddup::to_pseudobulk(aaa$expr,aaa$meta,
                       replicate_col = "orig.ident",
                       cell_type_col = "TU",
                       label_col = "IM",
                       min_cells = 3,
                       min_reps = 5,
                       min_features = 1,
                       external = F)
wula_sum=as.data.frame(wula_sum[["T"]])
"PTGDS"%in%rownames(wula_sum)
wula_sum["PTGDS",]

?to_pseudobulk

# _T$RNA_snn_res.0.3

diff.expr.genesU <- FindMarkers(
  object = sce_sub,
  group.by = "TN",
  ident.1 = "T",
  only.pos = T,
  min.pct = 0,
  min.cells.feature = 0,
  min.cells.group = 0,

  logfc.threshold = 2,
  verbose = T
)


library(SCpubr)

sample <- sce
assay <- "RNA"

# Retrieve prior knowledge network.
network <- decoupleR::get_dorothea(organism = "human",
                                   levels = c("A", "B", "C"))
decoupleR#bulk scRNA
# Run weighted means algorithm.
activities <- decoupleR::run_wmean(mat = as.matrix(sample@assays[[assay]]@data),
                                   network = network,
                                   .source = "source",
                                   .targe = "target",
                                   .mor = "mor",
                                   times = 100,
                                   minsize = 5)
acts <- decoupleR::run_wmean(mat=mat, net=net, .source='source', .target='target',
                  .mor='weight', times = 100, minsize = 5)
save(activities,file="activities.Rdata")
out <- SCpubr::do_TFActivityPlot(sample = sample,
                                 activities = activities)
p <- out$heatmaps$average_scores
p

out <- SCpubr::do_TFActivityPlot(sample = sample,
                                 activities = activities,
                                 split.by = "phase",
                                 min.cutoff = 0.1,
                                 max.cutoff = 0.7,
                                 plot_FeaturePlots = TRUE,
                                 plot_GeyserPlots = TRUE)
p <- out$heatmaps$average_scores
p

##
order.genes <- order.genes[!grepl("ENSMPUG", order.genes)]
pdf("Fig5e.MP_M1_pseudotime_heatmap.pdf", 5/2.54*2,9/2.54*1.5)
plot_pseudotime_heatmap(cds[order.genes,], num_clusters = 4, cores = 1, show_rownames = F, return_heatmap = T)
dev.off()



在以前的文章中，我们已经介绍了两种细胞注释方法：软件注释和利用 marker 基因注释。即便使用 CellAssign 等单细胞注释软件做注释，也需要提供细胞类型及其 marker 基因信息，而且无论用哪种软件进行注释，最终都要通过 marker 基因的表达情况来进行检验。因此，特定细胞类型的 marker 基因表达数据库就显得格外重要。 



现在有不少已发表的 marker 基因数据库，使用比较广泛地有 CellMarker、PanglaoDB 等，但科研人员在分析的过程中使用这些经典数据库也常发现注释结果不好甚至出现错误的情况。博奥晶典针对这个难点，提出了自己的解决方案，通过已发表的 2000+ 篇高质量单细胞文献数据进行总结，构建了多个物种和组织的细胞类群 marker 数据库，称为 SingleCellBase。 