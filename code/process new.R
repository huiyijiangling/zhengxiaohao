CRA001160_meta=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\scCancer\CRA001160\CRA001160 meta.xls)")#发表
library(Seurat)
library(clustree)
sce=readRDS(r"(C:\Users\zxh\Desktop\R\scCancer\results\snuc_NCPPAAD\expr.RDS)")#发表

# sNUC 我直接合并吧
Tosti_meta=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\scCancer\Tosti_data\source data\oriarticle\Tosti meta.xlsx)")#发表
GSE202051_meta=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\sourcedata\GSE202051 meta.xlsx)")
sNUC_meta=data.table::rbindlist(list(Tosti_meta,GSE202051_meta),use.names = T,fill = T)
table(sce@meta.data$sample.ident)
sNUC_meta$Sample[!sNUC_meta$Sample%in%sce@meta.data$sample.ident]
meta_new_rownames=rownames(sce@meta.data)

meta_new=dplyr::left_join(
  sce@meta.data,
  sNUC_meta,
  by = c("sample.ident"="Sample"),
  
)
rownames(meta_new)=meta_new_rownames
sce@meta.data=meta_new  
head(sce@meta.data)
sce@meta.data$dataset = ifelse(sce@meta.data$sample.ident%in%list.files(r"(C:\Users\zxh\Desktop\R\scCancer\Tosti_data\adultpancreas2020\results\)",pattern = "_"),"Tosti_Adult_Pancreas",
                               ifelse(sce@meta.data$sample.ident%in%  list.files(r"(C:\Users\zxh\Desktop\R\scCancer\Tosti_data\chronicpancreatitis\results\)",pattern = "_"),"Tosti_Chronic_Pancreatitis",
                                      ifelse(sce@meta.data$sample.ident%in% list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\results\)",pattern = "T|U"),"Hwang_pancreatic_cancer","ERROR")))
table(sce@meta.data$dataset)
colnames(sce)==rownames(sce@meta.data)
sce@meta.data$cellname=rownames(sce@meta.data)
table(sce@meta.data$OS.y)
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
# subset 不同的间隔室
# sc


sce.all=sce
for (res in c(0.3)) {#0.1, 0.2, 0.3,0.5,0.6,0.7,0.8,1
  sce.all=FindClusters(sce.all,  
                       resolution = res,
                       algorithm = 1)}
apply(sce.all@meta.data[,grep("RNA_snn_res",#RNA_snn_res.
                              colnames(sce.all@meta.data))],2,table)
sce=sce.all
sce.all=NULL
Idents(sce)=sce$RNA_snn_res.0.3

colnames(sce)==rownames(sce@meta.data)
rownames(sce)
# # GSE202051_phe=subset(GSE202051_phe,TU=="U")
# wula0=AggregateExpression(object = sce,group.by = c("sample.ident"),slot = "counts")#counts
# wula0=AverageExpression(object = sce,group.by = c("sample.ident"),slot = "counts")#counts
# # wula=AverageExpression(object = sce_T,group.by = c("seurat_clusters"),slot = "counts")#counts
# # wula=AverageExpression(object = sce_T,group.by = c("orig.ident","seurat_clusters"),slot = "counts")#counts
# wula=as.data.frame(wula0)
# colnames(wula)=stringr::str_split(colnames(wula),"\\.",simplify = T)[,2]
# wula=wula[,colnames(wula)%in%GSE202051_phe$ID_standard]
# GSE202051_phe=subset(GSE202051_phe,GSE202051_phe$ID_standard%in%colnames(wula))
# 
# wula=subset(wula,select=match(GSE202051_phe$ID_standard,colnames(wula)))
# 
# # wula=t(wula)
# # wula=as.data.frame(wula)
# sce_avg=AverageExpression(sce)
# wula=as.data.frame(sce_avg)
# #Plot()
# scCancer::runSurvival(
#   features=c("ITGA3"),
#   data=wula,
#   surv.time=GSE202051_phe$OS.time,
#   surv.event=GSE202051_phe$OS,
#   cut.off = 0.5,
#   savePath = NULL
# )
# plot_pseudotime_heatmap


if(T){
#to
exprSet3$sample=rownames(exprSet3)
dat=merge(phe,exprSet3,by="sample")
# dat$gender=factor(dat$gender)
# dat$stage=factor(dat$stage)
library(survival)
library(survminer)
dat$event=as.numeric(dat$event)
dat$time=as.numeric(dat$time)
s=Surv(time, event) ~   KIF15+FEN1+ZFP69B+SP6+SPARC+TTF2+MSI2+KYNU+ACLY+KIF21B+SLC12A7+ZNF823

#MSI2+KYNU+METTL2B+ACLY+CXCR4+KIF21B+VSNL1+SERPINE1+TNFAIP2
#KIF15+FEN1+TTF2+MSI2+KYNU+ZNF562+ACLY+CXCR4+KIF21B+SLC12A7+ZNF823+NHS  

model <- coxph(s, data =dat)
summary(model,data=dat)
options(scipen=1)

ggforest(model, data =dat, 
         main = "Forest plot (Hazard ratio)", 
         cpositions = c(0.10, 0.22, 0.4), 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 4)+
  ggsave('Forest.pdf',width = 8,height=6,dpi = 600)

}

# sce@meta.data$state = ifelse(sce@meta.data$sample.ident%in%list.files(r"(C:\Users\zxh\Desktop\R\scCancer\Tosti_data\adultpancreas2020\results\)",pattern = "_"),"Normal P",
#                                ifelse(sce@meta.data$sample.ident%in%  list.files(r"(C:\Users\zxh\Desktop\R\scCancer\Tosti_data\chronicpancreatitis\results\)",pattern = "_"),"CP",
#                                       ifelse(sce@meta.data$sample.ident%in% list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\results\)",pattern = "U"),"Untreat_PDAC",
#                                              ifelse(sce@meta.data$sample.ident%in% list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\results\)",pattern = "T"),"Treated","ERROR"))))
# #







match(Tosti_meta$sample,sce@meta.data$dataset)
table(sce@meta.data$orig.ident)
table(sce@meta.data$sample.ident)
library(SCpubr)
source(r"(C:\Users\zxh\Desktop\R\rsetup\theme_big_simple1.R)")
pdf("dimplot umap all.pdf",height = 4,width = 4)
p3=SCpubr::do_DimPlot(sce,reduction = "umap",label=T,raster=T,legend.position="none")#legend.position = "right",
print(p3)
dev.off()

pdf("dimplot umap dataset.pdf",height = 8,width = 8)
p4=SCpubr::do_DimPlot(sce,group.by = "dataset",reduction = "umap",label=T,legend.position = "top",raster=T,plot.axes = T,pt.size = 0.1)
print(p4)
dev.off()

pdf("dimplot umap highlight.pdf",height = 8,width = 8)
Idents(sce)=sce$dataset
p4=SCpubr::do_DimPlot(sce,reduction = "umap",label=T,legend.position = "top",raster=T,plot.axes = T,pt.size = 0.1,idents.highlight = "Tosti_Adult_Pancreas")#,group.by = "dataset"
Idents(sce)=sce$seurat_clusters
print(p4)
dev.off()
table(sce$dataset)

pdf("dimplot umap highlight cp.pdf",height = 8,width = 8)
Idents(sce)=sce$dataset
p4=SCpubr::do_DimPlot(sce,reduction = "umap",label=T,legend.position = "top",raster=T,plot.axes = T,pt.size = 0.1,idents.highlight = "Tosti_Adult_Pancreas")#,group.by = "dataset"
Idents(sce)=sce$seurat_clusters
print(p4)
dev.off()
table(sce$dataset)

pdf("VlnPlot sce FTO .pdf",height = 8,width = 8)
p6=VlnPlot(sce,group.by = "seurat_clusters", features = c("FTO"),log = F,pt.size=0)+theme(legend.position = "none")
print(p6)
dev.off()

pdf("VlnPlot sce FTO .pdf",height = 8,width = 8)
p6=VlnPlot(sce,group.by = "seurat_clusters", features = c("FTO"),log = F,pt.size=0)+theme(legend.position = "none")
print(p6)
dev.off()

pdf("FeaturePlot KRT19.pdf",height = 8,width = 8)
p5=SCpubr::do_FeaturePlot(sce, features = "KRT19",enforce_symmetry = F)#为了好看
print(p5)
dev.off()


pdf("FeaturePlot KRT19.pdf",height = 8,width = 8)
p5=SCpubr::do_FeaturePlot(sce, features = "KRT19",enforce_symmetry = F)#为了好看
print(p5)
dev.off()




  sce_T@meta.data$FTO_group=wula[match(sce_T@meta.data$orig.ident,wula$sample),"FTO_group"]
  sce_T@meta.data$orig.ident=factor(sce_T@meta.data$orig.ident)
  
  sce@meta.data$assignment=wula[match(sce@meta.data$orig.ident,wula$sample),"FTO_group"]
  table(sce@meta.data$assignment,useNA="ifany")
  
#   
#   
# C:\Users\zxh\Desktop\R\scCancer\Tosti_data\adultpancreas2020
# single.savePaths <- c(
#   list.files(r"(C:\Users\zxh\Desktop\R\scCancer\Tosti_data\adultpancreas2020\results\)",full.names = T,pattern = "_"),
#   list.files(r"(C:\Users\zxh\Desktop\R\scCancer\Tosti_data\chronicpancreatitis\results\)",full.names = T,pattern = "_"),
#   list.files(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\results\)",full.names = T,pattern = "T|U")
# )

  library(ggplot2)
  library(aplot)
  # ori
  aa=readr::read_tsv(r"(C:\Users\zxh\Downloads\table (6).tsv)")
  bb=readr::read_tsv(r"(C:\Users\zxh\Downloads\table (7).tsv)")
  aa=aa[aa$`q-Value`<0.05,]
  bb=bb[bb$`q-Value`<0.05,]
  cc=merge(aa,bb,by="Correlated Gene")
  cc=subset(cc,sign(cc$`Spearman's Correlation.x`)==sign(cc$`Spearman's Correlation.y`))
  write.csv(cc,file="common_fto_rna_protein.csv")
  gene_cluster$
  gene_cluster=read.csv(r"(common_fto_rna_protein.csv)")
  dot_plot <- gene_cluster %>%
    # mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>%
    # filter(count > 0, `% Expressing` > 1) %>%
    ggplot(aes(x=Spearman.s.Correlation.x, y = Correlated.Gene, color = q.Value.x, size = q.Value.x)) +
    geom_point() +
    cowplot::theme_cowplot() +
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab(NULL) +
    theme(axis.ticks = element_blank()) +
    scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)') +
    scale_y_discrete(position = "right")
  dot_plot
  # dot_plot %>%
  #   insert_left(ggtree_plot, width=.2) %>%
  #   insert_top(labels, height=.02) %>%
  #   insert_top(ggtree_plot_col, height=.1)

  
  
  
  
  
  
  pdf("VlnPlot sce T FTO2.pdf",height = 8,width = 16)
  library(ggpubr)
  p7=VlnPlot(sce,group.by = "seurat_clusters", features = c("FTO"),log = F,pt.size=0,split.by = "dataset",split.plot=T,cols =color[1:3] )+theme(legend.position = "top")+ stat_compare_means( method = "wilcox.test", paired = F, label = "p.signif")+    stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar",
    width = 0.6,
    position = position_dodge(width = .70)
  ) # + stat_compare_means(size = your_font_size, label = "p.signif")#ggpubr::stat_compare_means(label = "p.signif")
  print(p7)
  dev.off()
  #