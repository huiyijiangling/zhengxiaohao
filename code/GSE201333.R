
#pdx 看不到
#no reduction

library(reticulate)
library(Seurat)
reticulate::use_condaenv("C:/Users/zxh/AppData/Local/r-miniconda/envs/py37")
# py_install("scanpy")
# sc <- import("scanpy")
# 
# atlas.data <- sc$read_h5ad("~/kidney_cell_atlas_fetal/Fetal_full_v3.h5ad")
# 
# counts <- t(atlas.data$layers["counts"])
# colnames(counts) <-  atlas.data$obs_names$to_list()
# rownames(counts) <-  atlas.data$var_names$to_list()
# counts <- Matrix::Matrix(as.matrix(counts), sparse = T)
# 
# seurat <- CreateSeuratObject(counts)
# seurat <- AddMetaData(seurat,  atlas.data$obs)
# 直接使用scanpy 来解决问题

load(file= "matrix.rdata")
scRNA@meta.data$free_annotation=stringr::str_to_upper(scRNA@meta.data$free_annotation)
scRNA@meta.data$cell_ontology_class=stringr::str_to_upper(scRNA@meta.data$cell_ontology_class)
library(Seurat)
# scRNA$csfgroup=ifelse(scRNA$group=="HC","HC","PD-DLB")
# 名字非常不标准 但是可以倾向于 cell_ontology_class 而 free_annotation 更细有点不标准
if(F){
table(scRNA@meta.data$donor,useNA="ifany")
# TSP1   TSP2   TSP3   TSP4   TSP5   TSP6   TSP7   TSP8   TSP9  TSP10  TSP11  TSP12  TSP13  TSP14  TSP15 
# 39431 115557   2447  20930   2592   7978  57127  15067   6699  23463   2488  11505    538 171719   5611 
table(scRNA@meta.data$organ_tissue,useNA="ifany")
# Bladder           Blood     Bone_Marrow             Eye             Fat           Heart          Kidney 
# 24583           50115           12297           10650           20263           11505            9641 
# Large_Intestine           Liver            Lung      Lymph_Node         Mammary          Muscle        Pancreas 
# 13680            5007           35682           53275           11375           30746           13497 
# Prostate  Salivary_Gland            Skin Small_Intestine          Spleen          Thymus          Tongue 
# 16375           27199            9424           12467           34004           33664           15020 
# Trachea          Uterus     Vasculature 
# 9522            7124           16037 
table(scRNA@meta.data$method,useNA="ifany")
# 10X smartseq2 
# 456101     27051 

table(scRNA@meta.data$anatomical_information)
# Abdomen           Anterior              Aorta      AortaVeneCava              Atria              Chest 
# 7732               4805               7801               2522               8090               1971 
# Conjunctiva         Cornea-etc   CoronaryArteries          Diaphragm             Distal          Endocrine 
# 2084               3505               4867               2252              13769               1520 
# Endometrium           Exocrine           Inguinal                MAT      MedialDistal          Myometrium 
# 3813               5278               8040              10371                307               3025 
# Neuron    Neuroretina-etc            Orbital            Parotid          Posterior           Proximal 
# 292                675                  2              20171               8824              11780 
# SCAT         Sclera-etc      Submandibular  Supradiaphagmatic          Ventricle            abdomen 
# 9892               1431               6196               4181               3138                417 
# anterior              aorta              atria              chest          diaphragm             distal 
# 562                583                189                427               5573               1195 
# exocrine           inguinal      lacrimalgland                nan           noCornea            parotid 
# 6699                643                248             270984                 36                832 
# posterior           proximal   proxmedialdistal    rectusabdominus supradiaphragmatic          ventricle 
# 643               2126              17368              11898                519                 88 
# vertebralbody              whole 
# 3655                133 
scRNA@meta.data$free_annotation
# 320 Levels: AREG High T-Cell AREG Low T-Cell Adventitial Fibroblast ... vein_endothelial_cell
scRNA@meta.data$cell_ontology_class
# 177 Levels: acinar cell of salivary gland adipocyte adventitial cell alveolar fibroblast artery endothelial cell ... vein endothelial cell
scRNA@meta.data$compartment
# Levels: endothelial epithelial germ line immune stromal
# levels:
}
if(F){
  table(scRNA@meta.data[which(scRNA@meta.data$organ_tissue=="Pancreas"),]$donor,useNA="ifany")
  # TSP1  TSP2  TSP3  TSP4  TSP5  TSP6  TSP7  TSP8  TSP9 TSP10 TSP11 TSP12 TSP13 TSP14 TSP15 
  # 6798     0     0     0     0     0     0     0  6699     0     0     0     0     0     0 
  table(scRNA@meta.data[which(scRNA@meta.data$organ_tissue=="Pancreas"),]$method,useNA="ifany")
  # 10X smartseq2 
  # 12501       996 
  table(scRNA@meta.data[which(scRNA@meta.data$organ_tissue=="Pancreas"),]$anatomical_information)
  # Endocrine Exocrine  exocrine
  table(as.character(scRNA@meta.data[which(scRNA@meta.data$organ_tissue=="Pancreas"),]$free_annotation))
  # 320 Levels: AREG High T-Cell AREG Low T-Cell Adventitial Fibroblast ... vein_endothelial_cell
  table(as.character(scRNA@meta.data[which(scRNA@meta.data$organ_tissue=="Pancreas"),]$cell_ontology_class))
  # 177 Levels: acinar cell of salivary gland adipocyte adventitial cell alveolar fibroblast artery endothelial cell ... vein endothelial cell
  scRNA@meta.data[which(scRNA@meta.data$organ_tissue=="Pancreas"),]$compartment
  write.csv(scRNA@meta.data[which(scRNA@meta.data$free_annotation!=scRNA@meta.data$cell_ontology_class),],file="qwer.csv")
  # scRNA@meta.data$organ_tissue=="Pancreas"&
  # scRNA@meta.data$organ_tissue=="Pancreas"&
  # scRNA@meta
}
Idents(scRNA)
# sce.markers <- FindAllMarkers(object = scRNA, only.pos = TRUE, min.pct = 0.25, 
#                               thresh.use = 0.25)
# scRNA$csfgroup=ifelse(scRNA$group=="HC","HC","PD-DLB")

scRNA@meta.data$celltype=scRNA@meta.data$cell_ontology_class
Idents(scRNA)=scRNA$celltype
# names(Idents(scRNA))==names(scRNA$celltype)
# ct = names(table(scRNA$celltype)) 
ct = names(table(scRNA$celltype)) 
#统计各个亚群PD-DLB vs HC 差异基因
scRNA_Pancreas=subset(scRNA,organ_tissue=="Pancreas")
save(scRNA_Pancreas,file = "GSE201333_Normal_Pancreas.Rdata")
load(file = "GSE201333_Normal_Pancreas.Rdata")
# scRNA_pan2= scRNA[,scRNA@meta.data$organ_tissue %in% c("Pancreas")]
# scRNA= scRNA[, Idents(scRNA) %in% c( "T cell","B cell" )]



table(scRNA_pan$donor)
scRNA$compartment
scRNA_nogem$organ_tissue
scRNA_nogem=subset(scRNA,compartment!="germ line")
scRNA_nogem=subset(scRNA_nogem,organ_tissue%in%unique(wu$scRNA_nogem.meta.data.organ_tissue))
scRNA_nogem@meta.data$compartment=as.factor(scRNA_nogem@meta.data$compartment)
wu=data.frame(scRNA_nogem@meta.data$compartment,scRNA_nogem@meta.data$organ_tissue)
wu=distinct(wu)
wu=subset(wu,!(wu$scRNA_nogem.meta.data.organ_tissue%in%c("Fat","Blood","Bone_Marrow","Kidney","Muscle","Lymph_Node","Spleen")))
unique(wu$scRNA_nogem.meta.data.organ_tissue)
as.factor(wu$scRNA_nogem.meta.data.organ_tissue)
table(wu$scRNA_nogem.meta.data.organ_tissue)
p2=VlnPlot(scRNA_nogem, group.by = "compartment", features = "FTO") + #cell_ontology_class
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()+facet_wrap(vars(scRNA_nogem$organ_tissue),ncol = 1)
p2
ggsave(filename="Vlnplot_normal_FTO.png",plot=p2)
table(scRNA_nogem$compartment)
ggsave(filename="Vlnplot_normal_FTO.png",plot=p2)
all_sig_markers = lapply(ct, function(x){
  x = ct[1]
  print(x)
  markers <- FindMarkers(scRNA_pan, #subset.ident = "fibroblast",
                         ident.1 = 'pancreatic ductal cell',#上位的
                         ident.2 = "fibroblast",#下位的
                         assay = 'RNA',slot = 'counts',
                         logfc.threshold =0,min.pct = 0)
                         # group.by = 'gender')
  markers <- FindMarkers(scRNA_pan,ident.1 ="fibroblast",
                            only.pos = TRUE, 
                            features = "VCAN",
                         assay = 'RNA',slot = 'counts',

  # group.by = 'gender')
  features = "GIMAP4"
  markers["EPCAM",]#6
  markers["KRT19",]#2
  markers["VCAN",]#-1.8
  markers["ACTA2",]#-1.8
  markers["FTO",]#-1.8
  markers_sig <- subset(markers, p_val_adj < 0.05)
  return(markers_sig)
})
ggsave(filename="Vlnplot2.png",plot=p2)

names(all_sig_markers)=ct
save(all_sig_markers,file="DEGs_in_all_celltype.Rdata")

library(Libra)


DE = run_de(hagai_toy)

DE = run_de(hagai_toy, de_family = 'pseudobulk', de_method = 'DESeq2', de_type = 'LRT', n_threads = 16)

DE = run_de(hagai_toy, de_family = 'mixedmodel')

DE = run_de(hagai_toy, de_family = 'mixedmodel', de_method = 'linear', de_type = 'LRT', n_threads = 16)

DE = run_de(hagai_toy, de_family = 'mixedmodel', de_method = 'linear', de_type = 'LRT', n_threads = 8)

DV = calculate_delta_variance(hagai_toy)


cc=Libra::to_pseudobulk(
  scRNA,
  # meta = NULL,
  replicate_col = "donor",#生物学复制
  cell_type_col = "cell_ontology_class",#细胞类型
  label_col = "organ_tissue",#实验标签
  min_cells = 3,
  min_reps = 2,
  min_features = 0#保留基因的最小计数数。默认为0
)

pseudobulk sum
pseudobulk mean
#assay数据提取 
=GetAssayData( scRNA, slot = "counts")#

scRNA= scRNA[,scRNA@meta.data$seurat_clusters %in% c(0,2)]
scRNA= scRNA[, Idents(scRNA) %in% c( "T cell" , "B cell" )]



Idents(scRNA) <- scRNA$Majory_type
subset(x = scRNA, idents = c("CD4 T cells", "CD8 T cells"))
subset(x = scRNA, subset = nFeatures > 500 & PC1 > 5, idents = "B cells")
subset(x = scRNA, subset = orig.ident == "Replicate1")
subset(x = scRNA, downsample = 100)
subset(x = scRNA, features = VariableFeatures(object = scRNA))
scRNA= scRNA[,scRNA@meta.data$seurat_clusters %in% c(0,2)]
scRNA= scRNA[, Idents(scRNA) %in% c( "T cell","B cell" )]



#assay数据提取
GetAssayData( scRNA, slot = "counts")
scRNA<- SetAssayData(scRNA, slot = "scale.data", new.data = new.data)
#embeddings 数据提取
Embeddings(object = scRNA, reduction = "pca")
# FetchData can pull anything from expression matrices, cell embeddings, or metadata
FetchData(object = scRNA, vars = c("PC_1", "percent.mito"))
exprs <- data.frame(FetchData(object = scRNA, vars =  VariableFeatures(object = scRNA)))
exprs <- t(exprs)  #行列变换
write.csv(exprs,file = 'exprs.csv')