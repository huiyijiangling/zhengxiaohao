wula=AverageExpression(object = sce,group.by = c("seurat_clusters"),slot = "counts")#counts
# wula=AverageExpression(object = sce_T,group.by = c("orig.ident","seurat_clusters"),slot = "counts")#counts
wula=as.data.frame(wula)

# sce@meta.data$anno=""
pre_ipmn_primary_meta=readxl::read_excel(r"(C:\Users\zxh\Desktop\R\scCancer\pre_ipmn_primary_meta.xlsx)")                    
levels(sce@meta.data$sample.ident)
meta_new_ori=sce@meta.data


meta_new_rownames=rownames(sce@meta.data)

meta_new=dplyr::left_join(
  sce@meta.data,
  pre_ipmn_primary_meta,
  by=c("sample.ident"="V1")
)
rownames(meta_new)=meta_new_rownames
table(sce@meta.data$nCount_RNA==meta_new$nCount_RNA)
sce@meta.data=meta_new  

# sc epi 
celltype=data.frame(ClusterID=0:31,
                    celltype='Unknown')
celltype[celltype$ClusterID %in% c(0), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(1), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(2), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(3), 2] = 'Epithelial'#mixed
celltype[celltype$ClusterID %in% c(4), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(5), 2] = 'Acinar'
celltype[celltype$ClusterID %in% c(6), 2] = 'Epithelial'#T
celltype[celltype$ClusterID %in% c(7), 2] = 'Acinar'#B
celltype[celltype$ClusterID %in% c(8), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(9), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(10), 2] = 'Epithelial'#pro
celltype[celltype$ClusterID %in% c(11), 2] = 'Epithelial'#Plasma
celltype[celltype$ClusterID %in% c(12), 2] = 'NKT'
celltype[celltype$ClusterID %in% c(13), 2] = 'Epithelial'#中性粒S8 s9，monocyte
celltype[celltype$ClusterID %in% c(14), 2] = 'Epithelial'#mixed
celltype[celltype$ClusterID %in% c(15), 2] = 'B'#
celltype[celltype$ClusterID %in% c(16), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(17), 2] = 'Islet'
celltype[celltype$ClusterID %in% c(18), 2] = 'Epithelial'#
celltype[celltype$ClusterID %in% c(19), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(20), 2] = 'Mast'
celltype[celltype$ClusterID %in% c(21), 2] = 'Tuft'
celltype[celltype$ClusterID %in% c(22), 2] = 'RBC'
celltype[celltype$ClusterID %in% c(23), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(24), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(25), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(26), 2] = 'Epithelial'#T
celltype[celltype$ClusterID %in% c(27), 2] = 'Epithelial'#
celltype[celltype$ClusterID %in% c(28), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(29), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(30), 2] = 'hepatocellular_epithelial'#pro
celltype[celltype$ClusterID %in% c(31), 2] = 'Biliary_epithelial'#

# sc imm 
celltype=data.frame(ClusterID=0:22,
                    celltype='Unknown')
celltype[celltype$ClusterID %in% c(0), 2] = 'γδT'
celltype[celltype$ClusterID %in% c(1), 2] = 'Macrophages'
celltype[celltype$ClusterID %in% c(2), 2] = 'Effector CD8T'
celltype[celltype$ClusterID %in% c(3), 2] = 'TEM'#mixed
celltype[celltype$ClusterID %in% c(4), 2] = 'B'
celltype[celltype$ClusterID %in% c(5), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(6), 2] = 'CD14 Monocytes'#T
celltype[celltype$ClusterID %in% c(7), 2] = 'CD16 Monocytes'#B
celltype[celltype$ClusterID %in% c(8), 2] = 'cDC2'
celltype[celltype$ClusterID %in% c(9), 2] = 'CD16 NK'
celltype[celltype$ClusterID %in% c(10), 2] = 'Treg'#pro
celltype[celltype$ClusterID %in% c(11), 2] = 'CD56 NK'#Plasma
celltype[celltype$ClusterID %in% c(12), 2] = 'Plasma'#mixed
celltype[celltype$ClusterID %in% c(13), 2] = 'Effector Memory CD8T'#中性粒S8 s9，monocyte
celltype[celltype$ClusterID %in% c(14), 2] = 'Macrophages'#mixed
celltype[celltype$ClusterID %in% c(15), 2] = 'B'#
celltype[celltype$ClusterID %in% c(16), 2] = 'Naive Monocytes'
celltype[celltype$ClusterID %in% c(17), 2] = 'CD16 Monocytes'
celltype[celltype$ClusterID %in% c(18), 2] = 'Naive CD4 T'#
celltype[celltype$ClusterID %in% c(19), 2] = 'cDC1'
celltype[celltype$ClusterID %in% c(20), 2] = 'Macrophages'
celltype[celltype$ClusterID %in% c(21), 2] = 'B'
celltype[celltype$ClusterID %in% c(22), 2] = 'Epithelial'

# sc sto 
celltype=data.frame(ClusterID=0:22,
                    celltype='Unknown')
celltype[celltype$ClusterID %in% c(0), 2] = 'Endothelial'
celltype[celltype$ClusterID %in% c(1), 2] = 'Fibroblasts'
celltype[celltype$ClusterID %in% c(2), 2] = 'Fibroblasts'
celltype[celltype$ClusterID %in% c(3), 2] = 'Fibroblasts'#mixed
celltype[celltype$ClusterID %in% c(4), 2] = 'Fibroblasts'
celltype[celltype$ClusterID %in% c(5), 2] = 'Endothelial'
celltype[celltype$ClusterID %in% c(6), 2] = 'Endothelial'#T
celltype[celltype$ClusterID %in% c(7), 2] = 'Epithelial'#B
celltype[celltype$ClusterID %in% c(8), 2] = 'Fibroblasts'
celltype[celltype$ClusterID %in% c(9), 2] = 'Endothelial'
celltype[celltype$ClusterID %in% c(10), 2] = 'Fibroblasts'#pro
celltype[celltype$ClusterID %in% c(11), 2] = 'Fibroblasts'#Plasma
celltype[celltype$ClusterID %in% c(12), 2] = 'Macrophages'#mixed
celltype[celltype$ClusterID %in% c(13), 2] = 'Fibroblasts'#中性粒S8 s9，monocyte
celltype[celltype$ClusterID %in% c(14), 2] = 'NKT'#mixed
celltype[celltype$ClusterID %in% c(15), 2] = 'Epithelial'#
celltype[celltype$ClusterID %in% c(16), 2] = 'B'
celltype[celltype$ClusterID %in% c(17), 2] = 'NKT'
celltype[celltype$ClusterID %in% c(18), 2] = 'Fibroblasts'#
celltype[celltype$ClusterID %in% c(19), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(20), 2] = 'Schwann'
celltype[celltype$ClusterID %in% c(21), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(22), 2] = 'Mixed'


# snuc epi 
celltype=data.frame(ClusterID=0:42,
                    celltype='Unknown')
celltype[celltype$ClusterID %in% c(0), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(1), 2] = ' '
celltype[celltype$ClusterID %in% c(2), 2] = 'Acinar'
celltype[celltype$ClusterID %in% c(3), 2] = ' '#mixed
celltype[celltype$ClusterID %in% c(4), 2] = ' '
celltype[celltype$ClusterID %in% c(5), 2] = 'Acinar'
celltype[celltype$ClusterID %in% c(6), 2] = 'Epithelial'#T
celltype[celltype$ClusterID %in% c(7), 2] = 'Acinar'#B
celltype[celltype$ClusterID %in% c(8), 2] = 'Endocrine'
celltype[celltype$ClusterID %in% c(9), 2] = 'Storma'
celltype[celltype$ClusterID %in% c(10), 2] = 'Epithelial'#pro
celltype[celltype$ClusterID %in% c(11), 2] = 'Epithelial'#Plasma
celltype[celltype$ClusterID %in% c(12), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(13), 2] = 'Endocrine'#中性粒S8 s9，monocyte
celltype[celltype$ClusterID %in% c(14), 2] = 'Epithelial'#mixed
celltype[celltype$ClusterID %in% c(15), 2] = 'Epithelial'#
celltype[celltype$ClusterID %in% c(16), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(17), 2] = 'Endocrine'
celltype[celltype$ClusterID %in% c(18), 2] = 'Epithelial'#
celltype[celltype$ClusterID %in% c(19), 2] = 'Immune'
celltype[celltype$ClusterID %in% c(20), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(21), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(22), 2] = 'Acinar'
celltype[celltype$ClusterID %in% c(23), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(24), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(25), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(26), 2] = 'Epithelial'#T
celltype[celltype$ClusterID %in% c(27), 2] = ' '#
celltype[celltype$ClusterID %in% c(28), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(29), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(30), 2] = ' '#pro
celltype[celltype$ClusterID %in% c(31), 2] = ' '#
celltype[celltype$ClusterID %in% c(32), 2] = ' '
celltype[celltype$ClusterID %in% c(33), 2] = ' '
celltype[celltype$ClusterID %in% c(34), 2] = ' '
celltype[celltype$ClusterID %in% c(35), 2] = ' '
celltype[celltype$ClusterID %in% c(36), 2] = ' '#T
celltype[celltype$ClusterID %in% c(37), 2] = 'Acinar'#
celltype[celltype$ClusterID %in% c(38), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(39), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(40), 2] = 'Epithelial'#pro
celltype[celltype$ClusterID %in% c(41), 2] = ' '#
celltype[celltype$ClusterID %in% c(42), 2] = ' '#
# snuc imm 
celltype=data.frame(ClusterID=0:27,
                    celltype='Unknown')
celltype[celltype$ClusterID %in% c(0), 2] = 'CD8T'
celltype[celltype$ClusterID %in% c(1), 2] = 'Macrophages'
celltype[celltype$ClusterID %in% c(2), 2] = 'Macrophages'
celltype[celltype$ClusterID %in% c(3), 2] = 'CD4T'#mixed
celltype[celltype$ClusterID %in% c(4), 2] = ' '
celltype[celltype$ClusterID %in% c(5), 2] = ' '
celltype[celltype$ClusterID %in% c(6), 2] = 'B' #'CD14 Monocytes'#T
celltype[celltype$ClusterID %in% c(7), 2] = 'cDC2'#B
celltype[celltype$ClusterID %in% c(8), 2] = 'Monocytes'
celltype[celltype$ClusterID %in% c(9), 2] = 'CD16 NK'
celltype[celltype$ClusterID %in% c(10), 2] = 'Epithelial'#pro
celltype[celltype$ClusterID %in% c(11), 2] = 'Fibroblasts'#Plasma
celltype[celltype$ClusterID %in% c(12), 2] = 'Plasma'#mixed
celltype[celltype$ClusterID %in% c(13), 2] = 'Macrophages'#中性粒S8 s9，monocyte
celltype[celltype$ClusterID %in% c(14), 2] = 'Mixed'#mixed
celltype[celltype$ClusterID %in% c(15), 2] = 'Plasma'#
celltype[celltype$ClusterID %in% c(16), 2] = 'Naive Monocytes'
celltype[celltype$ClusterID %in% c(17), 2] = 'Naive T'
celltype[celltype$ClusterID %in% c(18), 2] = 'Naive CD4 T'#
celltype[celltype$ClusterID %in% c(19), 2] = ' '
celltype[celltype$ClusterID %in% c(20), 2] = 'DC'
celltype[celltype$ClusterID %in% c(21), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(22), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(23), 2] = 'CD4T'#中性粒S8 s9，monocyte
celltype[celltype$ClusterID %in% c(24), 2] = 'cDC3'#mixed
celltype[celltype$ClusterID %in% c(25), 2] = 'CD4T'#
celltype[celltype$ClusterID %in% c(26), 2] = 'B'
celltype[celltype$ClusterID %in% c(27), 2] = 'Epithelial'

# snuc sto 
celltype=data.frame(ClusterID=0:22,
                    celltype='Unknown')
celltype[celltype$ClusterID %in% c(0), 2] = 'Endothelial'
celltype[celltype$ClusterID %in% c(1), 2] = 'Fibroblasts'
celltype[celltype$ClusterID %in% c(2), 2] = 'Fibroblasts'
celltype[celltype$ClusterID %in% c(3), 2] = 'Fibroblasts'#mixed
celltype[celltype$ClusterID %in% c(4), 2] = 'Fibroblasts'
celltype[celltype$ClusterID %in% c(5), 2] = 'Endothelial'
celltype[celltype$ClusterID %in% c(6), 2] = 'Endothelial'#T
celltype[celltype$ClusterID %in% c(7), 2] = 'Epithelial'#B
celltype[celltype$ClusterID %in% c(8), 2] = 'Fibroblasts'
celltype[celltype$ClusterID %in% c(9), 2] = 'Endothelial'
celltype[celltype$ClusterID %in% c(10), 2] = 'Fibroblasts'#pro
celltype[celltype$ClusterID %in% c(11), 2] = 'Fibroblasts'#Plasma
celltype[celltype$ClusterID %in% c(12), 2] = 'Macrophages'#mixed
celltype[celltype$ClusterID %in% c(13), 2] = 'Fibroblasts'#中性粒S8 s9，monocyte
celltype[celltype$ClusterID %in% c(14), 2] = 'NKT'#mixed
celltype[celltype$ClusterID %in% c(15), 2] = 'Epithelial'#
celltype[celltype$ClusterID %in% c(16), 2] = 'B'
celltype[celltype$ClusterID %in% c(17), 2] = 'NKT'
celltype[celltype$ClusterID %in% c(18), 2] = 'Fibroblasts'#
celltype[celltype$ClusterID %in% c(19), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(20), 2] = 'Schwann'
celltype[celltype$ClusterID %in% c(21), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(22), 2] = 'Mixed'
??????????????????????????????
table(celltype$celltype)
sce@meta.data$celltype = "NA"
# for(i in 1:nrow(celltype)){
#   sce@meta.data[which(sce@meta.data$RNA_snn_res.1 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
output=sce@meta.data
output$cellname=stringr::str_split(rownames(output),"_",n=2,simplify = T)[,2]
output=output[,c("cellname","celltype")]
write.csv(output,file="pre_ipmn_primary_meta_epi.csv")
write.csv(output,file="pre_ipmn_primary_meta_imm.csv")
write.csv(output,file="pre_ipmn_primary_meta_sto.csv")

a1=data.table::fread(file="pre_ipmn_primary_meta_epi.csv",header = T)
a1$from="epi"
a2=data.table::fread(file="pre_ipmn_primary_meta_imm.csv",header = T)
a2$from="imm"
a3=data.table::fread(file="pre_ipmn_primary_meta_sto.csv",header = T)
a3$from="sto"
aa=rbind(a1,a2,a3)
sce@meta.data$cellname=rownames(sce@meta.data)
meta_new_ori=sce@meta.data
meta_new_rownames=rownames(sce@meta.data)
meta_new=dplyr::left_join(
  sce@meta.data,
  aa,
  by=c("cellname"="cellname")
)
rownames(meta_new)=meta_new_rownames
table(sce@meta.data$nCount_RNA==meta_new$nCount_RNA)
sce@meta.data=meta_new  
# sce@meta.data$=ifelse(sce@meta.data$seurat_clusters==10
table(sce@meta.data$celltype,useNA = "ifany")
SCpubr::do_DimPlot(raster=T,sce,group.by = "from",reduction = "umap",label=T,pt.size = 0.2)#,cols=ggsci::pal_d
# SCpubr::do_DimPlot(raster=T,sce,split.by = "celltype",reduction = "umap",label=T,pt.size = 0.2)#,cols=ggsci::pal_d
SCpubr::do_DimPlot(raster=T,sce,group.by = "celltype",reduction = "umap",label=T,pt.size = 0.2)#,cols=ggsci::pal_d
SCpubr::do_DimPlot(raster=T,sce,split.by = "celltype",reduction = "umap",label=T,pt.size = 0.2)#,cols=ggsci::pal_d
DimPlot(raster=T,sce,reduction = "umap",label=T,pt.size = 0.2,label.size = 6)#,cols=ggsci::pal_d3("category20")(20)
DimPlot(raster=T,sce,group.by = "celltype",reduction = "umap",label=T,pt.size = 0.2,label.size = 6)
# DimPlot(raster=T,sce,group.by = "sample.ident",reduction = "umap",label=T,pt.size = 0.2,label.size = 6)
DimPlot(raster=T,sce,split.by = "Location",reduction = "umap",label=T,pt.size = 0.2,label.size = 6)
kk=sce@meta.data
kk0=kk[kk$Location=="Pancreas",]
kk1=kk[kk$Location=="Liver",]
setdiff(kk0$seurat_clusters,kk1$seurat_clusters)
setdiff(kk1$seurat_clusters,kk0$seurat_clusters)#肝脏特有的
kk2=kk[kk$Location=="Lung",]
setdiff(kk0$seurat_clusters,kk2$seurat_clusters)
setdiff(kk2$seurat_clusters,kk0$seurat_clusters)#肺特有的
kk3=kk[kk$Location=="Seeding",]
setdiff(kk0$seurat_clusters,kk3$seurat_clusters)
setdiff(kk3$seurat_clusters,kk0$seurat_clusters)#seeding特有的 30
  