# 自己安装  mojaveazure/seurat-disk 这个GitHub包：
#remotes::install_github("mojaveazure/seurat-disk")

library(SeuratDisk)
library(patchwork)
#～～～～～开始读数据～～～～～
##h5ad是python的Scanpy读取文件格式，需要转换
#～～～～读取adipose～～～～
Convert(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\GSE202051_adata_010nuc_10x.h5ad)", "h5seurat",
        overwrite = TRUE,assay = "RNA")
scRNA <- LoadH5Seurat(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\GSE202051_adata_010nuc_10x.h5seurat)",meta.data = TRUE,misc = F,assays = "RNA")
#pdx 看不到
Convert(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\GSE202051_adata_010orgCRT_10x.h5ad)", "h5seurat",
        overwrite = TRUE,assay = "RNA")
scRNA <- LoadH5Seurat(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\GSE202051_adata_010orgCRT_10x.h5seurat)",meta.data = TRUE,misc = F,assays = "RNA")
#pdx 看不到
Convert(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\GSE202051_totaldata-final-toshare.h5ad)", "h5seurat",
        overwrite = TRUE,assay = "RNA",)
GSE202051_totaldatatoshare <- LoadH5Seurat(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\GSE202051_totaldata-final-toshare.h5seurat)",meta.data = TRUE,misc = F,assays = "RNA")

save(GSE202051_totaldatatoshare,file= "GSE202051_totaldatatoshare.rdata")
load(file= r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051_totaldatatoshare.rdata)")
sce=GSE202051_totaldatatoshare

sce.all.list <- SplitObject(sce , split.by = "pid")
sce.all.list
length(names(sce.all.list))
GSE202051_meta=sce@meta.data

colnames(sce)==rownames(sce@meta.data)


setwd(r"(C:\Users\zxh\Desktop\R\scCancer\GSE202051\)")
library(foreach)
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
load(r"(C:\Users\zxh\Desktop\R\meta collect\ensg_symbol_relation.Rdata)")

# stopCluster(cl)
foreach(i=names(sce.all.list))%dopar%{
  epi_mat = sce.all.list[[i]]@assays$RNA@counts
  if(T){
    epi_phe = sce.all.list[[i]]@meta.data
    epi_mat=as.matrix(epi_mat)
    rownames(epi_phe)==colnames(epi_mat)
    epi_mat=expm1(epi_mat)
    epi_mat=t(epi_mat)*(epi_phe$n_counts)/10000#ori tp10k
    epi_mat=t(epi_mat)
    table(epi_mat["FTO",])
    epi_mat=round(epi_mat)#手动检查后再用
}
  ct=as.matrix(epi_mat)
  dir.create(file.path(i))
  dir.create(file.path(paste0(i,'/filtered_feature_bc_matrix/')))#filtered_feature_bc_matrix 
  # feature=data.frame(rownames(ct))
  # colnames(feature)="gene_name"
  # feature=merge(feature,ensg_symbol_trans,by="gene_name",all.x=T)
  # feature[is.na(feature$id),]$id=feature[is.na(feature$id),]$gene_name
  # feature$id=gsub("\\.","",feature$id)
  # feature$gene_type="Gene Expression"
  # data.table::fwrite(feature,file = paste0(i,'/filtered_feature_bc_matrix/features.tsv.gz'),
  #                    quote = F,sep = '\t',
  #                    col.names = F,row.names = F)
  data.table::fwrite(data.frame(rownames(ct),rownames(ct),"Gene Expression"),file = paste0(i,'/filtered_feature_bc_matrix/features.tsv.gz'),
                     quote = F,sep = '\t',
                     col.names = F,row.names = F)
  data.table::fwrite(data.frame(colnames(ct)),file = paste0(i,'/filtered_feature_bc_matrix/barcodes.tsv.gz'),quote = F,sep = '\t',col.names = F,row.names = F)
  file=paste0(i,'/filtered_feature_bc_matrix/matrix.mtx')
  sink(file)
  cat("%%MatrixMarket matrix coordinate integer general\n")
  cat("%\n")
  cat(paste(nrow(ct),ncol(ct),sum(ct>0),"\n")) 
  sink()
  #再写入表达量信息
  tmp=data.table::rbindlist(lapply(1:ncol(ct),function(i){
    return(data.frame(row=1:nrow(ct),
                      col=i,
                      exp=ct[,i]))
  }))
  tmp=tmp[tmp$exp>0,]
  head(tmp)
  data.table::fwrite(tmp,file = paste0(i,'/filtered_feature_bc_matrix/matrix.mtx'),quote = F,sep = '\t',
                     col.names = F,row.names = F,append = T)
  R.utils::gzip(file = paste0(i,'/filtered_feature_bc_matrix/matrix.mtx'),remove=TRUE,overwrite=T)
}
stopCluster(cl)



# GSE201333
Convert('./GSE201333/GSM6058681_TabulaSapiens.h5ad', "h5seurat",
        overwrite = TRUE,assay = "RNA")
scRNA <- LoadH5Seurat("GSE201333/GSM6058681_TabulaSapiens.h5seurat",meta.data = TRUE,misc = F)
GSE202051_seurat=scRNA

# ooooo <- LoadH5Seurat("GSE201333/GSM6058681_TabulaSapiens.h5seurat",meta.data = TRUE,misc = T)

#


