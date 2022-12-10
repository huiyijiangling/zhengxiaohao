#create filtered counts
source(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\scripts\filter_f1000.R)")
filter_seq <- function(x){
library(edgeR)
expr = DGEList(counts = x)
keepALL <- rowSums(cpm(expr) > 1) >= 0.5*ncol(x)
nGenes <- as.numeric(summary(keepALL)[2]) + 
  as.numeric(summary(keepALL)[3])
nKeep <- summary(keepALL)[3]
cat (paste('Total Number of genes: ', nGenes, '\n', sep=''))
cat (paste('Number of genes for downstream analysis: ', nKeep, 
           '\n', sep=''))
exprSet=x[,substr(colnames(x),14,15)=="01"][keepALL,]
return(exprSet)
}

#tcga
tcga_filtered_rownames=rownames(filter_seq(rnaCounts))
#ICGC_PACA_AU_seq
ICGC_PACA_AU_seq_filtered_rownames=rownames(filter_seq(ICGC_PACA_AU_seq_raw_dornor))
#ICGC_PACA_AU_array
if(F){
flist=filterEx(ICGC_PACA_AU_array_norm_dornor,0.1,rep("A",ncol(ICGC_PACA_AU_array_norm_dornor)))
# ICGC_PACA_AU_array_norm_dornor=ICGC_PACA_AU_array_norm_dornor[flist,]
}
ICGC_PACA_AU_array_filtered_rownames=rownames(ICGC_PACA_AU_array_norm_dornor)
#ICGC_PACA_CA
ICGC_PACA_CA_seq_raw_dornor[is.na(ICGC_PACA_CA_seq_raw_dornor)] <- 0
ICGC_PACA_CA_seq_filtered_rownames=rownames(filter_seq(ICGC_PACA_CA_seq_raw_dornor))
#EMTAB6134
if(F){
flist=filterEx(genes_expr_mean_EMTAB6134,0.1,rep("A",ncol(genes_expr_mean_EMTAB6134)))
}
#GSE71729
if(F){
flist=filterEx(genes_expr_mean_GSE71729_rna,0.1,rep("A",ncol(genes_expr_mean_GSE71729_rna)))
}
#GSE21501
if(F){
  flist=filterEx(genes_expr_mean_GSE21501_rna,0.1,rep("A",ncol(genes_expr_mean_GSE21501_rna)))
}
#final keep
keep_forhigh=Reduce(intersect,list(tcga_filtered_rownames,ICGC_PACA_CA_seq_filtered_rownames,ICGC_PACA_AU_array_filtered_rownames,ICGC_PACA_AU_seq_filtered_rownames,rownames(genes_expr_mean_EMTAB6134),rownames(genes_expr_mean_GSE71729_rna),rownames(genes_expr_mean_GSE21501_rna)))
save(keep_forhigh,file="keep_forhigh.Rdata")
