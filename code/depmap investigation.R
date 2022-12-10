source("./updateName.R")# ncbi与biomart相比多map了200，等于是，ncbi》biomart》org.hs
library(readxl)
probe2gene=read.csv("C:/Users/zxh/Desktop/耐药文献/crisper/文献/两大crispr数据库一致性好/原始数据/StronglySelectiveDependencies.csv",header = T)

probe2gene=as.data.frame(probe2gene)
colnames(probe2gene)="symbol"
forthetree=alias2SymbolUsingNCBImodify(alias=probe2gene$symbol)
probe2gene=cbind(probe2gene,forthetree)
table(is.na(probe2gene$ENSEMBL))
probe2gene=probe2gene[!is.na(probe2gene$ENSEMBL),]#我们不研究没有在ENSEMBL里面没有
writexl::write_xlsx(probe2gene,"StronglySelectiveDependencies ENSG.xlsx")

probe2gene=readxl::read_excel("C:/Users/zxh/Desktop/耐药文献/crisper/文献/两大crispr数据库一致性好/原始数据/我总结的commun.xls")

probe2gene=as.data.frame(probe2gene)
probe2gene=lapply(as.list(probe2gene),na.omit)
commondep=Reduce(intersect,list(probe2gene[[1]],probe2gene[[2]]))

probe2gene=readxl::read_excel("C:/Users/zxh/Desktop/耐药文献/crisper/文献/两大crispr数据库一致性好/原始数据/fuckyou.xls")
probe2gene=commondep
forthetree=alias2SymbolUsingNCBImodify(alias=commondep)
probe2gene=cbind(probe2gene,forthetree)
table(is.na(probe2gene$ENSEMBL))
probe2gene=probe2gene[!is.na(probe2gene$ENSEMBL),]#我们不研究没有在ENSEMBL里面没有

writexl::write_xlsx(probe2gene,"983.xlsx")



#depmap 90%
common_dep_nc1=read.csv("C:/Users/zxh/Desktop/耐药文献/crisper/文献/两大crispr数据库一致性好/原始数据/GeneScores_Broad.csv",header = T,row.names = "gene")
common_dep_nc11=lapply(1:ncol(common_dep_nc1),function(x) rownames(common_dep_nc1[common_dep_nc1[[x]]<quantile(common_dep_nc1[[x]],0.7),]))
common_dep_nc11=Reduce(intersect,common_dep_nc11)
common_dep_nc2=read.csv("C:/Users/zxh/Desktop/耐药文献/crisper/文献/两大crispr数据库一致性好/原始数据/GeneScores_Sanger.csv",header = T,row.names = "gene")

########################

source("C:/Users/zxh/Desktop/R/20201025 pancreatic cancer validation/updateName.R")# ncbi与biomart相比多map了200，等于是，ncbi》biomart》org.hs
library(readxl)
load(file="C:/Users/zxh/Desktop/R/paad-tcga-gtex/iCSDB_PANC_top5.Rdata")
probe2gene=ansT2fenzhi3

probe2gene=as.data.frame(probe2gene)
colnames(probe2gene)="symbol"
forthetree=alias2SymbolUsingNCBImodify(alias=probe2gene$symbol)
probe2gene=cbind(probe2gene,forthetree)
table(is.na(probe2gene$ENSEMBL))
probe2gene=probe2gene[!is.na(probe2gene$ENSEMBL),]#我们不研究没有在ENSEMBL里面没有


probe2gene_mrna_5=probe2gene
save(probe2gene_mrna_5,file="probe2gene_mrna_5.Rdata")
writexl::write_xlsx(probe2gene_mrna_5,"probe2gene_mrna_5.xlsx")

load(file="C:/Users/zxh/Desktop/R/paad-tcga-gtex/iCSDB_PANC_top10.Rdata")
probe2gene=ansT2fenzhi3

probe2gene=as.data.frame(probe2gene)
colnames(probe2gene)="symbol"
forthetree=alias2SymbolUsingNCBImodify(alias=probe2gene$symbol)
probe2gene=cbind(probe2gene,forthetree)
table(is.na(probe2gene$ENSEMBL))
probe2gene=probe2gene[!is.na(probe2gene$ENSEMBL),]#我们不研究没有在ENSEMBL里面没有

# lnchighall %in%probe2gene$ENSEMBL
# mhighall %in%probe2gene$ENSEMBL
# probe2gene_mrna_10=probe2gene[probe2gene$ENSEMBL%in%mhighall,]
probe2gene_mrna_10=probe2gene
save(probe2gene_mrna_10,file="probe2gene_mrna_10.Rdata")
writexl::write_xlsx(probe2gene_mrna_10,"probe2gene_mrna_10.xlsx")


# probe2gene=unique(probe2gene[,c("probe_id","ENSEMBL")])#bioc
probe2gene=unique(probe2gene[,c("ID","ENSEMBL")])#soft
# probe2gene=unique(probe2gene[,c("ENSEMBL","probe_id")])
colnames(probe2gene)=c("probe_id","symbol")
genes_expr <- filterEM(probes_expr,probe2gene)
