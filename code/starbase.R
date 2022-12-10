rm(list = ls())  ## 魔幻操作，一键清空~ fpkm等于已经标准化完成了效果同TMM voom
options(stringsAsFactors = F)
library(dplyr)
library(data.table)
data.table::setDTthreads(threads =4 )#, percent = 99
mm <- data.table::fread(file ="mm.txt",header = T,sep = "\t",
                        verbose = F, integer64 = 'numeric')
table(mm$geneType)
mm=filter(mm,geneType=="protein_coding")
lm <- data.table::fread(file ="lm.txt",header = T,sep = "\t",
                        verbose = F, integer64 = 'numeric')
table(lm$geneType)
library(miRBaseConverter)
miRNANames=lm$miRNAname
version=checkMiRNAVersion(miRNANames, verbose = FALSE)
version#21
miRNANames=mm$miRNAname
version=checkMiRNAVersion(miRNANames, verbose = FALSE)
version#21
rnainter_homo_sig_cernalnc=as.data.frame(unique(lm[,c("miRNAname","geneID")]))
rnainter_homo_sig_cernamrna=as.data.frame(unique(mm[,c("miRNAname","geneID")]))
save(rnainter_homo_sig_cernalnc,rnainter_homo_sig_cernamrna,file = "starbase_homo_sig_cerna_all.Rdata")
#理论上用ENSG更准确，但是因为普通的array转换可能会导致错误，故不转换
#in gastric cancer 2020.3.2日发现删除下面两句以适应其他肿瘤
# load("STAD_GDCRNATOOLS.Rdata")
# rnainter_homo_sig_cernalnc=rnainter_homo_sig_cernalnc[rnainter_homo_sig_cernalnc$V2 %in% rownames(mirCounts),]
# rnainter_homo_sig_cernamrna=rnainter_homo_sig_cernamrna[rnainter_homo_sig_cernamrna$V2 %in% rownames(mirCounts),]

#2020.12.18 重叠cerna时
# rnainter_homo_sig_cernalnc=rnainter_homo_sig_cernalnc[rnainter_homo_sig_cernalnc$miRNAname %in% common_mir,]
# rnainter_homo_sig_cernamrna=rnainter_homo_sig_cernamrna[rnainter_homo_sig_cernamrna$miRNAname %in% common_mir,]
# rnainter_homo_sig_cernalnc=rnainter_homo_sig_cernalnc[rnainter_homo_sig_cernalnc$geneID %in% common_ensg_v36$gene_id,]
# rnainter_homo_sig_cernamrna=rnainter_homo_sig_cernamrna[rnainter_homo_sig_cernamrna$geneID %in% common_ensg_v36$gene_id,]

#按理说应该重叠需要分析的RNA
#lnc_inter_target=split(rnainter_homo_sig_cernalnc[,1], as.character(rnainter_homo_sig_cernalnc$geneID))
#mrna_inter_target=split(rnainter_homo_sig_cernamrna[,1], as.character(rnainter_homo_sig_cernamrna$geneID))