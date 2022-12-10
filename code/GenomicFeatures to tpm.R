library("GenomicFeatures")
## 导入gff3文件
txdb <- makeTxDbFromGFF("gencode.v22.annotation.gtf.gz",format="gtf")#makeTxDbFromGFF("gencode.v36.annotation.gtf.gz",format="gtf")
## 获取外显子位置
exons_gene <- exonsBy(txdb, by = "gene")
## 去除外显子重叠部分，计算外显子长度
exons_gene_len <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
exons_gene_len2 <- do.call(rbind,lapply(exons_gene_len, data.frame))
exons_gene_len2$gene_id=stringr::str_split(rownames(exons_gene_len2),"\\.",simplify=T)[,1]
exons_gene_len3=exons_gene_len2[!duplicated(exons_gene_len2$gene_id),]
rownames(exons_gene_len3)=exons_gene_len3$gene_id
exons_gene_len3$eff_length=exons_gene_len3$X..i..
exons_gene_len3=exons_gene_len3[,-1]
exons_gene_len=exons_gene_len3
save(exons_gene_len,file="exonlengthv22.Rdata")






kb <- mycounts$Length / 1000
write.csv(kb, file="kb.csv")
##筛选出样本的表达数据
# 这一步骤是筛选出所有样本的表达量数据，第一列是length，所以不要。
countdata <- mycounts[,2:1000]
##计算每个值的rpk
rpk <- countdata / kb

##计算TPM
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
head(tpm)