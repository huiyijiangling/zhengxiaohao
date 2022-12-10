library(DESeq2)
rnaCounts_Pancreas=round(rnaCounts_Pancreas)
coldata <- data.frame(TcgaTargetGTEX_phenotype_Pancreas624$sample_type)
colnames(coldata)="groups"
dds <- DESeqDataSetFromMatrix(countData = as.matrix(rnaCounts_Pancreas), colData = coldata, 
                              design = ~groups)
dds$groups <- relevel(dds$groups, ref = "SolidTissueNormal")
dds <- DESeq(dds)
sizeFactors(dds)
res <- results(dds)
knitr::kable(head(res))
res <- as.data.frame(res)
res <- cbind(rownames(res), res)
colnames(res) <- c("gene_id", "baseMean", "log2FoldChange", 
                   "lfcSE", "stat", "pval", "padj")
