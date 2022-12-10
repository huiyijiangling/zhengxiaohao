


uuuuuu=tibble()


eset=as.matrix(probes_expr_GSE32688_mir)

library(genefilter)
eset.filter <- nsFilter(eset, require.entrez=F, remove.dupEntrez=F,var.cutoff=0.5) #不使用函数的合并探针功能，50%很好有文献支持
eset.filter$filter.log #查看每一步筛掉多少探针
eset.filter <- eset.filter$eset #只留下ExpressionSet对象
probes_expr_eset.filter=exprs(eset.filter)
probes_expr_eset.filter=as.data.frame(probes_expr_eset.filter)
pdf("E-MTAB-6134_probes_expr_eset.filter.pdf",width=100)
par(mfrow = c(2,1));
boxplot(probes_expr_eset.filter,las=2)
dev.off()
#选择 注意每次要改 开始修剪了!!!!!!!!!
# probes_expr=probes_expr_rma
probes_expr=probes_expr[rownames(probes_expr) %in% rownames(probes_expr_eset.filter),]





?nsFilter
