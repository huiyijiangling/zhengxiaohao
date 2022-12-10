library(AnnoProbe)
library(ggpubr)
library(GEOquery)
# check the ExpressionSet

eSet=gset#eSet=gset[[1]]#如果是list的状态,取expression

# extract the expression matrix and phenotype data
probes_expr <- exprs(eSet);dim(probes_expr)
head(probes_expr[,1:4])
boxplot(probes_expr,las=2)
# probes_expr=log2(probes_expr+1)
# boxplot(probes_expr,las=2)
## pheno info
phenoDat <- pData(eSet)
head(phenoDat[,1:4])

## check GPL and annotate the probes to genes.
(gpl=eSet@annotation)
checkGPL(gpl)
printGPLInfo(gpl)
head(probe2gene)
probe2gene=read.table(file="./GPL13607.soft",sep="\t",stringsAsFactors = F,
           fill = TRUE,encoding = "UTF-8",skip=3516,comment.char = "!",header=T,quote ="")
#probe2gene=idmap(gpl,type = 'soft')#原
probe2gene=probe2gene[,c("ID","GeneName")]
colnames(probe2gene)=c("probeid","symbol")
genes_expr <- filterEM(probes_expr,probe2gene)
ooooooooooooooooooooooooooooooooooo=as.data.frame(rownames(genes_expr))
head(genes_expr)
#注意我们需要通过boxplot看看是否需要进行标准化。而如果要用包就需要raw count。这个集合不用








# do DEG
## define the group
group_list=factor(c(rep('Control',3),rep('Diabetes',3)))
table(group_list)
library(limma)
design=model.matrix(~factor(group_list))
design
fit=lmFit(genes_expr,design)
fit=eBayes(fit)
DEG=topTable(fit,coef=2,n=Inf)
head(DEG)

## visualization
need_deg=data.frame(symbols=rownames(DEG), logFC=DEG$logFC, p=DEG$P.Value)
deg_volcano(need_deg,1)
deg_volcano(need_deg,2)

deg_heatmap(DEG,genes_expr,group_list)
deg_heatmap(DEG,genes_expr,group_list,30)

check_diff_genes('PLCE1',genes_expr,group_list)
check_diff_genes('MPP6',genes_expr,group_list)
library(KEGGREST)
cg <- KEGGREST::keggGet("hsa03410")[[1]]$GENE
cg=as.character(sapply(cg[seq(2,length(cg),by=2)], function(x) strsplit(x,';')[[1]][1]))
check_diff_genes( cg ,genes_expr,group_list)
