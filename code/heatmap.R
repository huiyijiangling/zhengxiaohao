EXPR=rnaExpr
table(group_list)
deg=nrDEG3
x=deg$logFC
names(x)=rownames(deg)
class(x)
x
cg=c(names(head(sort(x),100)),names(tail(sort(x),100)))
head(cg)
library(pheatmap)
pheatmap([cg,],show_colnames =F,show_rownames = F)
n=t(scale(t([cg,])))#scale是对列进行操作，使基因在不同样本可比较，标准化。
# 下面两行可以尝试不运行和运行后的区别~
# 请注意~
n[n>2]=2
n[n<-2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac=data.frame(g=group_list)
rownames(ac)=colnames(n)

pheatmap(n,show_colnames =F,show_rownames = F, cluster_cols = T,
         annotation_col=ac,filename = "7.pheatmap_group.png",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50) # 增加color
)
dev.off()