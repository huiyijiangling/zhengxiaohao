rm(list = ls())  ## 魔幻操作，一键清空~ fpkm等于已经标准化完成了效果同TMM voom
options(stringsAsFactors = F)
library(dplyr)
rnainter=read.table(file="RNA-RNA.txt",sep="\t",stringsAsFactors = F,
                    fill = TRUE,encoding = "UTF-8",header=FALSE)
rnainter_homo=filter(rnainter,V5=="Homo sapiens")
rnainter_homo_sig=filter(rnainter_homo,V10>=0.5)#score0.3-0.7都行
# rnainter_homo_sig=rnainter_homo#score0.3-0.7都行
#1
rnainter_homo_sigmrna=filter(rnainter_homo_sig,V8=='mRNA')#score0.3-0.7都行
rnainter_homo_sig_cernamrna=rnainter_homo_sigmrna[rnainter_homo_sigmrna$V4=='miRNA',]
#2
rnainter_homo_siglnc=filter(rnainter_homo_sig,V8=='lncRNA')#score0.3-0.7都行
rnainter_homo_sig_cernalnc=rnainter_homo_siglnc[rnainter_homo_siglnc$V4=='miRNA',]

load(file="lncRNA_genecodev32.Rdata")

a=as.data.frame(rnainter_homo_sig_cernalnc)
ids=b
ids=ids[ids$gene_name %in%  a$V6,]
ids=unique(ids[,c("gene_id","gene_name")])

a=a[a$V6 %in% ids$gene_name,]
#ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids=ids[order(ids$gene_id,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids=ids[!duplicated(ids$gene_id),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
a=a[a$V6 %in% ids$gene_name,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
rnainter_homo_sig_cernalnc=merge(a ,ids ,by.x ="V6",by.y= "gene_name", all=TRUE,sort=TRUE) 

load(file="mRNA_genecodev32.Rdata")

a=as.data.frame(rnainter_homo_sig_cernamrna)
ids=b
ids=ids[ids$gene_name %in%  a$V6,]
ids=unique(ids[,c("gene_id","gene_name")])
a=a[a$V6 %in% ids$gene_name,]
#ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids=ids[order(ids$gene_id,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids=ids[!duplicated(ids$gene_id),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
a=a[a$V6 %in% ids$gene_name,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
rnainter_homo_sig_cernamrna=merge(a ,ids ,by.x ="V6",by.y= "gene_name", all.x = T,sort=TRUE) 

rnainter_homo_sig_cernalnc=unique(rnainter_homo_sig_cernalnc[,c(3,11)])
rnainter_homo_sig_cernamrna=unique(rnainter_homo_sig_cernamrna[,c(3,11)])

#in gastric cancer 2020.3.2日发现删除下面两句以适应其他肿瘤
load("STAD_GDCRNATOOLS.Rdata")
rnainter_homo_sig_cernalnc=rnainter_homo_sig_cernalnc[rnainter_homo_sig_cernalnc$V2 %in% rownames(mirCounts),]
rnainter_homo_sig_cernamrna=rnainter_homo_sig_cernamrna[rnainter_homo_sig_cernamrna$V2 %in% rownames(mirCounts),]


lnc_inter_target=split(rnainter_homo_sig_cernalnc[,1], rnainter_homo_sig_cernalnc$gene_id)
mrna_inter_target=split(rnainter_homo_sig_cernamrna[,1], rnainter_homo_sig_cernamrna$gene_id)
save(rnainter_homo_sig_cernalnc,rnainter_homo_sig_cernamrna,lnc_inter_target,mrna_inter_target,file = "rnainter_homo_sig_cerna_all.Rdata")
load(file="rnainter_homo_sig_cerna_all.Rdata")
