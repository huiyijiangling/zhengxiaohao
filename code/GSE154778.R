## 系统报错改为英文
Sys.setenv(LANGUAGE="en")
## 禁止转化为因子
options(stringsAsFactors = F)
## 清空环境
rm(list = ls())
## 安装包
library(remotes)
url='https://gitee.com/jmzeng/GEOmirror.git'
class(url)
install_git(url)
## 加载r包
library(GEOquery)
library(GEOmirror)
library(dplyr)
library(data.table)
##设置工作目录
getwd()
# setwd("F:/002DEG")
##下载数据getGPL=F这个代表不下载平台注释文件，因为有时候网络不稳定，后面我们会在网页中下载，然后读取。
gest = getGEO('GSE154778',destdir = ".",getGPL = T)
save(gest,file = "GSE154778_gset.Rdata")
## 获取ExpressionSet对象，包括的表达矩阵和分组信息
expr <- gest[["GSE154778_series_matrix.txt.gz"]]@assayData[["exprs"]]
## fread函数读取txt文档更快更友好。GPL96的数据平台，用GPL96的探针
anno=fread("GPL96-57554.txt",sep = "\t",header=T,data.table=F)
## 看一下anno的列名，方便取出ID和Gene Symbpl
colnames(anno)
## 取出2列的两种方法
gene=anno[,c(1,11)]
gene.1 <- anno[,c("ID","Gene Symbol")]
## 取出Gene Symbol 方便将一个探针对应多个基因名的给分隔开。
x <- gene$`Gene Symbol`
## 字符串切割
a1=strsplit(x,split = "///",fixed = T)
## 把分割后的字符串的第一个元素提取出来，合并成为一个新的向量
## sapply 对每列，保留第一个数。
gene.all <- sapply(a1,function(x){x[1]})
## 向量是有序的排列，所以可以直接合并。
a3=data.frame(anno$ID,gene.all)
###到此为止，我们将注释好的gene symbol 文件与表达矩阵都已经准备好了。
###进行探针的转化
###将矩阵转化为数据框的目的是为了更好的对表达矩阵进行操作。
exp <- as.data.frame(expr)
###
###用merge函数探针转化的两种方法。1是根据列数，2是根据列名
exp1=merge(x=a3,y=exp,by.x = 1,by.y = 0)
exp1 <- merge(x=a3,y=exp,by.x = "anno.ID",by.y = 0)
###整理表达矩阵
##这里根据报错一步一步往下走，其实可以不用报错直接出来，但是还是先学习一下。
rownames(exp1) = exp1$gene.all
###报错的原因就是有重复的基因名，多个探针对应同一个基因名。
###对列去重复
exp2 <- distinct(exp1,gene.all,.keep_all=T)
rownames(exp2) <- exp2$gene.all
###发现还是报错 报错的原因是有NA 缺失值
exp3 <- na.omit(exp2)
###13436 obs 变成了13435obs
rownames(exp3)=exp3$gene.all
###去除第一列和第二列
exp4 <- exp3[,-c(1,2)]
write.table(exp4,file = "exp4.txt",quote = F,sep = '\t',row.names = T,col.names = T)