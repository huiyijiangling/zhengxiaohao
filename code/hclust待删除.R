导读
当用到多个数据集合并分析时不可避免要处理批次效应，一篇2011年的批次效应处理工具测评文章中说，sva的ComBat函数是表现最好的~所以赶紧学习下它。

1.准备R包
if(!require(BiocManager)) install.packages("BiocManager")
if(!require(sva)) BiocManager::install("sva")
if(!require(bladderbatch)) BiocManager::install("bladderbatch")
library(sva)
library(bladderbatch)
2. 了解数据
示例数据取自bladderbatch包，用data加载，和GEO下载的数据一样，可直接用函数提取表达矩阵和临床信息。

data(bladderdata)
edata <- exprs(bladderEset) 
pheno <- pData(bladderEset) 
edata=rnaCountslog21_stomach
pheno=TcgaTargetGTEX_phenotype_stomach624
pheno$batch <- ifelse(pheno$`_study`=="GTEX",
                                                       "G",
                                                       ifelse(pheno$`_study`=="TCGA"
                                                              &substr(rownames(pheno),14,15)=="01",
                                                              "T",
                                                              ifelse(pheno$`_study`=="TCGA"
                                                                     &substr(rownames(pheno),14,15)=="11",
                                                                     "T",
                                                                     "WRONG")
                                                       )
)
pheno$cancer <- ifelse(pheno$`_study`=="GTEX",
                       "N1",
                       ifelse(pheno$`_study`=="TCGA"
                              &substr(rownames(pheno),14,15)=="01",
                              "T",
                              ifelse(pheno$`_study`=="TCGA"
                                     &substr(rownames(pheno),14,15)=="11",
                                     "N2",
                                     "WRONG")
                       )
)

dim(edata);head(pheno)
## [1] 22283    57
##              sample outcome batch cancer
## GSM71019.CEL      1  Normal     3 Normal
## GSM71020.CEL      2  Normal     2 Normal
## GSM71021.CEL      3  Normal     2 Normal
## GSM71022.CEL      4  Normal     3 Normal
## GSM71023.CEL      5  Normal     3 Normal
## GSM71024.CEL      6  Normal     3 Normal
table(pheno$cancer)
## Biopsy Cancer Normal 
##   9     40      8 
edata[1:4,1:4]
##           GSM71019.CEL GSM71020.CEL GSM71021.CEL GSM71022.CEL
## 1007_s_at    10.115170     8.628044     8.779235     9.248569
## 1053_at       5.345168     5.063598     5.113116     5.179410
## 117_at        6.348024     6.663625     6.465892     6.116422
## 121_at        8.901739     9.439977     9.540738     9.254368
3.设置model（可选）
mod = model.matrix(~as.factor(cancer), data=pheno)
4.校正其实就一步
combat_edata <- ComBat(dat = edata, batch = pheno$batch, mod = mod)
## Found5batches
## Adjusting for2covariate(s) or covariate level(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
5.聚类看下批次效应处理前后对比
## before
dist_mat <- dist(t(edata))
clustering <- hclust(dist_mat, method = "complete")
## after1
dist_mat_combat <- dist(t(combat_edata))
clustering_combat <- hclust(dist_mat_combat, method = "complete")
可视化

pdf(file="a.pdf",width = 120)
par(mfrow = c(2,2))
plot(clustering, labels = pheno$batch)
# pdf(file="a.pdf",width = 30)
plot(clustering, labels = pheno$cancer)
# dev.off()
plot(clustering_combat, labels = pheno$batch)
plot(clustering_combat, labels = pheno$cancer)
dev.off()
作者：小洁忘了怎么分身
链接：https://www.jianshu.com/p/2f20e3ef1183
来源：简书
著作权归作者所有。商业转载请联系作者获得授权，非商业转载请注明出处。