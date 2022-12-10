# d
# 
# 提供的需要聚类的数据矩阵，其中列是样本，行是features，可以是基因表达矩阵。
# 
# maxK
# 
# 聚类结果中分类的最大数目，必须是整数。
# 
# reps
# 
# 重抽样的次数
# 
# pItem
# 
# 样品的抽样比例，如 pItem=0.8 表示采用重抽样方案对样本的80%抽样，经过多次采样，找到稳定可靠的亚组分类。
# 
# pFeature
# 
# Feature的抽样比例
# 
# clusterAlg
# 
# 使用的聚类算法，“hc”用于层次聚类，“pam”用于PAM(Partioning Around Medoids)算法，“km”用于K-Means算法，也可以自定义函数。
# 
# title
# 
# 设置生成的文件的路径
# 
# distance	
# 计算距离的方法，有pearson、spearman、euclidean、binary、maximum、canberra、minkowski。
# 
# tmyPal
# 
# 可以指定一致性矩阵使用的颜色，默认使用白-蓝色
# seed
# 
# 设置随机种子。
# 
# plot
# 
# 不设置时图片结果仅输出到屏幕，也可以设置输出为'pdf', 'png', 'pngBMP' 。
# 
# writeTable
# 
# 若为TRUE，则将一致性矩阵、ICL、log输出到CSV文件
# 
# weightsItem
# 
# 样品抽样时的权重
# 
# weightsFeature
# 
# Feature抽样时的权重
# 
# verbose
# 
# 若为TRUE，可输出进度信息在屏幕上
# 
# corUse
# 
# 设置如何处理缺失值：
# 
# all.obs：假设不存在缺失数据——遇到缺失数据时将报错
# 
# everything：遇到缺失数据时，相关系数的计算结果将被设为missing
# 
# complete.obs：行删除
# 
# pairwise.complete.obs：成对删除，pairwisedeletion
# 
# 
# （2）calcICL函数：
# 用法：
# calcICL(res,title='untitled_consensus_cluster',plot=NULL,writeTable=FALSE)
# 
# 参数：
# res
# 
# consensusClusterPlus的结果
# 
# title
# 
# 设置生成的文件的路径
# 
# plot
# 
# 不设置时图片结果仅输出到屏幕，也可以设置输出为'pdf', 'png', 'pngBMP' 。
# 
# writeTable
# 
# 若为TRUE，则将一致性矩阵、ICL、log输出到CSV文件



### 导入数据进行预处理
##使用ALL示例数据
# library(ALL)
# data(ALL)
# d=exprs(ALL)
# d[1:5,1:5]  #共128个样品，12625个探针数据
# 
# # 01005 01010 03002 04006 04007
# # 1000_at 7.597323 7.479445 7.567593 7.384684 7.905312
# # 1001_at 5.046194 4.932537 4.799294 4.922627 4.844565
# # 1002_f_at 3.900466 4.208155 3.886169 4.206798 3.416923
# # 1003_s_at 5.903856 6.169024 5.860459 6.116890 5.687997
# # 1004_at 5.925260 5.912780 5.893209 6.170245 5.615210
# 
# ##对这个芯片表达数据进行简单的normalization，取在各个样品差异很大的那些gene或者探针的数据来进行聚类分析
# mads=apply(d,1,mad)   #计算每个基因的标准差
# d=d[rev(order(mads))[1:5000],]
# 
# #sweep函数减去中位数进行标准化
# d = sweep(d,1, apply(d,1,median,na.rm=T))

#前面我门不管啦


# boxplot(d)
#也可以对这个d矩阵用DESeq的normalization 进行归一化，取决于具体情况

### 一致性聚类
library(IOBR)
PAAD_signature<-readxl::read_excel("C:/Users/zxh/Desktop/subtype/subtype clsut gene.xlsx") 
PAAD_signature=PAAD_signature[,c(2,5,9)]
PAAD_signature<-format_signatures(PAAD_signature)

load("C:/Users/zxh/Desktop/R/paad-tcga-gtex/dataset5.Rdata")
library(ConsensusClusterPlus)
title="Paadclsut"
d <- list()
for(i in length(PAAD_signature)){
d[[i]]=ui[rownames(ui) %in% PAAD_signature[[i]],]}
PAAD_signature<-readxl::read_excel("C:/Users/zxh/Desktop/R/meta collect/icgc_publish/Genomic analyses identify molecular subtypes of pancreatic cancer/有用的表/sam list/genelist2000.xlsx") 
#结果将会输出k从2-6各个情况下的分型情况，聚类的方法用的是 hc ，抽样比例为0.8，最后输出png图
results = ConsensusClusterPlus(d[[3]],maxK=6,reps=100,pItem=0.8,pFeature=1,title=title,
                               clusterAlg='hc',distance="pearson",seed=1262118388,plot='pdf',writeTable=T)
#euclidean bailey
#distance not show in 2011
#distance correlation in moffit
#1-p TCGA 甲基化 binary
#无所谓
#这里设置的maxK=6、reps=50，但是实际上需要更高的reps（如1000）和更高的maxK（如20）

# 生成的图片如下：
# 
# （1）图例：矩阵热图的颜色图例
# 
# 
# 
# （2）k = 2, 3, 4, 5, 6 时的矩阵热图：矩阵的行和列表示的都是样本，一致性矩阵的值按从0（不可能聚类在一起）到1（总是聚类在一起）用白色到深蓝色表示，一致性矩阵按照一致性分类（热图上方的树状图）来排列。树状图和热图之间的长条即分出来的类别。
# 
# 
# 
# 
# 
# （3）一致性累积分布函数（CDF）图：此图展示了k取不同数值时的累积分布函数，用于判断当k取何值时，CDF达到一个近似最大值，此时的聚类分析结果最可靠。即考虑CDF下降坡度小的k值。
# 
# 
# 
# （4）Delta Area Plot：此图展示的是 k 和 k-1 相比CDF曲线下面积的相对变化。当k=2时，因为没有k=1，所以第一个点表示的是k=2时CDF曲线下总面积，而非面积的相对变化值。当k=6时，曲线下面积仅小幅增长，故5为合适的k值。



###查看结果

#ConsensusClusterPlus输出的是一个列表，其中列表对应于来自Kth集群的结果，例如，results[[2]]是k=2的结果。
# View(results)



#输出K=2时的一致性矩阵，consensusMatrix输出一致矩阵。
results[[2]][['consensusMatrix']][1:5,1:5]

# [,1] [,2] [,3] [,4] [,5]
# [1,] 1.0000000 1.0000000 0.9655172 1.0000000 1.0000000
# [2,] 1.0000000 1.0000000 0.8857143 1.0000000 1.0000000
# [3,] 0.9655172 0.8857143 1.0000000 0.9166667 0.8823529
# [4,] 1.0000000 1.0000000 0.9166667 1.0000000 1.0000000
# [5,] 1.0000000 1.0000000 0.8823529 1.0000000 1.0000000

# hclust选项
results[[2]][['consensusTree']]

# Call:
# hclust(d = as.dist(1 - fm), method = finalLinkage)
# Cluster method   : average
# Number of objects: 128

#consensusClass可以看到各个样品被分到了哪个类别里面去
results[[2]][['consensusClass']][1:5]

# 01005 01010 03002 04006 04007
# 1 1 1 1 1

# 计算聚类一致性 (cluster-consensus) 和样品一致性 (item-consensus)
icl <- calcICL(results, title = title,plot = 'pdf')

# 结果：
# （1）Tracking Plot：此图下方的黑色条纹表示样品，展示的是样品在k取不同的值时，归属的分类情况，不同颜色的色块代表不同的分类。取不同k值前后经常改变颜色分类的样品代表其分类不稳定。若分类不稳定的样本太多，则说明该k值下的分类不稳定。
# 
# 
# 
# （2）Cluster-Consensus Plot：此图展示的是不同k值下，每个分类的cluster-consensus value（该簇中成员pairwise consensus values的均值）。该值越高（低）代表稳定性越高（低）。可用于判断同一k值下以及不同k值之间cluster-consensus value的高低。
# 
# 
# 
# （3）item-Consensus Plot：纵坐标代表Item-consensus values。k值不同时，每个样本都会有一个对应不同簇的item-consensus values。竖条代表每一个样本，竖条的高度代表该样本的总item-consensus values。每个样本的上方都有一个小叉叉，小叉叉的颜色代表该样本被分到了哪一簇。从这张图，可以看到每个样本的分类是否足够“纯净”，从而帮助决定k值，例如当k=6时，样本的分类变得没有那么纯净，说明k=5才是合适的。


###############################################################这
## 返回了具有两个元素的list，分别查看一下
dim(icl[['clusterConsensus']])
# [1] 20 3

icl[['clusterConsensus']]

# k cluster clusterConsensus
# [1,] 2 1 0.9079483
# [2,] 2 2 0.7584326
# [3,] 3 1 0.6246200
# [4,] 3 2 0.9111359
# [5,] 3 3 0.9864123
# [6,] 4 1 0.8908356
# [7,] 4 2 0.8869606
# [8,] 4 3 0.6663949
# [9,] 4 4 0.9829523
# [10,] 5 1 0.8612347
# [11,] 5 2 0.8848722
# [12,] 5 3 0.5568284
# [13,] 5 4 0.8390983
# [14,] 5 5 1.0000000
# [15,] 6 1 0.8256498
# [16,] 6 2 0.9377737
# [17,] 6 3 0.6496445
# [18,] 6 4 0.7267928
# [19,] 6 5 0.6982017
# [20,] 6 6 1.0000000

dim(icl[['itemConsensus']]) #128*（2+3+4+5+6）=2560
# [1] 2560 4

icl[['itemConsensus']][1:5,]

# k cluster item itemConsensus
# 1 2 1 01007 0.05261929
# 2 2 1 01003 0.05551604
# 3 2 1 02020 0.04554248
# 4 2 1 04018 0.06059102
# 5 2 1 09002 0.06779347
# 
# 结尾：确定亚型后，接着可以基于各个亚型来分析：比如绘制不同亚型的表达模型热图、看看某个分类下不同亚型的表达高低差异、做不同亚型之间基因表达的显著性差异以及结合PCA或者共表达网络等等。
