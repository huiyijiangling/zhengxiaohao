library(utils)
library(estimate)
help(package="estimate")
cos (0.6049872018+0.0001467884* ESTIMATE_score)



转录组数据计算，介绍了转录组count数据如何得到三个score和肿瘤纯度

4.操练起来
# library(utils)
# rforge <- "http://r-forge.r-project.org"
# if(!require("estimate"))install.packages("estimate", repos=rforge, dependencies=TRUE)
# library(estimate)
# #help(package="estimate")
4.1 从count数据计算estimate score
找了TCGA的ACC count数据作为示例数据。如果你想要我的示例数据，请在生信星球公众号后台回复“est766”。用自己的count数据也可以噢。

load("exprSet.Rdata")
exprSet[1:3,1:3]
##         TCGA-OR-A5JP-01A TCGA-OR-A5JG-01A TCGA-OR-A5K1-01A
## MT-RNR2           810396          1190259          1206077
## MT-CO1            579888          1298037          1400198
## MT-ND4            623896           768059          1050890
# 这是曾老板写的函数，转录组数据与芯片数据计算过程不同的地方是platform是illumina。

dat=log2(edgeR::cpm(exprSet)+1)
library(estimate)
estimate <- function(dat,pro){
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(dat,file = input.f,sep = '\t',quote = F)
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")   ## 注意platform
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}
pro='ACC'
scores=estimate(dat,pro)
## [1] "Merged dataset includes 10221 genes (191 mismatched)."
## [1] "1 gene set: StromalSignature  overlap= 139"
## [1] "2 gene set: ImmuneSignature  overlap= 141"
head(scores)
##                  StromalScore ImmuneScore ESTIMATEScore
## TCGA.OR.A5JP.01A    -773.8226  -1143.9749    -1917.7975
## TCGA.OR.A5JG.01A    -878.7773   -685.5286    -1564.3059
## TCGA.OR.A5K1.01A    -663.8511   -360.2218    -1024.0729
## TCGA.OR.A5JR.01A    -931.0601   -344.3306    -1275.3907
## TCGA.OR.A5KU.01A    -925.6045  -1222.4672    -2148.0717
## TCGA.OR.A5L9.01A    -247.1255    404.5509      157.4254
4.2 发现输出结果里没有TumorPurity列
affy芯片输出结果是有这一列的。

我对比了一下15年的那篇NC的方法部分，他们计算使用的是 level 3 RNA-seq profiles (RNAseqV2 normalized RSEM)数据，用estimate包计算了scores，用13年NC文章中的公式计算了肿瘤纯度。

公式是：

# Tumour purity=cos (0.6049872018+0.0001467884 × ESTIMATE score)

不要忘了R语言是个好计算器

TumorPurity = cos(0.6049872018+0.0001467884 * scores[,3])
head(TumorPurity)
## TCGA.OR.A5JP.01A TCGA.OR.A5JG.01A TCGA.OR.A5K1.01A TCGA.OR.A5JR.01A 
##        0.9481360        0.9303738        0.8984081        0.9139941 
## TCGA.OR.A5KU.01A TCGA.OR.A5L9.01A 
##        0.9583367        0.8091481
4.3 拿原文的rsem数据来计算
15年的文章给出了计算结果，我复现一下他的计算。RNAseqV2 normalized RSEM 数据不好找，我是从firehouse找到的，并进行了一些整理，让它变成了规范的表达矩阵。

load("exprSet2.Rdata")
exprSet2[1:3,1:3]
##       TCGA-OR-A5J1-01A TCGA-OR-A5J2-01A TCGA-OR-A5J3-01A
## A1BG           16.3305           9.5987          20.7377
## A1CF            0.0000           0.0000           0.5925
## A2BP1          17.2911           5.6368           8.8876
dat2=log2(exprSet2+1)
scores2=estimate(dat2,pro)
## [1] "Merged dataset includes 10412 genes (0 mismatched)."
## [1] "1 gene set: StromalSignature  overlap= 141"
## [1] "2 gene set: ImmuneSignature  overlap= 141"
head(scores2)
##                  StromalScore ImmuneScore ESTIMATEScore
## TCGA.OR.A5J1.01A   -1161.6834  -524.37956    -1686.0629
## TCGA.OR.A5J2.01A    -569.1191  -765.96407    -1335.0831
## TCGA.OR.A5J3.01A   -1295.4628 -1070.18777    -2365.6506
## TCGA.OR.A5J5.01A   -1710.5108  -918.61256    -2629.1234
## TCGA.OR.A5J6.01A    -730.8294    64.28074     -666.5487
## TCGA.OR.A5J7.01A   -1191.5230 -1013.01517    -2204.5382
TumorPurity2 = cos(0.6049872018+0.0001467884 * scores2[,3])
head(TumorPurity2)
## TCGA.OR.A5J1.01A TCGA.OR.A5J2.01A TCGA.OR.A5J3.01A TCGA.OR.A5J5.01A 
##        0.9367771        0.9175140        0.9669692        0.9761016 
## TCGA.OR.A5J6.01A TCGA.OR.A5J7.01A 
##        0.8741344        0.9606713
我把这个计算结果与15年的NC做了比较，一毛不差，开心。

4.4 比较cpm和rsem结果
我把两个数据处理得到的结果组成一个表格来比较一下：

TumorPurity2 = TumorPurity2[match(names(TumorPurity),names(TumorPurity2))]
compare = cbind(TumorPurity,TumorPurity2)
round(compare,4)[1:10,]
##                  TumorPurity TumorPurity2
## TCGA.OR.A5JP.01A      0.9481       0.9493
## TCGA.OR.A5JG.01A      0.9304       0.9344
## TCGA.OR.A5K1.01A      0.8984       0.9003
## TCGA.OR.A5JR.01A      0.9140       0.9133
## TCGA.OR.A5KU.01A      0.9583       0.9603
## TCGA.OR.A5L9.01A      0.8091       0.8106
## TCGA.OR.A5JQ.01A      0.8303       0.8306
## TCGA.OR.A5K4.01A      0.9829       0.9852
## TCGA.OR.A5JL.01A      0.9416       0.9434
## TCGA.OR.A5LS.01A      0.9748       0.9774
其实count数据里还少了两个基因(有mismatch)。不过计算结果也相差无几咯。非常完美的结果。

5.一点争议
illumina输出结果不带有Tumorpurity列，这是包自身的设置。

在biostars上面看到一个讨论，有人认为estimate score 计算肿瘤纯度的公式是根据Affymetrix的芯片数据得出的，是专门针对芯片数据使用，因此不可以用于转录组。建议只计算出estimate score，用这个分数来代替肿瘤纯度的绝对数值用于后续分析。

原帖讨论见：https://www.biostars.org/p/279853/
  
  然而NC 15年就已经发了这篇文章，五年来没人反对，可以认为人家做的是可用的，用就是了呗。