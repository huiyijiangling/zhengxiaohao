library(NMF)
ranks <- 2:6
f = "rank.rdata"
seed <- 2019

if(T){
# nneg(x, method='posneg')#加最低值
# nmf.input=nneg(d, method='min')#正负分解
d=ui_clust[rownames(ui_clust)%in% mhighall25_left,]
nmf.input=nneg(d, method='posneg')

result <- nmf(nmf.input,rank=2:6,
              nrun = 5,
              seed = seed)
}

if(F){
#第二种 scale 不合适
my.algorithm <- function(x, seed, scale.factor=1){
  # do something with starting point
  # ...
  # for example:
  # 1. compute principal components
  pca <- prcomp(t(x), retx=TRUE)
  # 2. use the absolute values of the first PCs for the metagenes
  # Note: the factorization rank is stored in object 'start'
  factorization.rank <- nbasis(seed)
  basis(seed) <- abs(pca$rotation[,1:factorization.rank])
  # use the rotated matrix to get the mixture coefficient
  # use a scaling factor (just to illustrate the use of extra parameters)
  coef(seed) <- t(abs(pca$x[,1:factorization.rank])) / scale.factor
  # return updated data
  return(seed)
}
result1 <-nmf(nmf.input, c(50,100,150), my.algorithm, mixed=TRUE, scale.factor=10,nrun=5)
}

plot(result)#最大下降的前一个

result <- nmf(nmf.input,rank=4,nrun = 10,seed = seed,method = "Brunet")

index <- extractFeatures(result,"max") # 能给出选了哪些关键的基因，可以用这些基因再次聚类
# kkk=unique(unlist(index))
# kkk=kkk[(kkk+81)%in%kkk]

nmf.input3 = d[unique(unlist(index_r)),]
nmf.input2 = nmf.input[unique(unlist(index)),]

# nmf.input2=d[kkk,]
# nmf.input2=nneg(nmf.input2, method='posneg')
#check
# s <- featureScore(result)
# s
# s[rownames(nmf.input2)]
# sample.order <- names(group[order(group)])
# result <- nmf(nmf.input,
#                rank = 4,
#                seed = seed)
# group <- predict(result) # 提出亚型

result2 <- nmf(nmf.input2,
               rank = 4,
               nrun = 10,
               seed = seed)
group <- predict(result2) # 提出亚型
write.csv(group,file="group.csv")
table(group)


## group
##   1   2   3 
## 124  54 340
jco <- c("#2874C5","#EABF00","#C6524A","#868686")

# png(file = "consensusmap.png")
consensusmap(result2,
             labRow = NA,
             labCol = NA,
             annCol = data.frame("cluster"=group[colnames(nmf.input)]),
             annColors = list(cluster=c("1"=jco[1],"2"=jco[2],"3"=jco[3],"4"=jco[4])))
nmf_name_chek=as.data.frame(nmf.input2)

d=d[unique(stringr::str_split(rownames(nmf_name_chek),"\\.",simplify = T)[,1]),]

#顺序问题我解决不了

basismap(result2,
         cexCol = 1.5,
         cexRow = 1,
         Rowv = T,
         scale="none",
         annColors=list(c("#2874C5","#EABF00","#C6524A","#868686")))

result_3=rposneg(result)
index_r <- extractFeatures(result_3,"max") 

rownames(ui_clust[rownames(ui_clust)%in% mhighall25_left,][all_posneg,])%in%unique(stringr::str_split(rownames(as.data.frame(nmf.input[unique(unlist(index)),])),"\\.",simplify = T)[,1])

huluwa=unique(stringr::str_split(rownames(as.data.frame(nmf.input[unique(unlist(index)),])),"\\.",simplify = T)[,1])


consensusmap(result2,
             labRow = NA,
             labCol = NA,
             annCol = data.frame("cluster"=group[colnames(nmf.input)]),
             annColors = list(cluster=c("1"=jco[1],"2"=jco[2],"3"=jco[3],"4"=jco[4])))
basismap(result2,
         cexCol = 1.5,
         subsetRow =extractFeatures(result_3,"max","combine") ,
         cexRow = 1,
         Rowv = F,
         scale="none",
         annColors=list(c("#2874C5","#EABF00","#C6524A","#868686")))
coefmap(result_3)

index_r <- extractFeatures(result_3,"max") 
all_posneg=unique(unlist(index_r))

d=d[all_posneg,]
group <- predict(result_3) # 提出亚型
write.csv(group,file="group.csv")


index_r <- extractFeatures(result,"max") 

table(group)
