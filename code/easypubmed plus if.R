#lengths list长度
rm(list=ls())
options(stringsAsFactors = F)
gc()
library(clusterProfiler)
library(org.Hs.eg.db)
filenames_1 <- list.files('./datatmp/',pattern = ".xls",full.name=TRUE)
filenames_1
resultALL1 <- lapply(filenames_1, function(fl) 
  read.delim(fl,header = F))
gene_commonlist=Reduce(intersect,resultALL1)#0
for (i in 1:length(resultALL1)){
  names(resultALL1)[i] <- stringr::str_split(filenames_1,"/",simplify = T)[,3][i]
}


dfd <- list()
for (i in 1:length(resultALL1)){
dfd[i] <- bitr(resultALL1[[i]][[1]], fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)[2]
}
#直接用
ck <- compareCluster(dfd,OrgDb='org.Hs.eg.db',fun = "enrichKEGG",pvalueCutoff=1)
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
#再造变量
# if(F){
mydf=cbind(unlist(dfd),unlist(lapply(1:length(dfd),function(x) rep(names(resultALL1[x] ),length(dfd[[x]] )))))
mydf=as.data.frame(mydf)


xx.formula <- compareCluster(V1~V2, data=mydf,
                             fun='enrichGO', 
                             OrgDb='org.Hs.eg.db',
                             keyType = "ENTREZID",
                             ont = "ALL",
                             pvalueCutoff = 0.05,
                             maxGSSize = 500,
                             minGSSize = 1,
                             readable = T,
                             qvalueCutoff = 0.2)

enrichplot::dotplot(xx.formula, split="ONTOLOGY",showCategory = 5000,font.size =6,
                    x="V2")+  facet_grid(ONTOLOGY~., scale="free") 
ggsave('gene_up_all_GO_dotplot.pdf',width = 8,height=40,dpi = 600)
write.csv(xx.formula@compareClusterResult,'gene.up.all.go.csv',quote = T)


xx.formula <- compareCluster(V1~V2, data=mydf,minGSSize = 1,
                             fun='enrichKEGG', 
                             organism = "hsa",
                             pvalueCutoff = 0.05,
                             maxGSSize = 500,
                             qvalueCutoff = 0.2)
enrichplot::dotplot(xx.formula, showCategory = 1500,font.size =6)
ggsave('gene_up_all_KEGG_dotplot.pdf',width = 16,height=12,dpi = 600)
kk=setReadable(xx.formula, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(kk,'kk.up.all.csv',quote = T)

# pdf(file='upset_mrna_1.pdf',height = 6,width = 8,onefile = T)
# upset(fromList(mrna_1), #输入的数据集
#       sets = c("Gene_WGCNA","Gene_down","Gene_up","Gene_predict"),#,"Survival"
#       # sets = c("Gene_WGCNA","Gene_up","Gene_down","Gene_predict"),
#       # sets = c("PITA", "miRanda","TargetScan"),
#       nsets = 4, #想要可视化的数据集数量,也可以用sets选项自定义
#       nintersects = 20, #要绘制的交点数量 2^n
#       keep.order = T, #是否保持输入顺序,否则按各集合大小排序
#       main.bar.color = 'black', #主条图的颜色
#       mainbar.y.label = 'Number', #y轴标签
#       sets.bar.color = 'blue', #设置集合条图的颜色
#       sets.x.label = 'Number of elements', #设置集合调图的坐标轴标签
#       mb.ratio = c(0.7,0.3), #条形图点图的比例
#       order.by = c("freq","degree"), #交集如何排序,
#       decreasing = c(TRUE,TRUE),
#       # 以上排序是否降序,FALSE
#       # boxplot.summary = c('ReleaseDate','Comedy') #添加箱线图展示数据分布,最多展示两个数据集
#       queries = list(list(query = intersects,
#                           params = list("Gene_up","Gene_predict", "Gene_WGCNA"),#,"Survival"
#                           color = 'red',
#                           active = T)
#                      ,
#                      list(query = intersects,
#                           params = list("Gene_down","Gene_predict", "Gene_WGCNA"),#,"Survival"
#                           color = 'green',
#                           active = T)
#       )
# )
# dev.off()

# dami_query=paste(ab1," ",collapse = " ")
dami_query='"Human Genome Structural Variation Consortium"'
library(easyPubMed)
library(parallel)
library(foreach)
library(doParallel)
library(readxl)
if(T){
impactfactor=read_excel("C:\\Users\\zxh\\Desktop\\2020-IFs ALL.xlsx",sheet=1,col_names=T)
impactfactor$`Journal Impact Factor`=as.numeric(impactfactor$`Journal Impact Factor`)
table(impactfactor$`Journal Impact Factor`>500)
impactfactor$`JCR Abbreviated Title`=gsub("[[:punct:]]"," ", impactfactor$`JCR Abbreviated Title`)
impactfactor$`Full Journal Title`=gsub("[[:punct:]]"," ",impactfactor$`Full Journal Title`)
impactfactor$`JCR Abbreviated Title`=toupper(impactfactor$`JCR Abbreviated Title`)
impactfactor$`Full Journal Title`=toupper(impactfactor$`Full Journal Title`)
}

# dami_query[2]
# dami_query  <-c('"vcan" AND "gastric cancer"','(EMB AND embb306 AND susceptibility AND Mycobacterium) OR (tuberculosis embb306)'                )

dami_query  <-'("ca cancer j clin"[Journal]) AND (2020:2021[pdat])'
dami_on_pubmed <- get_pubmed_ids(dami_query)
dami_abstracts_xml <- fetch_pubmed_data(dami_on_pubmed,encoding="ASCII",retmax = 4999)
dami_abstracts_list <- articles_to_list(dami_abstracts_xml,encoding="ASCII")

if(T){
  cl <- makeCluster(8)
  registerDoParallel(cl)
  fullDF <- tryCatch(
    {foreach(x=dami_abstracts_list, 
             .packages = 'easyPubMed',
             .combine = rbind) %dopar% article_to_df(pubmedArticle = x, 
                                                     autofill = T, 
                                                     max_chars = 5000, 
                                                     getKeywords = T, 
                                                     getAuthors = F)}, 
    error = function(e) {NULL},
    finally = {stopCluster(cl)})
}
fullDF$jabbrv=gsub("[[:punct:]]"," ",fullDF$jabbrv )
fullDF$journal=gsub("[[:punct:]]"," ",fullDF$journal)
fullDF$jabbrv=toupper(fullDF$jabbrv )
fullDF$journal=toupper(fullDF$journal)

fullDF1=merge(fullDF,impactfactor,by.x = "jabbrv",by.y = "Full Journal Title",all.x =  T)
fullDF1=fullDF1[,c(1:14,17)]
fullDF2=merge(fullDF,impactfactor,by.x = "jabbrv",by.y = "JCR Abbreviated Title",all.x =T)
fullDF2=fullDF2[,c(1:14,17)]
fullDF3=merge(fullDF,impactfactor,by.x = "journal",by.y = "Full Journal Title",all.x =T)
fullDF3=fullDF3[,c(1:14,17)]
fullDF4=merge(fullDF,impactfactor,by.x = "journal",by.y = "JCR Abbreviated Title",all.x =T)
fullDF4=fullDF4[,c(1:14,17)]
fullDFsum=data.table::rbindlist(list(fullDF1,fullDF2,fullDF3,fullDF4),use.names = F)
fullDFsum=fullDFsum[order(fullDFsum$`Journal Impact Factor`,decreasing = T,na.last = T),]
fullDFsum=subset(fullDFsum,!duplicated(fullDFsum$pmid))
writexl::write_xlsx(fullDFsum, "fullDFsum.xlsx")

# x=paste(ab2," AND cancer") #太长了

# split(x, sort(rep_len(1:n, length(x))))
####################

multiple_easypubmed <- function(x){
  dami_query=x
  dami_on_pubmed <- get_pubmed_ids(dami_query)
  dami_abstracts_xml <- fetch_pubmed_data(dami_on_pubmed,encoding="ASCII",retmax = 4999)
  dami_abstracts_list <- articles_to_list(dami_abstracts_xml,encoding="ASCII")
  # cl <- makeCluster(4)
  # registerDoParallel(cl)
  fullDF <- tryCatch(
    {foreach(x=dami_abstracts_list, 
             .packages = 'easyPubMed',
             .combine = rbind) %dopar% article_to_df(pubmedArticle = x, 
                                                     autofill = T, 
                                                     max_chars = 5000, 
                                                     getKeywords = T, 
                                                     getAuthors = F)}#, 
    # error = function(e) {NULL}   
    # ,    finally = {stopCluster(cl)}
  )
fullDF$jabbrv=gsub("[[:punct:]]"," ",fullDF$jabbrv )
fullDF$journal=gsub("[[:punct:]]"," ",fullDF$journal)
fullDF$jabbrv=toupper(fullDF$jabbrv )
fullDF$journal=toupper(fullDF$journal)

fullDF1=merge(fullDF,impactfactor,by.x = "jabbrv",by.y = "Full Journal Title",all.x =  T)
fullDF1=fullDF1[,c(1:14,17)]
fullDF2=merge(fullDF,impactfactor,by.x = "jabbrv",by.y = "JCR Abbreviated Title",all.x =T)
fullDF2=fullDF2[,c(1:14,17)]
fullDF3=merge(fullDF,impactfactor,by.x = "journal",by.y = "Full Journal Title",all.x =T)
fullDF3=fullDF3[,c(1:14,17)]
fullDF4=merge(fullDF,impactfactor,by.x = "journal",by.y = "JCR Abbreviated Title",all.x =T)
fullDF4=fullDF4[,c(1:14,17)]
fullDFsum=data.table::rbindlist(list(fullDF1,fullDF2,fullDF3,fullDF4),use.names = F)
fullDFsum=fullDFsum[order(fullDFsum$`Journal Impact Factor`,decreasing = T,na.last = T),]
fullDFsum=subset(fullDFsum,!duplicated(fullDFsum$pmid)) 
fullDFsum$dami_query=x
return(fullDFsum)
}
#1
ab2=unlist(resultALL1)
ab2=unique(ab2)
x=ab2
qwer=lapply(split(x, rep_len(1:180, length(x))),function(x) paste(unlist(x),"", collapse = "OR "))
qwer=lapply(qwer,function(x) paste(unlist(x),"AND cancer"))
qwer=unlist(qwer)
#2
qwer=lapply(ab2,function(x) paste(x,"AND cancer"))

library(doParallel) #
cl <- makeCluster(30)
registerDoParallel(cl)
fullDFsum_ALL <- foreach(i=c("fuckyouyouyou","VCAN AND cancer"),.packages = c(library(easyPubMed),
                            library(parallel),
                            library(foreach),
                            library(doParallel),
                            library(readxl))) %dopar% multiple_easypubmed(i)
stopCluster(cl)#.combine=rbind
library(data.table)
save(fullDFsum_ALL,file="fullDFsum_ALL.Rdata")
load(file="fullDFsum_ALL.Rdata")
fullDFsum_ALL_unlist=data.table::rbindlist(fullDFsum_ALL,use.names = F)
fullDFsum_ALL_unlist=dplyr::distinct(fullDFsum_ALL_unlist)
fullDFsum_ALL_unlist=subset(fullDFsum_ALL_unlist,fullDFsum_ALL_unlist$`Journal Impact Factor`>14)
fullDFsum_ALL_unlist$year=as.numeric(fullDFsum_ALL_unlist$year)
fullDFsum_ALL_unlist=subset(fullDFsum_ALL_unlist,fullDFsum_ALL_unlist$year>=2010)
writexl::write_xlsx(fullDFsum_ALL_unlist, "fullDFsum_ALL.xlsx")
# 循环返回NA
#Error in multiple_easypubmed(i) : task 1 failed - "参数长度为零" 尚无法解决