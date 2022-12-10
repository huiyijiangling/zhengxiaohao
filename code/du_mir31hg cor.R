

TcgaTargetGTEX_phenotype_Pancreas624=subset(TcgaTargetGTEX_phenotype_Pancreas624,TcgaTargetGTEX_phenotype_Pancreas624$sample_type=="PrimaryTumor")


common_sample=Reduce(intersect,list(rownames(TcgaTargetGTEX_phenotype_Pancreas624),colnames(rnaCounts),colnames(mirCounts)))

rnaCounts=as.data.frame(rnaCounts)
rnaCounts=subset(rnaCounts,select=common_sample)
rnaCounts=log2(rnaCounts+1)
mirCounts=as.data.frame(mirCounts)
mirCounts=subset(mirCounts,select=common_sample)
mirCounts=log2(mirCounts+1)

cor=rbind(rnaCounts,mirCounts)

new_row=unique(c("ENSG00000215417",rownames(cor)))
cor=cor[new_row,]
cor=as.data.frame(t(cor))

if(T){
  #杜老师任务
  #"hsa-miR-19a-5p", 几乎不表达
  cor_20220227=subset(cor,select=c("ENSG00000215417","hsa-miR-17-5p","hsa-miR-17-3p","hsa-miR-18a-5p","hsa-miR-18a-3p","hsa-miR-19a-3p","hsa-miR-20a-5p","hsa-miR-20a-3p","hsa-miR-19b-1-5p","hsa-miR-92a-1-5p"))
  write.csv(cor_20220227,"cor_20220227.csv")
  write.csv(box,"cor_20220227_long.csv")
  pdf(file =paste0("Figure A. mir17hg vs others left.pdf"),height = 30,width = 30)
  box=tibble(cor_20220227)
  library(tidyr)
  library(ggpubr)
  box=tidyr::pivot_longer(
    data=box,
    cols=!ENSG00000215417,
    names_to = c("SYMBOL"),
    names_prefix = NULL,
    # names_sep = "_",
    names_pattern = NULL,
    # names_ptypes = list(),
    # names_transform = list(),
    names_repair = "check_unique",
    values_to = "log2(value+1)",
    values_drop_na = T,
    # values_ptypes = list(),
    # values_transform = list()
  )
  p1=ggplot(box, aes(x = ENSG00000215417, y = `log2(value+1)`)) + 
    ylab("")+xlab("")+
    geom_point(shape = 21, colour = "#4682B4", fill = "#87CEFA", size = 10, stroke = .5,alpha=0.8)+ geom_smooth(method="glm",formula = y ~ x,linetype=2,color="#6495ED",fill="#D3D3D3") + theme_bw()+stat_cor(method = 'spearman', aes(),na.rm = T,size = 20)+facet_wrap(vars(box$SYMBOL),ncol = 3,scales ="free_y",strip.position = "left")

  # p2=ggExtra::ggMarginal(p1, type = "density", xparams = list(fill = "#FFE4B5"),yparams = list(fill = "#90EE90"))
  p1
  dev.off()
  
}
library(corrplot)
library(tidyverse)
library(ggplot2)
# cor=as.numercoric(cor)
# cor=subset(cor,select=c("ENSG00000140718",rownames(TCGA_immune)))
# library(tidyverse)
# library(corrplot)
# # p.mat=res1$p,
# 
# # cor=t(cor)
# cor=as.data.frame(cor)
library(doParallel)
cl <- makeCluster(30)
registerDoParallel(cl)
cor_mir17hg <- list()
cor_mir17hg=foreach(i=1:ncol(cor))%dopar%{
  wula=cor.test(cor$ENSG00000215417,cor[[i]],conf.level = .95, method = "pearson")
  wula1=wula[["estimate"]][["cor"]]
  wula2=wula[["p.value"]]
  return(list(wula1,wula2))
}
cor_mir17hg_wuhua=data.table::rbindlist(cor_mir17hg)
cor_mir17hg_wuhua=as.data.frame(cor_mir17hg_wuhua)
rownames(cor_mir17hg_wuhua)=colnames(cor)
cor_mir17hg_wuhua$V2=as.numeric(cor_mir17hg_wuhua$V2)
cor_mir17hg_wuhua=subset(cor_mir17hg_wuhua,V2<0.05)

source(r"(C:\Users\zxh\Desktop\R\20201025 pancreatic cancer validation\updateName.R)")# ncbi与biomart相比多map了200，等于是，ncbi》biomart》org.hs
forthetree=ensg2SymbolUsingNCBImodify(alias=rownames(cor_mir17hg_wuhua),gene.info.file = "Homo_sapiens.gene_info.gz")
cor_mir17hg_wuhua=cbind(cor_mir17hg_wuhua,forthetree)

cor_mir17hg_wuhua=subset(cor_mir17hg_wuhua,gene_type%in%c("protein_coding","lncRNA"))
write.csv(cor_mir17hg_wuhua,file = "mir17hg cor anno.csv")
cor_mir17hg=foreach(i=1:ncol(cor))%dopar%{
  cor_mir17hg[[i]][[2]]
}
stopCluster(cl)

food=cor[,colnames(cor)%in%rownames(cor_mir17hg_wuhua)][,1:10]
res1 <- cor.mtest(cor,conf.level = .95, method = "pearson")
pdf("CoR mir31hg.pdf",width = 10,height=10)
pp1<- cor(food, method ="pearson") %>%  corrplot(type="lower",method = "square",
                                                insig="label_sig",p.mat = res1$p,outline = "white",
                                                sig.level=c(.05,.01,.001),  ,#
                                                pch.cex=.9,pch.col="green",tl.srt=45,tl.col="black")
dev.off()
# order = "AOE", cl.pos = "b", tl.pos = "d", tl.srt = 60) 

pp1=pp1%>%  corrplot(type="lower",method = "square",
                  insig="label_sig",p.mat = res1$p,outline = "white",
                  sig.level=c(.05,.01,.001),  ,#
                  pch.cex=.9,pch.col="green",tl.srt=45,tl.col="black")
# pdf("Construction_risk_score.pdf",width = 8,height=8)
# cor(cor) %>%   corrplot(method = "square",type = "lower",tl.srt = 45,tl.col = "black",)
# dev.off()
pp1

pp2=res1[["p"]]
pp2=as.data.frame(pp2)
pp1=as.data.frame(pp1)
pp1
