#20210118 code wyz split 50%50% high low dePC
#一定0:cpm这里能顺便filter,是1先Q=TMM，2可以或无需CPM(如果是选了TMM,就不能用cpm了，cpm的输入数据是Q后的),3再log或voom，4最后combat（cbcbSEQ）
#一定是0可以或无需单用CPM来filter,1先Q=TMM，2再log或voom，3最后combat（hpgltools？）待证实，以第一行为准。
#' This just calls normalize expt with the most common arguments except log2
#' transformation, but that may be appended with 'transform=log2', so I don't
#' feel bad.  Indeed, it will allow you to overwrite any arguments if you wish.
#' In our work, the most common normalization is: quantile(cpm(low-filter(data))).
# Step 1: not doing count filtering.(cpm法)
# Step 2: normalizing the data with quant.
# Step 3: not converting the data.(convert=Conversion to perform? (raw, cpm, rpkm, cp_seq_m))
# Step 4: transforming the data with log2.
# transform_counts: Found 81436 values equal to 0, adding 1 to the matrix.
# Step 5: not doing batch correction.
library(GDCRNATools)
rm(list=ls())
options(stringsAsFactors = F)
gc()
#最重要的ajust是TSS67,center2627
load("TcgaTargetGTEX_phenotype_stomach.Rdata")
load("STAD_GDCRNATOOLS.Rdata")
load("rnaCountslog21_stomach.Rdata")
load("rnatpmlog2001_stomach.Rdata")
load("rnainter_homo_sig_cerna_all.Rdata")
# load("metaMatrix.RNA.ESCA.Rdata") #没有
#TCGA-SW-A7EB-01 有信息但是没有count数据
# metaMatrix.RNA.STAD=metaMatrix.RNA
# metaMatrix.RNA=rbind(metaMatrix.RNA.STAD,metaMatrix.RNA.ESCA)#写错了
#数据处理RSEM expected count可以当做raw count处理，并且被TMM+voom，
#RSEM 但是找DE时，ENseq>DEseq2>≈edgeR
TcgaTargetGTEX_phenotype_stomach624=TcgaTargetGTEX_phenotype_stomach[-which(rownames(TcgaTargetGTEX_phenotype_stomach) %in% c("TCGA-SW-A7EB-01")),]
table(TcgaTargetGTEX_phenotype_stomach$sample_type)
rownames(rnatpmlog2001_stomach)=substr(rownames(rnatpmlog2001_stomach),1,15)
rownames(rnaCountslog21_stomach)=substr(rownames(rnaCountslog21_stomach),1,15)
rnatpm_stomach=2^rnatpmlog2001_stomach-0.001#625
rnaCounts_stomach=2^rnaCountslog21_stomach-1#624
# save(rnatpm_stomach,file="1.Rdata")
# # 如果输入DESeq2,必须取整,函数round
# rnaCounts_stomach <-round(rnaCounts_stomach,digits = 0)
# DEGAll <- gdcDEAnalysis(counts     = rnaCounts_stomach, 
#                         group      = TcgaTargetGTEX_phenotype_stomach624$sample_type, 
#                         comparison = 'PrimaryTumor-SolidTissueNormal', 
#                         method     = 'DESeq2',
#                         n.cores    = 8,
#                         filter =F)#默认是true，不筛选则选择false

# Normalization of RNAseq data


# rnaExpr<- gdcVoomNormalization(counts = rnaCounts_stomach, filter = F)#limma filter
#我这里改了啊
rnaCounts_stomach_rownames <- rownames(rnaCounts_stomach)
rnaCounts_stomach_colnames <- colnames(rnaCounts_stomach)
rnaCounts_quant <- preprocessCore::normalize.quantiles(
  as.matrix(rnaCounts_stomach))
rownames(rnaCounts_quant) <- rnaCounts_stomach_rownames
colnames(rnaCounts_quant) <- rnaCounts_stomach_colnames
rnaExpr_quant=log2(rnaCounts_quant+1)
mod = model.matrix(~as.factor(sample_type), data=TcgaTargetGTEX_phenotype_stomach624)
# ,mod = NULL,
rnaExpr_quant=ComBat(dat=rnaExpr_quant,batch=TcgaTargetGTEX_phenotype_stomach624$`_study`,mod=mod)
rnaExpr=as.data.frame(rnaExpr_quant)
# pdf("a.pdf")
# boxplot(rnaExpr)
# dev.off()
# Normalization of miRNAs data


# mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = F)#limma filter

mirCounts_rownames <- rownames(mirCounts)
mirCounts_colnames <- colnames(mirCounts)
mirCounts_quant <- preprocessCore::normalize.quantiles(
  as.matrix(mirCounts))
rownames(mirCounts_quant) <- mirCounts_rownames
colnames(mirCounts_quant) <- mirCounts_colnames
mirCounts_quant=as.data.frame(mirCounts_quant)
mirExpr=log2(mirCounts_quant+1)
DEGAll <- gdcDEAnalysis(counts     = rnaCounts_stomach, 
                        group      = TcgaTargetGTEX_phenotype_stomach624$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma',#edgeR
                        filter =F)#默认是true，不筛选则选择false
LNC_0 <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = 0, pval = 100)#
PC_0 <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = 0, pval = 100)#其实放到0.05000也可以
# WGCNA
# rnatpm_stomach_lnc=rnatpm_stomach[rownames(LNC_0),]
# rnatpm_stomach_mrna=rnatpm_stomach[rownames(PC_0),]
# save(rnatpm_stomach_lnc,rnatpm_stomach_mrna,file="rnatpm_stomach_WGCNA.Rdata")
#
DEGAll <- gdcDEAnalysis(counts     = rnaCounts_stomach, 
                        group      = TcgaTargetGTEX_phenotype_stomach624$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma',#edgeR
                        filter =T)#默认是true，不筛选则选择false
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all', fc = 1.0000, pval = 0.05000)#
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = 1.0000, pval = 0.05000)#
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = 1.0000, pval = 0.05000)#其实放到0.05000也可以
LNC_0_F <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = 0, pval = 100)#
PC_0_F <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = 0, pval = 100)
#wyz pn01
phenotype_stomach1=TcgaTargetGTEX_phenotype_stomach624[TcgaTargetGTEX_phenotype_stomach624$sample_type=="PrimaryTumor",]
phenotype_stomach2=readr::read_tsv(file="C:/Users/zxh/Desktop/R/wyz_gc/stad-tcga-gtex/TCGA-STAD.GDC_phenotype.tsv.gz",guess_max=min(1000000, Inf))
phenotype_stomach2=as.data.frame(phenotype_stomach2)
rownames(phenotype_stomach2)=substr(phenotype_stomach2$submitter_id.samples,1,15)
phenotype_stomach2=phenotype_stomach2[rownames(phenotype_stomach2)%in%rownames(phenotype_stomach1),]
table(is.na(phenotype_stomach2$number_of_lymphnodes_positive_by_he))
table(is.na(phenotype_stomach2$pathologic_N))
table(phenotype_stomach2$pathologic_N)
phenotype_stomach2=phenotype_stomach2[!is.na(phenotype_stomach2$pathologic_N),]
phenotype_stomach2=phenotype_stomach2[!phenotype_stomach2$pathologic_N=="NX",]

phenotype_stomach2$pn01=ifelse(phenotype_stomach2$pathologic_N=="N0","N","P")

rnaCounts_T=rnaCounts_stomach[,colnames(rnaCounts_stomach)%in%rownames(phenotype_stomach2)]
phenotype_stomach2=phenotype_stomach2[rownames(phenotype_stomach2)%in%colnames(rnaCounts_T),]
DEN01_tcga <- gdcDEAnalysis(counts     = rnaCounts_T, 
                                  group      = phenotype_stomach2$pn01, 
                                  comparison = 'P-N', 
                                  method     = 'limma',#edgeR
                                  filter =T)#默认是true，不筛选则选择false
save(DEN01_tcga,phenotype_stomach2,rnaCounts_T,file = "DEN01_tcga.Rdata")



######################################这里是传新的按照某个基因split
TcgaTargetGTEX_phenotype_stomach624$coln=rownames(TcgaTargetGTEX_phenotype_stomach624)
if(T){

  TumorOnlysplitinto5050 <- function(datasets,selectgroup,phenoDat,coln,traitV,trait){
    new_phenoDat=subset(phenoDat,phenoDat[[traitV]]==trait)
    new_datasets=as.data.frame(datasets)
    new_datasets=subset(new_datasets,select=colnames(new_datasets) %in% new_phenoDat[[coln]])
    try(if(anyNA(new_datasets[selectgroup,])) stop("Please impute NA value first!"))
    # new_datasets$grouphl=ifelse(new_datasets[[selectgroup]]>quantile(new_datasets[[selectgroup]])[3],1,
    #                             ifelse(new_datasets[[selectgroup]]<=quantile(new_datasets[[selectgroup]])[3],0,NA))
    # df[,-which(names(df)%in%c("a","b")]
    new_datasetsA=subset(new_datasets,select=new_datasets[selectgroup,]>(quantile(new_datasets[selectgroup,])[[3]]))
    new_phenoDatA=new_phenoDat[match(colnames(new_datasetsA),new_phenoDat[[coln]]),]
    new_datasetsB=subset(new_datasets,select=new_datasets[selectgroup,]<=(quantile(new_datasets[selectgroup,])[[3]]))
    new_phenoDatB=new_phenoDat[match(colnames(new_datasetsB),new_phenoDat[[coln]]),]
    new_datasets=cbind(new_datasetsA,new_datasetsB)
    new_phenoDat=rbind(new_phenoDatA,new_phenoDatB)
    new_phenoDat$grouphl=c(rep("H",nrow(new_phenoDatA)),rep("L",nrow(new_phenoDatB)))
    return(list(list(new_phenoDat,new_datasets),list(new_phenoDatA,new_datasetsA),list(new_phenoDatB,new_datasetsB)))
  }
  tcga_5050=TumorOnlysplitinto5050(rnaCounts_stomach,"ENSG00000113504",TcgaTargetGTEX_phenotype_stomach624,"coln","sample_type","PrimaryTumor")
  DEGAll_tcga_5050 <- gdcDEAnalysis(counts     = tcga_5050[[1]][[2]], 
                          group      = tcga_5050[[1]][[1]]$grouphl, 
                          comparison = 'H-L', 
                          method     = 'limma',#edgeR
                          filter =T)#默认是true，不筛选则选择false
  }
save(dePC,PC_0,DEGAll,tcga5050,DEGAll_tcga_5050,file = "tcga.Rdata")
write.csv(dePC,quote = T,file = "dePC_filtered_tcga.csv")
write.csv(PC_0,quote = T,file = "all_gene_weknow.csv")
write.csv(DEGAll,quote = T,file = "DEGAll_filtered_tcga.csv")
write.csv(DEGAll_tcga_5050,quote = T,file = "DEGAll_tcga_filtered_5050.csv")
# write.csv(tcga_5050[[1]][[2]],quote = T,file = "expr_tcga_5050.csv")
# write.csv(tcga_5050[[1]][[1]],quote = T,file = "ph_tcga_5050.csv")
#合并别人的结果 但事实上geo没有filter
load("C:/Users/zxh/Desktop/R/gastric cancer validation/DEGAll_GSE66229.Rdata")
deALL_GSE66229 <- gdcDEReport(deg = DEGAll_GSE66229, gene.type = 'all', fc = 1.0000, pval = 0.05000)#
dePC=dePC[which(dePC$symbol %in% rownames(deALL_GSE66229)),]
deLNC=deLNC[which(deLNC$symbol %in% rownames(deALL_GSE66229)),]
DEGAll=DEGAll[which(DEGAll$symbol %in% rownames(DEGAll_GSE66229)),]
LNC_0_F_over <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = 0, pval = 100)#
PC_0_F_over <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = 0, pval = 100)
# load("C:/Users/zxh/Desktop/R/gastric cancer validation/DEGAll_GSE29272.Rdata")
# deALL_GSE29272 <- gdcDEReport(deg = DEGAll_GSE29272, gene.type = 'all', fc = 1.0000, pval = 0.05)#
# dePC=dePC[which(dePC$symbol %in% rownames(deALL_GSE29272)),]
# deLNC=deLNC[which(deLNC$symbol %in% rownames(deALL_GSE29272)),]
# load("C:/Users/zxh/Desktop/R/gastric cancer validation/DEGAll_GSE26899.Rdata")
# deALL_GSE26899 <- gdcDEReport(deg = DEGAll_GSE26899, gene.type = 'all', fc = 1.0000, pval = 0.05)#
# dePC=dePC[which(dePC$symbol %in% rownames(deALL_GSE26899)),]
# deLNC=deLNC[which(deLNC$symbol %in% rownames(deALL_GSE26899)),]

MIRAll <- gdcDEAnalysis(counts     = mirCounts, 
                        group      = metaMatrix.MIR$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma',#edgeR
                        filter =T)#默认是true，不筛选则选择false
MIR_0 <- gdcDEReport(deg = MIRAll, gene.type = 'miRNAs', fc = 0, pval = 100)#
deMIR <- gdcDEReport(deg = MIRAll, gene.type = 'miRNAs', fc = 1.0000, pval = 0.05000)#
miall <- gdcDEAnalysis(counts     = mirCounts, 
                       group      = metaMatrix.MIR$sample_type, 
                       comparison = 'PrimaryTumor-SolidTissueNormal', 
                       method     = 'limma',#edgeR
                       filter =F)#默认是true，不筛选则选择false
MIRNAallremian <- gdcDEReport(deg = miall, gene.type = 'miRNAs', fc = 0, pval = 100)#

# survOutput <- gdcSurvivalAnalysis(gene     = mrna_2list,
#                                   method   = 'KM',
#                                   rna.expr = rnaExpr,
#                                   metadata = phe,
#                                   sep      = 'median')
# survOutput$pValue=as.character(survOutput$pValue)
# XXXXXXXXXXXXXXXXx=rownames(survOutput[survOutput$pValue<0.05,])
# 
# survOutput <- gdcSurvivalAnalysis(gene     = rownames(dePC), 
#                                   method   = 'coxph', 
#                                   rna.expr = rnaExpr, 
#                                   metadata = metaMatrix.RNA)
# XXXXXXXXXXXXXXXXx=rownames(survOutput[survOutput$pValue<0.05,])
# x2=dePC[XXXXXXXXXXXXXXXXx,]

Gene_up=dePC[dePC$logFC>0,]
Gene_down=dePC[dePC$logFC<0,]
MIR_up=deMIR[deMIR$logFC>0,]
MIR_down=deMIR[deMIR$logFC<0,]
LNC_up=deLNC[deLNC$logFC>0,]
LNC_down=deLNC[deLNC$logFC<0,]

#deMIR只标记#a
ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC),
                          deMIR       = rownames(deMIR),
                          lnc.targets = lnc_inter_target,
                          pc.targets  = mrna_inter_target,
                          rna.expr    = rnaExpr,
                          mir.expr    = mirExpr)#超几何分布及p值

#取出高中文表达的mirna，而不是找差异的miRNA，至少这篇不是#b
miRNA_1=mirCounts[order(apply(mirCounts,1,mad), decreasing = T)[500:0],]


######！！！！！！注意这种cerna只在tumor样本中计算

ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC),
                          deMIR       = rownames(miRNA_1),
                          lnc.targets = lnc_inter_target,
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)#超几何分布及p值

# Filter potential ceRNA interactions
ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.05 & 
                        ceOutput$corPValue<0.05 & ceOutput$regSim != 0,]


# group_list=ifelse(substr(colnames(rnaExpr),1,4)=="GTEX",'normal',
#                      ifelse(substr(colnames(rnaExpr),1,4)=="TCGA" & substr(colnames(rnaExpr),14,15)=="11","normal",
#                             ifelse(substr(colnames(rnaExpr),1,4)=="TCGA" & substr(colnames(rnaExpr),14,15)=="01","tumor","wrong")))
# table(group_list)
# exprSet=rnaExpr[,group_list=='tumor']
#算cerna的时候是带normal的，不应去掉正常。否者会少。
if(F){
  ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC),#unique(ceOutput2$)
                            pc          = rownames(dePC),#unique(ceOutput2)
                            deMIR       = rownames(deMIR),
                            lnc.targets = 'starBase',
                            pc.targets  = 'starBase',
                            rna.expr    = rnaExpr, 
                            mir.expr    = mirExpr)#超几何分布及p值
  # Filter potential ceRNA interactions
  ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.05 & 
                          ceOutput$corPValue<0.05 & ceOutput$regSim != 0,]
  edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
  nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')
  nodes_pc=nodes[nodes$type=='pc',]
  nodes_lnc=nodes[nodes$type=='lnc',]
  lnc_star_target_gdcrnatools=merge(nodes_lnc,edges,by.x="gene",by.y="fromNode")
  mrna_star_target_gdcrnatools=merge(nodes_pc,edges,by.x="gene",by.y="fromNode")
  save(lnc_star_target_gdcrnatools,mrna_star_target_gdcrnatools,file = "star_target_gdcrnatool.Rdata")
  ####造交集
  lnc_star_target=inner_join(lnc_star_target_gdcrnatools,rnainter_homo_sig_cernalnc,by=c("toNode" = "V2","gene" = "gene_id"))
  mrna_star_target=inner_join(mrna_star_target_gdcrnatools,rnainter_homo_sig_cernamrna,by=c("toNode" = "V2","gene" = "gene_id"))
  lnc_target=split(lnc_star_target_gdcrnatools[,5], lnc_star_target_gdcrnatools$gene)
  mrna_target=split(mrna_star_target_gdcrnatools[,5], mrna_star_target_gdcrnatools$gene)
}





#14000 857 buxing
wgcna_mrna_1=read.table(file = "./14000log/CytoscapeInput-nodes-turquoise.txt",sep="\t",stringsAsFactors = F,
                        fill = TRUE,encoding = "UTF-8",header=T)
wgcna_mrna2=read.table(file = "./14000log/red/CytoscapeInput-nodes-red.txt",sep="\t",stringsAsFactors = F,
                       fill = TRUE,encoding = "UTF-8",header=T)
# wgcna_mrna3=read.table(file = "./14000log/green/CytoscapeInput-nodes-green.txt",sep="\t",stringsAsFactors = F,
#                        fill = TRUE,encoding = "UTF-8",header=T)
wgcna_mrna_14000=rbind(wgcna_mrna1,wgcna_mrna2)#,wgcna_mrna3
#13000 794
wgcna_mrna1=read.table(file = "./13000log/CytoscapeInput-nodes-turquoise.txt",sep="\t",stringsAsFactors = F,
                       fill = TRUE,encoding = "UTF-8",header=T)
wgcna_mrna2=read.table(file = "./13000log/green/CytoscapeInput-nodes-green.txt",sep="\t",stringsAsFactors = F,
                       fill = TRUE,encoding = "UTF-8",header=T)
# wgcna_mrna3=read.table(file = "./13000log/black/CytoscapeInput-nodes-black.txt",sep="\t",stringsAsFactors = F,
#                        fill = TRUE,encoding = "UTF-8",header=T)
# wgcna_mrna4=read.table(file = "./13000log/red/CytoscapeInput-nodes-red.txt",sep="\t",stringsAsFactors = F,
#                        fill = TRUE,encoding = "UTF-8",header=T)
wgcna_mrna_13000=rbind(wgcna_mrna1,wgcna_mrna2)#,wgcna_mrna3,wgcna_mrna4
#12000
wgcna_mrna1=read.table(file = "./12000logsft080/black/CytoscapeInput-nodes-black.txt",sep="\t",stringsAsFactors = F,
                       fill = TRUE,encoding = "UTF-8",header=T)
wgcna_mrna2=read.table(file = "./12000logsft080/CytoscapeInput-nodes-turquoise.txt",sep="\t",stringsAsFactors = F,
                       fill = TRUE,encoding = "UTF-8",header=T)
wgcna_mrna3=read.table(file = "./12000logsft080/green/CytoscapeInput-nodes-green.txt",sep="\t",stringsAsFactors = F,
                       fill = TRUE,encoding = "UTF-8",header=T)
wgcna_mrna4=read.table(file = "./12000logsft080/yellow/CytoscapeInput-nodes-yellow.txt",sep="\t",stringsAsFactors = F,
                       fill = TRUE,encoding = "UTF-8",header=T)
wgcna_mrna_12000=rbind(wgcna_mrna1,wgcna_mrna2,wgcna_mrna3,wgcna_mrna4)
#12000 NON-LOG 500+
wgcna_mrna_12000NONLOG=read.table(file = "./70rna12000/turquoise/CytoscapeInput-nodes-turquoise.txt",sep="\t",stringsAsFactors = F,
                                  fill = TRUE,encoding = "UTF-8",header=T)
#16000
wgcna_mrna_16000NONLOG1=read.table(file = "./16000/CytoscapeInput-nodes-turquoise.txt",sep="\t",stringsAsFactors = F,
                                   fill = TRUE,encoding = "UTF-8",header=T)
# wgcna_mrna_16000NONLOG2=read.table(file = "./16000/16000black/CytoscapeInput-nodes-black.txt",sep="\t",stringsAsFactors = F,
#                                    fill = TRUE,encoding = "UTF-8",header=T)
# wgcna_mrna_16000NONLOG3=read.table(file = "./16000/grey/CytoscapeInput-nodes-grey.txt",sep="\t",stringsAsFactors = F,
#                                    fill = TRUE,encoding = "UTF-8",header=T)
wgcna_mrna_16000NONLOG=rbind(wgcna_mrna_16000NONLOG1)#,wgcna_mrna_16000NONLOG2,wgcna_mrna_16000NONLOG3
#
if(T){
  mrna_1=list(Gene_name=rownames(PC_0_F_over),
              Gene_up=rownames(Gene_up),
              Gene_down=rownames(Gene_down),
              # Gene_model=unique(ceOutput5$Genes),
              Gene_predict=unique(ceOutput2$Genes),
              Gene_WGCNA=wgcna_mrna$nodeName,
              Survival=rownames(x2))
  
  pdf(file='upset_mrna_1.pdf',height = 6,width = 8,onefile = T)
  upset(fromList(mrna_1), #输入的数据集
        sets = c("Gene_WGCNA","Gene_down","Gene_up","Gene_predict"),#,"Survival"
        # sets = c("Gene_WGCNA","Gene_up","Gene_down","Gene_predict"),
        # sets = c("PITA", "miRanda","TargetScan"),
        nsets = 4, #想要可视化的数据集数量,也可以用sets选项自定义
        nintersects = 20, #要绘制的交点数量 2^n
        keep.order = T, #是否保持输入顺序,否则按各集合大小排序
        main.bar.color = 'black', #主条图的颜色
        mainbar.y.label = 'Number', #y轴标签
        sets.bar.color = 'blue', #设置集合条图的颜色
        sets.x.label = 'Number of elements', #设置集合调图的坐标轴标签
        mb.ratio = c(0.7,0.3), #条形图点图的比例
        order.by = c("freq","degree"), #交集如何排序,
        decreasing = c(TRUE,TRUE),
        # 以上排序是否降序,FALSE
        # boxplot.summary = c('ReleaseDate','Comedy') #添加箱线图展示数据分布,最多展示两个数据集
        queries = list(list(query = intersects,
                            params = list("Gene_up","Gene_predict", "Gene_WGCNA"),#,"Survival"
                            color = 'red',
                            active = T)
                       ,
                       list(query = intersects,
                            params = list("Gene_down","Gene_predict", "Gene_WGCNA"),#,"Survival"
                            color = 'green',
                            active = T)
        )
  )
  dev.off()
}
wgcna_mrna=wgcna_mrna_13000
mrna_1up=list(Gene_up=rownames(Gene_up),
              Gene_predict=unique(unique(ceOutput2$Genes)),
              Gene_WGCNA=wgcna_mrna$nodeName)

mrna_1down=list(Gene_down=rownames(Gene_down),
                Gene_predict=unique(unique(ceOutput2$Genes)),
                Gene_WGCNA=wgcna_mrna$nodeName)

gene_uplist=Reduce(intersect,mrna_1up)
gene_downlist=Reduce(intersect,mrna_1down)

mrna_2list=c(gene_uplist,gene_downlist)

#14000 361 297(2)
#13000 417 387(2)
#12000 422 385(2)
#16000nonlog 384 370(1)
#12000NONLOG 200
# mod = model.matrix(~as.factor(sample_type), data=TcgaTargetGTEX_phenotype_stomach624)
# # ,mod = NULL,
# rnaExpr_quant=ComBat(dat=rnaExpr_quant,batch=TcgaTargetGTEX_phenotype_stomach624$`_study`,mod=mod)
# rnaExpr=as.data.frame(rnaExpr_quant)
survOutput1 <- gdcSurvivalAnalysis(gene     = mrna_2list, 
                                   method   = 'KM', 
                                   rna.expr = rnaExpr, 
                                   metadata = metaMatrix.RNA, 
                                   sep      = 'median')
survOutput1$pValue=as.character(survOutput1$pValue)
mrna_2listcox=rownames(survOutput1[survOutput1$pValue<0.05,])
# CoxPH analysis
survOutput2 <- gdcSurvivalAnalysis(gene     = mrna_2listcox, 
                                   method   = 'coxph', 
                                   rna.expr = rnaExpr, 
                                   metadata = metaMatrix.RNA)
mrna_3listcox=rownames(survOutput2[survOutput2$pValue<0.05,])

m3=dePC[mrna_3listcox,]
n4=deLNC[unique(ceOutput5$lncRNAs),]


ceOutput3=ceOutput2[ceOutput2$deMIRCounts>0,]
ceOutput3$Counts=ceOutput3$deMIRCounts
ceOutput3$miRNAs=ceOutput3$deMIRs#用de造出来的
ceOutput4=ceOutput3[ceOutput3$Genes %in% mrna_3listcox, ]#交mrna
ceOutput5=ceOutput4[ceOutput4$lncRNAs %in% wgcna_lnc$nodeName,]#交lncrna
#这里顺便弄回来了
ceOutput5=llllll3
sum(as.numeric(ceOutput5$Counts))
edges <- gdcExportNetwork(ceNetwork = ceOutput5, net = 'edges')
nodes <- gdcExportNetwork(ceNetwork = ceOutput5, net = 'nodes')

write.csv(edges, file='edges.csv',quote = F) ### Network of Cytoscape
write.csv(nodes, file='nodes.csv',quote = F) ### Table of Cytoscape

ceOutput_readable=merge(ceOutput5,nodes,by.x ="lncRNAs",by.y= "gene", all.x=TRUE,sort=TRUE) 
ceOutput_readable$lncRNAs=ceOutput_readable$symbol
ceOutput_readable=merge(ceOutput_readable,nodes,by.x ="Genes",by.y= "gene", all.x=TRUE,sort=TRUE) 
ceOutput_readable$Genes=ceOutput_readable$symbol.y
paste(colnames(ceOutput_readable),collapse = "','")
ceOutput_readable=ceOutput_readable[,c('Genes','lncRNAs','miRNAs','hyperPValue','cor','corPValue')]
write.csv(ceOutput_readable, file='ceOutput5.csv',quote = T)
paste(nodes$gene,collapse = ",")
paste(nodes$symbol,collapse = ",")


wgcna_lnc=read.table(file = "./2000lnc10best/CytoscapeInput-nodes-yellow.txt",sep="\t",stringsAsFactors = F,
                     fill = TRUE,encoding = "UTF-8",header=T)

# lncrna_1=list(lncRNA_up=rownames(LNC_up),
#               lncRNA_down=rownames(LNC_down),
#               lncRNA_predict=unique(ceOutput4$lncRNAs),
#               lncRNA_WGCNA=wgcna_lnc$nodeName)

lncrna_1=list(lncRNA_up=rownames(LNC_up),
              lncRNA_down=rownames(LNC_down),
              lncRNA_predict=unique(ceOutput5$lncRNAs))

pdf(file='upset_lncrna_1.pdf',height = 6,width = 8, onefile = F)
upset(fromList(lncrna_1), #输入的数据集
      sets = c("lncRNA_down","lncRNA_up","lncRNA_predict"),#"lncRNA_WGCNA",
      nsets = 5, #想要可视化的数据集数量,也可以用sets选项自定义
      nintersects = 12, #要绘制的交点数量 2^n
      keep.order = T, #是否保持输入顺序,否则按各集合大小排序
      main.bar.color = 'black', #主条图的颜色
      mainbar.y.label = 'Number', #y轴标签
      sets.bar.color = 'blue', #设置集合条图的颜色
      sets.x.label = 'Number of elements', #设置集合调图的坐标轴标签
      mb.ratio = c(0.7,0.3), #条形图点图的比例
      order.by = c("freq","degree"), #交集如何排序,
      decreasing = c(TRUE,TRUE),
      # 以上排序是否降序,FALSE
      # boxplot.summary = c('ReleaseDate','Comedy') #添加箱线图展示数据分布,最多展示两个数据集
      queries = list(list(query = intersects,
                          params = list("lncRNA_up","lncRNA_predict"),#, "lncRNA_WGCNA"
                          color = 'red',
                          active = T)
                     ,
                     list(query = intersects,
                          params = list("lncRNA_down","lncRNA_predict"),#"Gene_WGCNA"
                          color = 'green',
                          active = T)
      )
)
dev.off()

y2=deLNC[unique(ceOutput5$lncRNAs),]
paste(unique(ceOutput_readable$lncRNAs),collapse = ", ")

shinyKMPlot(gene=rownames(dePC), rna.expr=rnaExpr, 
            metadata=metaMatrix.RNA)

##################################################################################################
不是所有的都能筛选出好的risk图的

r3<-rbind(r1,r2)
r3_dePC=dePC
r3_dePC$nodeName=rownames(r3_dePC)
# r3_dePC=dePC[rownames(r3),]
r3_dePC=merge(r3_dePC,r3,all=F,sort=TRUE,by="nodeName") 
rownames(r3_dePC)=r3_dePC$nodeName
write.csv(r3_dePC, file='r3_dePC.csv',quote = T)

# survOutput1 <- gdcSurvivalAnalysis(gene     = mrna_2list, 
#                                    method   = 'coxph', 
#                                    rna.expr = rnaExpr, 
#                                    metadata = metaMatrix.RNA, 
#                                    sep      = 'median')
# mrna_2listcox=rownames(survOutput1[survOutput1$pValue<0.05,])
# 
# 
# 
cox_results <-apply(exprSet_input , 1 , function(gene){
  # gene= exprSet[1,]
  group=ifelse(gene>median(gene),'high','low')
  
})
slnc <- gdcSurvivalAnalysis(gene     = mrna_2list, 
                            method   = 'coxph', 
                            rna.expr = exprSet_input, 
                            metadata = metaMatrix.RNA, 
                            cox_results= cox_results)
cox_results <- data.frame(group=cox_results,stringsAsFactors = F)

slnc$pValue=as.character(slnc$pValue)
slnc=rownames(slnc[slnc$pValue<0.05,])
slnc=deLNC[slnc,]
x3=dePC[rownames(exprSet2[c(15,54,59,65,68,71,94,112,140,157),]),]


rev(strsplit('PrimaryTumor-SolidTissueNormal', "-", fixed = TRUE)[[1]])
