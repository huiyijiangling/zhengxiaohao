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
rm(list=ls())
options(stringsAsFactors = F)
gc()
#最重要的ajust是TSS67,center2627
load("TcgaTargetGTEX_phenotype_Pancreas.Rdata")
load("PAAD_GDCRNATOOLS.Rdata")
load("rnaCountslog21_Pancreas.Rdata")
load("rnafpkmlog2001_Pancreas.Rdata")
load("rnainter_homo_sig_cerna_all.Rdata")
# load("metaMatrix.RNA.ESCA.Rdata") #没有
#TCGA-SW-A7EB-01 有信息但是没有count数据
# metaMatrix.RNA.PAAD=metaMatrix.RNA
# metaMatrix.RNA=rbind(metaMatrix.RNA.PAAD,metaMatrix.RNA.ESCA)#写错了
#数据处理RSEM expected count可以当做raw count处理，并且被TMM+voom，
#RSEM 但是找DE时，ENseq>DEseq2>≈edgeR
#z注意删除 注意质控 在运行完之后回头删除
table(TcgaTargetGTEX_phenotype_Pancreas$sample_type)
rownames(TcgaTargetGTEX_phenotype_Pancreas[TcgaTargetGTEX_phenotype_Pancreas$sample_type=="WRONG",])#"TCGA-HZ-A9TJ-06"
TcgaTargetGTEX_phenotype_Pancreas624=TcgaTargetGTEX_phenotype_Pancreas[-which(rownames(TcgaTargetGTEX_phenotype_Pancreas) %in% c("TCGA-HZ-A9TJ-06")),]
table(TcgaTargetGTEX_phenotype_Pancreas624$sample_type)
rownames(rnafpkmlog2001_Pancreas)=substr(rownames(rnafpkmlog2001_Pancreas),1,15)
rownames(rnaCountslog21_Pancreas)=substr(rownames(rnaCountslog21_Pancreas),1,15)
# 删除正常样本
TcgaTargetGTEX_phenotype_Pancreas624=TcgaTargetGTEX_phenotype_Pancreas624[-which(TcgaTargetGTEX_phenotype_Pancreas624$`_sample_type`=="Solid Tissue Normal"),]
# 删除pca不合格的样本 
TcgaTargetGTEX_phenotype_Pancreas624=TcgaTargetGTEX_phenotype_Pancreas624[-which(rownames(TcgaTargetGTEX_phenotype_Pancreas624) %in% c("GTEX-W5X1-0226-SM-5CHTO","GTEX-U8XE-2026-SM-5CHQF","GTEX-146FH-1826-SM-5QGQ7","GTEX-13NZA-1726-SM-5J1NA")),]#,"TCGA−H6−A45N−11" 
# 删除pca不合格的样本 目视法
TcgaTargetGTEX_phenotype_Pancreas624=TcgaTargetGTEX_phenotype_Pancreas624[-which(rownames(TcgaTargetGTEX_phenotype_Pancreas624) %in% c("TCGA-F2-6880-01")),]#,"TCGA−H6−A45N−11"
#rnafpkmlog2001_Pancreas=subset(rnafpkmlog2001_Pancreas,select=c(-`GTEX-W5X1-0226-SM-5CHTO`,-`GTEX-U8XE-2026-SM-5CHQF`,-`GTEX-146FH-1826-SM-5QGQ7`,-`GTEX-13NZA-1726-SM-5J1NA`))#,-`TCGA-H6-A45N-11`
#rnaCountslog21_Pancreas=subset(rnaCountslog21_Pancreas,select=c(-`GTEX-W5X1-0226-SM-5CHTO`,-`GTEX-U8XE-2026-SM-5CHQF`,-`GTEX-146FH-1826-SM-5QGQ7`,-`GTEX-13NZA-1726-SM-5J1NA`))#,-`TCGA-H6-A45N-11`
rnafpkmlog2001_Pancreas=subset(rnafpkmlog2001_Pancreas,select=rownames(TcgaTargetGTEX_phenotype_Pancreas624))
rnaCountslog21_Pancreas=subset(rnaCountslog21_Pancreas,select=rownames(TcgaTargetGTEX_phenotype_Pancreas624))
rnafpkm_Pancreas=2^rnafpkmlog2001_Pancreas-0.001#625
rnaCounts_Pancreas=2^rnaCountslog21_Pancreas-1#624
# # 如果输入DESeq2,必须取整,函数round
# rnaCounts_Pancreas <-round(rnaCounts_Pancreas,digits = 0)
# DEGAll <- gdcDEAnalysis(counts     = rnaCounts_Pancreas, 
#                         group      = TcgaTargetGTEX_phenotype_Pancreas624$sample_type, 
#                         comparison = 'PrimaryTumor-SolidTissueNormal', 
#                         method     = 'DESeq2',
#                         n.cores    = 8,
#                         filter =F)#默认是true，不筛选则选择false

# Normalization of RNAseq data

# rnaCountslog21_Pancreas$
#    #
# rnaExpr<- gdcVoomNormalization(counts = rnaCounts_Pancreas, filter = F)#limma filter
#我这里改了啊
rnaCounts_Pancreas_rownames <- rownames(rnaCounts_Pancreas)
rnaCounts_Pancreas_colnames <- colnames(rnaCounts_Pancreas)
rnaCounts_quant <- preprocessCore::normalize.quantiles(
  as.matrix(rnaCounts_Pancreas))
rownames(rnaCounts_quant) <- rnaCounts_Pancreas_rownames
colnames(rnaCounts_quant) <- rnaCounts_Pancreas_colnames
rnaExpr_quant=log2(rnaCounts_quant+1)

# ,mod = NULL,这里的含义是协变量，告诉其本身有内在的分组，如果不删除tcga里癌旁正常，可以留着。
# mod = model.matrix(~as.factor(sample_type), data=TcgaTargetGTEX_phenotype_Pancreas624)
# rnaExpr_quant=ComBat(dat=rnaExpr_quant,batch=TcgaTargetGTEX_phenotype_Pancreas624$`_study`,mod=mod)
# 千万别一个不带mod的不然结果很难看的
# rnaExpr_quant=ComBat(dat=rnaExpr_quant,batch=TcgaTargetGTEX_phenotype_Pancreas624$`_study`)
rnaExpr=as.data.frame(rnaExpr_quant)
# pdf("a.pdf")
# boxplot(rnaExpr)
# dev.off()
# Normalization of miRNAs data
save(rnaExpr,TcgaTargetGTEX_phenotype_Pancreas624,rnafpkmlog2001_Pancreas,rnafpkm_Pancreas,file="rnafpkmlog2001_Pancreas_pca.Rdata")

# mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = F)#limma filter

mirCounts_rownames <- rownames(mirCounts)
mirCounts_colnames <- colnames(mirCounts)
mirCounts_quant <- preprocessCore::normalize.quantiles(
  as.matrix(mirCounts))
rownames(mirCounts_quant) <- mirCounts_rownames
colnames(mirCounts_quant) <- mirCounts_colnames
mirCounts_quant=as.data.frame(mirCounts_quant)
mirExpr=log2(mirCounts_quant+1)
DEGAll <- gdcDEAnalysis(counts     = rnaCounts_Pancreas, 
                        group      = TcgaTargetGTEX_phenotype_Pancreas624$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma',#edgeR
                        filter =F)#默认是true，不筛选则选择false
LNC_0 <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = 0, pval = 100)#
PC_0 <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = 0, pval = 100)#其实放到0.01也可以
# WGCNA 在用pca时必须打开
rnafpkm_Pancreas_lnc=rnafpkm_Pancreas[rownames(LNC_0),]
rnafpkm_Pancreas_mrna=rnafpkm_Pancreas[rownames(PC_0),]
save(rnafpkm_Pancreas_lnc,rnafpkm_Pancreas_mrna,file="rnafpkm_Pancreas_WGCNA.Rdata")
#
DEGAll <- gdcDEAnalysis(counts     = rnaCounts_Pancreas, 
                        group      = TcgaTargetGTEX_phenotype_Pancreas624$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma',#edgeR
                        filter =T)#默认是true，不筛选则选择false
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all', fc = 1.3, pval = 0.01)#
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = 1.3, pval = 0.01)#
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = 1.3, pval = 0.01)#其实放到0.01也可以
LNC_0_F <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = 0, pval = 100)#
PC_0_F <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = 0, pval = 100)
#合并别人的结果 但事实上geo没有filter
load("C:/Users/zxh/Desktop/R/gastric cancer validation/DEGAll_GSE66229.Rdata")
deALL_GSE66229 <- gdcDEReport(deg = DEGAll_GSE66229, gene.type = 'all', fc = 1.3, pval = 0.01)#
dePC=dePC[which(dePC$symbol %in% rownames(deALL_GSE66229)),]
deLNC=deLNC[which(deLNC$symbol %in% rownames(deALL_GSE66229)),]
DEGAll=DEGAll[which(DEGAll$symbol %in% rownames(DEGAll_GSE66229)),]
LNC_0_F_over <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = 0, pval = 100)#
PC_0_F_over <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = 0, pval = 100)
# load("C:/Users/zxh/Desktop/R/gastric cancer validation/DEGAll_GSE29272.Rdata")
# deALL_GSE29272 <- gdcDEReport(deg = DEGAll_GSE29272, gene.type = 'all', fc = 1.3, pval = 0.05)#
# dePC=dePC[which(dePC$symbol %in% rownames(deALL_GSE29272)),]
# deLNC=deLNC[which(deLNC$symbol %in% rownames(deALL_GSE29272)),]
# load("C:/Users/zxh/Desktop/R/gastric cancer validation/DEGAll_GSE26899.Rdata")
# deALL_GSE26899 <- gdcDEReport(deg = DEGAll_GSE26899, gene.type = 'all', fc = 1.3, pval = 0.05)#
# dePC=dePC[which(dePC$symbol %in% rownames(deALL_GSE26899)),]
# deLNC=deLNC[which(deLNC$symbol %in% rownames(deALL_GSE26899)),]

MIRAll <- gdcDEAnalysis(counts     = mirCounts, 
                        group      = metaMatrix.MIR$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma',#edgeR
                        filter =T)#默认是true，不筛选则选择false
MIR_0 <- gdcDEReport(deg = MIRAll, gene.type = 'miRNAs', fc = 0, pval = 100)#
deMIR <- gdcDEReport(deg = MIRAll, gene.type = 'miRNAs', fc = 1.3, pval = 0.01)#
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
# mod = model.matrix(~as.factor(sample_type), data=TcgaTargetGTEX_phenotype_Pancreas624)
# # ,mod = NULL,
# rnaExpr_quant=ComBat(dat=rnaExpr_quant,batch=TcgaTargetGTEX_phenotype_Pancreas624$`_study`,mod=mod)
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