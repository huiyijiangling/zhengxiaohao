# genes_expr_raw=2^(genes_expr)-1
# probes_expr_raw=2^(probes_expr)-1
# rnaExpr_v1 <- gdcVoomNormalization(counts = probes_expr_raw, filter = T)
# GSE51575 不用标准化
?gdcDEAnalysis
DEGAll_v1 <- gdcDEAnalysis(counts     = genes_expr, 
                        group      = gset$description, 
                        comparison = 'G0-G1', 
                        method     = 'limma',
                        filter =T)#默认是true，不筛选则选择false
LNC_0 <- gdcDEReport(deg = DEGAll_v1, gene.type = 'long_non_coding', fc = 0, pval = 100)#
PC_0 <- gdcDEReport(deg = DEGAll_v1, gene.type = 'protein_coding', fc = 0, pval = 100)#其实放到0.01也可以

deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all', fc = 1.5, pval = 0.01)#
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = 1.5, pval = 0.01)#
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = 1.5, pval = 0.01)#其实放到0.01也可以

MIRAll <- gdcDEAnalysis(counts     = mirCounts, 
                        group      = metaMatrix.MIR$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma',
                        filter =T)#默认是true，不筛选则选择false
MIR_0 <- gdcDEReport(deg = MIRAll, gene.type = 'miRNAs', fc = 0, pval = 100)#
deMIR <- gdcDEReport(deg = MIRAll, gene.type = 'miRNAs', fc = 1.5, pval = 0.01)#
miall <- gdcDEAnalysis(counts     = mirCounts, 
                       group      = metaMatrix.MIR$sample_type, 
                       comparison = 'PrimaryTumor-SolidTissueNormal', 
                       method     = 'limma',
                       filter =F)#默认是true，不筛选则选择false
MIRNAallremian <- gdcDEReport(deg = miall, gene.type = 'miRNAs', fc = 0, pval = 100)#
rnaall <- gdcDEAnalysis(counts     = rnaCounts, 
                        group      = metaMatrix.RNA$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma',
                        filter =F)#默认是true，不筛选则选择false
Geneall <- gdcDEReport(deg = rnaall, gene.type = 'protein_coding', fc = 0, pval = 100)#
LNCall <- gdcDEReport(deg = rnaall, gene.type = 'long_non_coding', fc = 0, pval = 100)#


survOutput <- gdcSurvivalAnalysis(gene     = rownames(dePC), 
                                  method   = 'KM', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA, 
                                  sep      = 'median')
survOutput$pValue=as.character(survOutput$pValue)
XXXXXXXXXXXXXXXXx=rownames(survOutput[survOutput$pValue<0.05,])

survOutput <- gdcSurvivalAnalysis(gene     = XXXXXXXXXXXXXXXXx, 
                                  method   = 'coxph', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA)
XXXXXXXXXXXXXXXXx=rownames(survOutput[survOutput$pValue<0.05,])
x2=dePC[XXXXXXXXXXXXXXXXx,]

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

ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC),
                          deMIR       = rownames(miRNA_1),
                          lnc.targets = lnc_inter_target,
                          pc.targets  = mrna_inter_target,
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)#超几何分布及p值

# Filter potential ceRNA interactions
ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.01 & 
                        ceOutput$corPValue<0.01 & ceOutput$regSim != 0,]

wgcna_mrna=read.table(file = "./17000rna92/CytoscapeInput-nodes-turquoise.txt",sep="\t",stringsAsFactors = F,
                      fill = TRUE,encoding = "UTF-8",header=T)
# wgcna_mrna=read.table(file = "./9000rna/CytoscapeInput-nodes-turquoise.txt",sep="\t",stringsAsFactors = F,
#                       fill = TRUE,encoding = "UTF-8",header=T)


mrna_1=list(Gene_name=rownames(Geneall),
            Gene_up=rownames(Gene_up),
            Gene_down=rownames(Gene_down),
            # Gene_predict=unique(starbasemrna_filter$geneID),
            Gene_predict=unique(ceOutput3$Genes),
            Gene_WGCNA=wgcna_mrna$nodeName,
            Survival=rownames(x2))

pdf(file='upset_mrna_1.pdf',height = 6,width = 8)
upset(fromList(mrna_1), #输入的数据集
      sets = c("Gene_WGCNA","Gene_down","Gene_up","Gene_predict","Survival"),
      # sets = c("Gene_WGCNA","Gene_up","Gene_down","Gene_predict"),
      # sets = c("PITA", "miRanda","TargetScan"),
      nsets = 5, #想要可视化的数据集数量,也可以用sets选项自定义
      nintersects = 12, #要绘制的交点数量 2^n
      keep.order = T, #是否保持输入顺序,否则按各集合大小排序
      main.bar.color = 'black', #主条图的颜色
      mainbar.y.label = 'Number', #y轴标签
      sets.bar.color = 'black', #设置集合条图的颜色
      sets.x.label = 'Number of elements', #设置集合调图的坐标轴标签
      mb.ratio = c(0.7,0.3), #条形图点图的比例
      order.by = c("freq","degree"), #交集如何排序,
      decreasing = c(TRUE,TRUE),
      # 以上排序是否降序,FALSE
      # boxplot.summary = c('ReleaseDate','Comedy') #添加箱线图展示数据分布,最多展示两个数据集
      queries = list(list(query = intersects,
                          params = list("Gene_up","Gene_predict", "Gene_WGCNA","Survival"),
                          color = 'red',
                          active = T)
                     # ,
                     # list(query = intersects,
                     #      params = list("Gene_down","Gene_predict", "Gene_WGCNA","Survival"),
                     #      color = 'green',
                     #      active = T)
      )
)
dev.off()

mrna_1up=list(Gene_up=rownames(Gene_up),
              Gene_predict=unique(unique(ceOutput3$Genes)),
              Gene_WGCNA=wgcna_mrna$nodeName)

mrna_1down=list(Gene_down=rownames(Gene_down),
                Gene_predict=unique(unique(ceOutput3$Genes)),
                Gene_WGCNA=wgcna_mrna$nodeName)

gene_uplist=Reduce(intersect,mrna_1up)
gene_downlist=Reduce(intersect,mrna_1down)

mrna_2list=c(gene_uplist,gene_downlist)

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

wgcna_lnc=read.table(file = "./2000lnc10best/CytoscapeInput-nodes-yellow.txt",sep="\t",stringsAsFactors = F,
                     fill = TRUE,encoding = "UTF-8",header=T)

lncrna_1=list(lncRNA_up=rownames(LNC_up),
              lncRNA_down=rownames(LNC_down),
              lncRNA_predict=unique(ceOutput4$lncRNAs),
              lncRNA_WGCNA=wgcna_lnc$nodeName)

pdf(file='upset_lncrna_1.pdf',height = 6,width = 8)
upset(fromList(lncrna_1), #输入的数据集
      sets = c("lncRNA_WGCNA","lncRNA_down","lncRNA_up","lncRNA_predict"),
      nsets = 5, #想要可视化的数据集数量,也可以用sets选项自定义
      nintersects = 12, #要绘制的交点数量 2^n
      keep.order = T, #是否保持输入顺序,否则按各集合大小排序
      main.bar.color = 'black', #主条图的颜色
      mainbar.y.label = 'Number', #y轴标签
      sets.bar.color = 'black', #设置集合条图的颜色
      sets.x.label = 'Number of elements', #设置集合调图的坐标轴标签
      mb.ratio = c(0.7,0.3), #条形图点图的比例
      order.by = c("freq","degree"), #交集如何排序,
      decreasing = c(TRUE,TRUE),
      # 以上排序是否降序,FALSE
      # boxplot.summary = c('ReleaseDate','Comedy') #添加箱线图展示数据分布,最多展示两个数据集
      queries = list(list(query = intersects,
                          params = list("lncRNA_up","lncRNA_predict", "lncRNA_WGCNA"),
                          color = 'red',
                          active = T)
                     # ,
                     # list(query = intersects,
                     #      params = list("Gene_down","Gene_predict", "Gene_WGCNA"),
                     #      color = 'green',
                     #      active = T)
      )
)
dev.off()

shinyKMPlot(gene=rownames(dePC), rna.expr=rnaExpr, 
            metadata=metaMatrix.RNA)

##################################################################################################
不是所有的都能筛选出好的risk图的

# survOutput1 <- gdcSurvivalAnalysis(gene     = mrna_2list, 
#                                    method   = 'coxph', 
#                                    rna.expr = rnaExpr, 
#                                    metadata = metaMatrix.RNA, 
#                                    sep      = 'median')
# mrna_2listcox=rownames(survOutput1[survOutput1$pValue<0.05,])
# 
# 
# 
# cox_results <-apply(exprSet , 1 , function(gene){
#    # gene= exprSet[1,]
#    group=ifelse(gene>median(gene),'high','low')
#    
# })
slnc <- gdcSurvivalAnalysis(gene     = rownames(deLNC), 
                            method   = 'KM', 
                            rna.expr = rnaExpr, 
                            metadata = metaMatrix.RNA, 
                            sep      = 'median')
slnc$pValue=as.character(slnc$pValue)
slnc=rownames(slnc[slnc$pValue<0.05,])
slnc=deLNC[slnc,]
