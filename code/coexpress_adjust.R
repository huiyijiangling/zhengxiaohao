boxplot(mirExpr_tcga)

mirExpr_tcga
mirExpr_GSE32688
mirExpr_GSE119794
mirExpr_GSE41372
mirExpr_GSE43797

pdf("box.pdf", height=10, width=9) 
par(mfrow=c(3,2))
boxplot(mirExpr_tcga)
boxplot(mirExpr_GSE119794)
boxplot(mirExpr_GSE32688)
boxplot(mirExpr_GSE41372)
boxplot(mirExpr_GSE41372)
dev.off() 

pdf("boxrna.pdf", height=10, width=9) 
par(mfrow=c(3,2))
boxplot(rnaExpr_tcga)
boxplot(rnaExpr_GSE119794)
boxplot(rnaExpr_GSE32688)
boxplot(rnaExpr_GSE41372)
boxplot(rnaExpr_GSE41372)
dev.off() 

pdf("bo.pdf", height=10, width=9) 
par(mfrow=c(3,2))
boxplot(rnaExpr_tcga)
boxplot(rnaCounts_quant)
dev.off() 



mirnaT=list(selectTrait(mirExpr_tcga,metaMatrix.MIR,"sample","sample_type","PrimaryTumor"),
            selectTrait(mirExpr_GSE32688,phenoDat_GSE32688_mir,"newname","disease status:ch1","pancreatic cancer"),
            selectTrait(mirExpr_GSE119794,phenoDat_GSE119794_mir,"newname","source_name_ch1","PC"),
            selectTrait(mirExpr_GSE41372,phenoDat_GSE41372_mir,"newname","tissue:ch1","pancreatic ductal adenocarcinoma"),
            selectTrait(mirExpr_GSE43797,ph
)




library(WGCNA) # (Section will take ~5-10 minutes to run)
# source("collapseRows_NEW.R") # ONLY uncomment this line if you get an error with it commented

# commonProbesA = intersect (rownames(mirExpr_tcga),rownames(mirExpr_GSE119794))
common_mir
common_ensg
datExprA1p = mirExpr_tcga[commonProbesA,]
datExprA2p = mirExpr_GSE32688[commonProbesA,]
datExprA3p = mirExpr_GSE119794[commonProbesA,]


rankExprA1= rank(rowMeans(datExprA1p))
rankExprA2= rank(rowMeans(datExprA2p))
rankExprA3= rank(rowMeans(datExprA3p))



datExprA1p = mirExpr_tcga[common_mir,]
datExprA2p = mirExpr_GSE32688[common_mir,]
datExprA3p = mirExpr_GSE119794[common_mir,]
datExprA4p = mirExpr_GSE41372[common_mir,]
datExprA5p = mirExpr_GSE43797[common_mir,]

rankExprA1= rank(rowMeans(datExprA1p))
rankExprA2= rank(rowMeans(datExprA2p))
rankExprA3= rank(rowMeans(datExprA3p))
rankExprA4= rank(rowMeans(datExprA4p))
rankExprA5= rank(rowMeans(datExprA5p))




pdf("generalNetworkProperties.pdf", height=10, width=9) 

par(mfrow=c(2,2))
verboseScatterplot(rankExprA1,rankExprA2, xlab="Ranked Expression (A1)",
                   ylab="Ranked Expression (A2)")
verboseScatterplot(rankExprA1,rankExprA3, xlab="Ranked Connectivity (A1)",
                   ylab="Ranked Connectivity (A2)")
verboseScatterplot(rankExprA1,rankExprA4, xlab="Ranked Expression (A1)",
                   ylab="Ranked Expression (A2)")
verboseScatterplot(rankExprA1,rankExprA5, xlab="Ranked Connectivity (A1)",
                   ylab="Ranked Connectivity (A2)")
dev.off() 

datExprA1p = rnaExpr_tcga[common_ensg,]
datExprA2p = rnaExpr_GSE32688[common_ensg,]
datExprA3p = rnaExpr_GSE119794[common_ensg,]
datExprA4p = rnaExpr_GSE41372[common_ensg,]
datExprA5p = rnaExpr_GSE43797[common_ensg,]

rankExprA1= rank(rowMeans(datExprA1p))
rankExprA2= rank(rowMeans(datExprA2p))
rankExprA3= rank(rowMeans(datExprA3p))
rankExprA4= rank(rowMeans(datExprA4p))
rankExprA5= rank(rowMeans(datExprA5p))
            
pdf("generalNetworkProperties_ensg.pdf", height=10, width=9) 

par(mfrow=c(2,2))
verboseScatterplot(rankExprA1,rankExprA2, xlab="Ranked Expression (A1)",
                   ylab="Ranked Expression (A2)")
verboseScatterplot(rankExprA1,rankExprA3, xlab="Ranked Connectivity (A1)",
                   ylab="Ranked Connectivity (A2)")
verboseScatterplot(rankExprA1,rankExprA4, xlab="Ranked Expression (A1)",
                   ylab="Ranked Expression (A2)")
verboseScatterplot(rankExprA1,rankExprA5, xlab="Ranked Connectivity (A1)",
                   ylab="Ranked Connectivity (A2)")
dev.off()        



write.csv(ans, file="ans.csv",quote=T)
