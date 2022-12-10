####===Figure 6===####
rm (list=ls())
gc()


####==loading data==####
library(dplyr)
library(ggsci)
library(ggradar)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(ggplot2)
load('data/CRISPR/crispr.Rdata')
load('data/sig/Stem.Sig.Rdata')
source('data/sig/sig_others.R')


####===Figure 6A===####

df <- crispr
df$mean <- apply(df,1,function(z){mean(z,na.rm=T)})
df <- df[order(df$mean),]
df <- df[c(1:8, (nrow(df)-7):nrow(df)),]
range(df,na.rm=T)
column_ha = HeatmapAnnotation(cohort = colnames(df))

Heatmap(as.matrix(df), name = "z scores",
        column_title = NULL,
        row_title = NULL,
        cluster_rows = F,
        cluster_columns = F,
        col = colorRamp2(c(-2, 0, 2), c("#3D7671", "grey", "#CB9B0C")),
        show_row_names = T, show_column_names = F,
        rect_gp = gpar(col = "black", lwd = 2),
        width = ncol(df)*unit(5, "mm"), 
        height = nrow(df)*unit(5, "mm"),
        na_col = 'white',
        column_names_side = c('top'),
        row_split = c(rep('a',8),rep('b',8)),
        column_split = c(rep('a',17),'b'),
        top_annotation = column_ha
)



###===ranking of genes ===###
rank <- apply(crispr,1,function(z){ mean(z,na.rm=T)})
rank <- data.frame(meanZ=rank,row.names = rownames(crispr))
rank$genes <- rownames(crispr)
rank <- rank[order(rank$meanZ),]
rank$order <- order(rank$meanZ)


####==Fig.6B==####

#signatures with immune resistant genes
num <- c('Stem.Sig', 
         'ImmuneCells.Sig',
         'TcellExc.Sig',
         'IMS.Sig',
         'LRRC15.CAF.Sig',
         'CRMA.Sig')


compare_crs <- data.frame(row.names = num)
row <- data.frame(row.names =  rank$genes)
ft <- data.frame()

for(i in c(0.01,0.02,0.03)){
  for(j in num){
    if ( j == "ImmuneCells.Sig"){
      sig <- readRDS('data/sig/ImSig.rds')
    }else if (j == 'TcellExc.Sig'){
      load("data/sig/FINAL_Tcell_Exclusion.sig.RData")
      sig <- exc.sig$exc.up
    } else if (j == 'TRS.Sig'){
      sig <- TRS.Sig
    } else if (j == 'LRRC15.CAF.Sig'){
      sig <- LRRC15.CAF.Sig
    } else if (j == 'IPRES.Sig'){
      sig <- IPRES.Sig
    } else if (j == "CRMA.Sig"){
      sig <- CRMA.Sig
    } else if (j == 'Stem.Sig') {
      sig <- Stem.Sig
    } 
    
    row[,j] <- ifelse(rank$genes %in% sig,'p','n')
    
    r <- rank[rank$genes %in% sig,]
    
    
    compare_crs[j,paste(i)] <- sum(r$order<=nrow(rank)*i) / nrow(r)
    
    ft[paste(i,'1'),j] <- sum(r$order<=nrow(rank)*i)
    ft[paste(i,'0'),j] <- nrow(r)-sum(r$order<=nrow(rank)*i)
    
    print(paste(j , i*100,'%' , sum(r$order<=nrow(rank)*i),nrow(r)))
    
  }
}

compare_crs <- as.data.frame(t(compare_crs))

compare_crs$sig <- as.character(paste('top',as.numeric(rownames(compare_crs))*100,'%','top-ranked genes'))


compare_crs$sig <- factor(compare_crs$sig,levels = compare_crs$sig)

compare_crs <- compare_crs[,c(ncol(compare_crs),1:(ncol(compare_crs)-1))]

compare_crs

ggradar(compare_crs,
        values.radar = c("0","3.6%","4.8%"),
        grid.min = 0,
        grid.mid = 0.036,
        grid.max = 0.048,
        group.colours = c('#83976B','#3D7671','#10332B')
        # )
        
) 



ft$all <- c(round(0.01*nrow(rank)),nrow(rank)-round(0.01*nrow(rank)),round(0.02*nrow(rank)),nrow(rank)-round(0.02*nrow(rank)),round(0.03*nrow(rank)),nrow(rank)-round(0.03*nrow(rank)))


fisher.test(ft[5:6,c(1,7)],alternative = "greater") 
# Immune-resistant genes (3% top-ranked genes) were significantly over-represented in Stem.Sig




##Fig. 6C
df_rank <- crispr

r <- rank[rank$genes %in% Stem.Sig,]
df_rank <- df_rank[rownames(r)[r$order<nrow(rank)*0.03],] #top 3% genes
df_rank$`mean Z Scores` <- rank[rownames(df_rank),'meanZ']
df_rank <- t(df_rank)

Heatmap(df_rank, name = "CRISPR immune Scores",
        row_title = NULL,
        cluster_rows = F,
        cluster_columns = F,
        col = colorRamp2(c(-2, 0, 2), c("#3D7671", "grey", "#CB9B0C")),
        show_row_names = T, show_column_names = T,
        width = ncol(df_rank)*unit(5, "mm"), 
        height = nrow(df_rank)*unit(5, "mm"),
        rect_gp = gpar(col = "black", lwd = 2),
        na_col = 'white',
        column_names_side = c('top')
)





















 

