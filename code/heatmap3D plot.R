load(file="resultTlist_frequency.Rdata")
resultTlist_frequency$sum=as.numeric(resultTlist_frequency$sum)
resultTlist_frequency$TCGA=as.numeric(resultTlist_frequency$TCGA)
mat=resultTlist_frequency[,c("e","miRNAs","sum")]
mat=tidyr::spread(mat,key = "miRNAs",value ="sum")
rownames(mat)=mat$e
uuuu=rownames(mat)
mat$e=NULL
mat[is.na(mat)] <- 0
# mat=as.matrix(mat)
mat=apply(mat,2,as.numeric)
rownames(mat)=uuuu
mat=as.data.frame(mat)
rowSums(mat)
mat=mat[order(rowSums(mat),decreasing = T),]
rowSums(mat)
colSums(mat)
mat=mat[,order(colSums(mat),decreasing = T)]
# mat=subset(mat,select=order(colSums(mat),decreasing = T))
colSums(mat)
rowSums(mat)
if(T){
  resultTlist_frequency_our3=resultTlist_frequency
  resultTlist_frequency_our3$sum=ifelse(resultTlist_frequency$TCGA==1&resultTlist_frequency_our3$sum>=3,1,0)
  eligible=resultTlist_frequency_our3[,c("e","miRNAs","sum")]
  eligible=tidyr::spread(eligible,key = "miRNAs",value ="sum")
  rownames(eligible)=eligible$e
  uuuu=rownames(eligible)
  eligible$e=NULL
  eligible[is.na(eligible)] <- 0
  # eligible=as.eligiblerix(eligible)
  eligible=apply(eligible,2,as.numeric)
  rownames(eligible)=uuuu
  eligible=as.data.frame(eligible)
  eligible=eligible[rownames(mat),colnames(mat)]
  identical(rownames(eligible),rownames(mat))
  identical(colnames(eligible),colnames(mat))

}
#155 miRNA eligble
#507 lcnRNA-mRNA eligble
# mat=mat[rowSums(eligible)>0,colSums(eligible)>0]
# eligible=eligible[rowSums(eligible)>0,colSums(eligible)>0]




# ha_column = HeatmapAnnotation(
# year = anno_text(year_text, rot = 0, location = unit(1, "npc"), just = "top")
# )
library(circlize)
library(ComplexHeatmap)
# col_fun = circlize::colorRamp2(c(0, 1, 3, 5), c("white", "cornflowerblue", "yellow","red"))
# col_fun = circlize::colorRamp2(c(0, 1,2, 3,4, 5), c("white", "cornflowerblue", "green","yellow","coral","red"))
col_fun = structure( c("white", "cornflowerblue", "green","yellow","coral","red"),names =as.character(0:5))
ht_opt$TITLE_PADDING = unit(c(15, 2), "mm")

# mat1 = readRDS(system.file("extdata", "measles.rds", package = "ComplexHeatmap"))
# mat=mat[1:1000,1:5]


library(tidyverse)
library(ComplexHeatmap)

# making the dataframe




if(T){
  m =mat
  tbl = lapply(m, table)
  mx = matrix(0, nrow =6, ncol =  ncol(mat))
  rownames(mx) = 0:5
  colnames(mx) = colnames(m)
  for(i in 1:ncol(mx)){
    mx[names(tbl[[i]]), i] = tbl[[i]]
  }
  Heatmap(mx, cluster_rows = FALSE, cluster_columns = FALSE)
  mx=mx[-1,]
  #
  tbl = lapply(as.data.frame(t(m)),table)
  my = matrix(0, nrow =6, ncol =  nrow(mat))
  rownames(my) = 0:5
  colnames(my) = rownames(m)
  for(i in 1:ncol(my)){
    my[names(tbl[[i]]), i] = tbl[[i]]
  }
  Heatmap(my, cluster_rows = FALSE, cluster_columns = FALSE)
  my=my[-1,]
}
col_fun15 = circlize::colorRamp2(c(0, 10,50, 100,500, 1000),  c("white", "cornflowerblue", "green","yellow","coral","red"))
col_fun135 = circlize::colorRamp2(c(0, 5,10,15 ,20, 30),  c("white", "cornflowerblue", "green","yellow","coral","red"))

#
mx_top = ComplexHeatmap::Heatmap3D(as.matrix(mx[rev(rownames(mx)),]), name = "Evidences per miRNA", col = col_fun15,
                                   cluster_rows = F,
                                   cluster_columns = FALSE, show_row_dend = F,# rect_gp = gpar(col= "white"),
                                   show_column_names = F,
                                   row_title = "Evidences per miRNA",
                                   row_title_side = c("left"),
                                   show_row_names = T,
                                   row_title_gp = gpar(fontsize = 13.2),
                                   column_labels = stringr::str_split(colnames(mat),"hsa-",simplify = T)[,2],
                                   row_names_side = "right", row_names_gp = grid::gpar(fontsize = 12),
                                   column_names_gp = gpar(fontsize = 4),
                                   use_raster=F,
                                   # column_title = 'Heatmap of mRNA-miRNA-lncRNA regulation pairs using in silico predictions',
                                   # top_annotation = ha1,
                                   # right_annotation = ha2,
                                   
                                   # bottom_annotation = ha_column,
                                   heatmap_legend_param = list(
                                     # at = c(0, 10,50, 100, 1000),
                                                               # labels =as.character(c("0", "10","50","100", "1000")),
                                                               direction = "horizontal"),
                                   bar_rel_width = 1, bar_rel_height = 1, bar_max_length = unit(1, "cm"),
                                   height=1.5
)
# draw(mx_top)
ha1 = HeatmapAnnotation(
  "Number of prediction per miRNA" = anno_barplot(
    colSums(mat), 
    bar_width = 1, 
    gp = gpar(col = "white", fill = "#FFE200"), 
    border = FALSE,
    axis_param = list(at = c(0, 2e2,5e2,1e3,2e3,3e3),
                      labels = c("0", "200", "500","1000","2000","3000")),
    height =  unit(1, "cm"),
  ), 
  "Eligible miRNAs"=colSums(eligible),
  show_annotation_name = T,
  annotation_legend_param = list(direction = "horizontal",nrow=1),
  col=list("Eligible miRNAs"=col_fun15,
           "Number of prediction per miRNA"="#FFE200")
  )

ha2 = rowAnnotation(
  "Evidences per LncRNA-mRNA pair"=as.matrix(t(my)),
  col=list("Evidences per LncRNA-mRNA pair"=col_fun135,
           "Eligible LncRNA-mRNA pairs"=col_fun135,
           "Number of predicted evidences"="#FFE200"),
  # name="Per LncRNA-mRNA relationships",
  "Eligible LncRNA-mRNA pairs"=rowSums(eligible),
  "Number of predicted evidences" = anno_barplot(
    rowSums(mat),
    # names="vvv",
    bar_width = 1, 
    gp = gpar(col = "NA", fill = "#FFE200"), 
    border = FALSE,
    axis_param = list(at = c(0, 10, 25,50,100,150),
                      labels = c("0",  "10", "25","50","100","150")),
    width = unit(1, "cm"),

  ),    annotation_legend_param = list(direction = "horizontal",nrow=1),show_annotation_name =T,annotation_name_side = "top")

# ha3 = HeatmapAnnotation(label = anno_mark(at = which(colSums(eligible)>0), labels = colnames(eligible)[colSums(eligible)>0],side="bottom"))

# ha4 = rowAnnotation(label = anno_mark(at = which(rowSums(eligible)>0), labels = rownames(eligible)[rowSums(eligible)>0]))
# ha3=HeatmapAnnotation(mx_right = ComplexHeatmap::Heatmap(as.matrix(t(my)), name = "cases", col = col_fun135,cluster_rows = F,
#                                    cluster_columns = FALSE, show_row_dend = F,# rect_gp = gpar(col= "white"),
#                                    show_column_names = F,
#                                    show_row_names = FALSE,
#                                    column_labels = stringr::str_split(colnames(mat),"hsa-",simplify = T)[,2],
#                                    row_names_side = "left", row_names_gp = grid::gpar(fontsize = 6),
#                                    column_names_gp = gpar(fontsize = 6),
#                                    use_raster=F,
#                                    # column_title = 'Heatmap of mRNA-miRNA-lncRNA regulation pairs using in silico predictions',
#                                    # top_annotation = ha1,
#                                    # right_annotation = ha2,
#                                    # bottom_annotation = ha_column,
#                                    heatmap_legend_param = list(at = c(0, 5,10, 20, 30),
#                                                                labels =as.character(c("0", "5","10","20", "30"))),
#                                    # bar_rel_width = 1, bar_rel_height = 1, bar_max_length = unit(1, "cm"),
#                                    # height=2
#                                    ))
# library(doParallel) #并行处理包
# cl <- makeCluster(30)#makeCluster(detectCores())
# registerDoParallel(cl)
# ht <- list()
# ht=lapply(splitIndices(ncol(mat), 30),unlist)
# ht=lapply(1:30,function(x) mat[,ht[[x]]])
# 
# ht <- foreach(i=30) %dopar% ComplexHeatmap::Heatmap3D(as.matrix(ht[[i]]), name = "cases", col = col_fun,
#                                                       cluster_columns = FALSE, show_row_dend = FALSE, #rect_gp = gpar(col= "white"), 
#                                                       show_column_names = FALSE,
#                                                       show_row_names = FALSE,
#                                                       row_names_side = "left", row_names_gp = grid::gpar(fontsize = 8),
#                                                       use_raster=T,
#                                                       column_title = 'Heatmap of mRNA-miRNA-lncRNA regulation pairs using in silico predictions',
#                                                       # top_annotation = ha1,
#                                                       # bottom_annotation = ha_column,
#                                                       heatmap_legend_param = list(at = 0:5, 
#                                                                                   labels = as.character(0:5)),
#                                                       # new arguments for Heatmap3D()
#                                                       bar_rel_width = 1, bar_rel_height = 1, bar_max_length = unit(2, "cm") 
# )
# stopCluster(cl)
# + ha2
# ht_list=eval(parse(text =paste0("ht","[[",1:length(ht),"]]",collapse = "+")))
# mat22=mat[1:100,1];mat22
# deparse(substitute(mat))
# pheatmap::pheatmap(mat = as.matrix(mat),cluster_rows = F,cluster_cols = F, show_colnames = F, show_rownames = F)

# frequencyHeatmap(as.matrix(mat), use_3d = T,ylim = c(0:5))


ht = ComplexHeatmap::Heatmap(as.matrix(mat), name = "Supported number of datasets", col = col_fun,cluster_rows = F,
                             cluster_columns = FALSE, show_row_dend = F,# rect_gp = gpar(col= "white"),
                             show_column_names = F,
                             show_row_names = FALSE,
                             column_labels = stringr::str_split(colnames(mat),"hsa-",simplify = T)[,2],
                             # column_names_side = "top",
                             row_names_side = "left", row_names_gp = grid::gpar(fontsize = 1),
                             column_names_gp = gpar(fontsize = 1),
                             use_raster=F,
                             # column_title = 'Heatmap of mRNA-miRNA-lncRNA regulation pairs using in silico predictions',
                             
                             # bottom_annotation = ha1,
                             right_annotation = ha2,
                             # left_annotation = ha4,
                             # bottom_annotation = ha_column,
                             heatmap_legend_param = list(at = 0:5,labels = as.character(0:5),direction = "horizontal",nrow=1),
                             width = 8,
                             height = 4
)
# ht = Heatmap(m)
# ha = HeatmapAnnotation(foo = 1:10)
# ht = attach_annotation(ht, ha)
# ht = attach_annotation(ht,ha,side = "right")
# ht=ht+ha2


# draw(mx_right, ht_gap = unit(0.3, "mm"))
#
# draw(ha1%v% wll+ha1, ht_gap = unit(0.3, "mm"))
wll= ha1%v%mx_top%v% ht#%v% ha3
pdf(file="if.pdf",width = 12,height = 10)
draw(wll, ht_gap = unit(0.2, "mm"),     
     annotation_legend_side = "bottom",
     heatmap_legend_side = "bottom",
     merge_legend = TRUE)
# draw(ht, ht_gap = unit(0.3, "mm"))
# dev.off()
# 
# draw(ht, ht_gap = unit(0.3, "mm"))
decorate_heatmap_body("Supported number of datasets", {
  
  i = which(colnames(mat) == colnames(mat)[max(1:quantile(1:ncol(mat))[2])])
  x = i/ncol(mat)
  grid.lines(x=c(x, x),y= c(0, 1), gp = gpar(lwd = 2, lty = 2))
  grid.text("75%", x, unit(1, "npc") + unit(5, "mm"))
})
decorate_heatmap_body("Supported number of datasets", {
  
  i = which(rownames(mat) == rownames(mat)[max(1:quantile(1:nrow(mat))[2])])
  y = i/nrow(mat)
  #这里有问题所以我手写了
  grid.lines(c(0, 1), c(0.75,0.75), gp = gpar(lwd = 2, lty = 2))
  grid.text("75%", y=0.75, unit(1, "npc") + unit(5, "mm"))
})
dev.off()

# new arguments for Heatmap3D()
# bar_rel_width = 1, bar_rel_height = 1, bar_max_length = unit(2, "cm")
pdf(file="if.pdf",width = 8,height = 8)
ht_list =Heatmap3D(as.matrix(mat), name = "cases", col = col_fun,
                   cluster_columns = FALSE, show_row_dend = FALSE, rect_gp = gpar(col= "white"),
                   show_column_names = FALSE,
                   show_row_names = FALSE,
                   row_names_side = "left", row_names_gp = gpar(fontsize = 8),
                   use_raster=T,
                   column_title = 'Heatmap of mRNA-miRNA-lncRNA regulation pairs using in silico predictions',
                   top_annotation = ha1, 
                   # bottom_annotation = ha_column,
                   heatmap_legend_param = list(at = 0:5, 
                                               labels = as.character(0:5)),
                   # new arguments for Heatmap3D()
                   bar_rel_width = 1, bar_rel_height = 1, bar_max_length = unit(2, "cm") 
)+ ha2
draw(ht_list, ht_gap = unit(0, "mm"))
# decorate_heatmap_body("cases", {
#   i = which(colnames(mat) == "1961")
#   x = i/ncol(mat)
#   grid.lines(c(x, x), c(0, 1), gp = gpar(lwd = 2, lty = 2))
#   grid.text("Vaccine introduced", x, unit(1, "npc") + unit(5, "mm"))
# })
dev.off()
pdf(file="if.pdf",width = 8,height = 8)
kkkkk
dev.off()
eval(parse(text = paste0(,collapse = "+")))
ht_opt$ht_gap = unit(0, "mm")
dev.off()



ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm, width = unit(4, "cm"))
ht2 = Heatmap(mat2, name = "runif", col = col_runif, width = unit(6, "cm"))
ht3 = Heatmap(le, name = "letters", col = col_letters, width = unit(1, "cm"))
ht1 + ht2 + ht3

