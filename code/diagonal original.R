library(ComplexHeatmap)
library(circlize)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
# 右下侧的数据 其实反了
up <- read.table(r"(C:\Users\zxh\Desktop\R\paad-tcga-gtex\easy_input_amp.txt)",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
up[1:3, 1:3]

##              KIRC      KIRP       KICH
## ALKBH5 0.04924242 0.5902778 0.00000000
## FTO    0.18560606 0.5243056 0.31818182
## ZC3H13 0.03977273 0.1006944 0.03030303

# 左上角的数据
dn <- read.table(r"(C:\Users\zxh\Desktop\R\paad-tcga-gtex\easy_input_del.txt)",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
dn[1:3, 1:3]

##              KIRC       KIRP       KICH
## ALKBH5 0.08522727 0.04861111 0.75757576
## FTO    0.04166667 0.02430556 0.07575758
## ZC3H13 0.15340909 0.08333333 0.66666667

# 检验两个矩阵列名和行名是一致的
identical(rownames(up), rownames(dn))

## [1] TRUE

identical(colnames(up), colnames(dn))

## [1] TRUE


RightOrder <- rev(rownames(up))
up <- up[RightOrder,]
dn <- dn[RightOrder,]
#其实具体看有没有空值，如果没有控制其实就很好看了 不一定非要用 viridis
UpColor <- colorRamp2(breaks = seq(-0.6,0.6, by = 0.1), colors = viridis::viridis(13))#"#AB221F"
DnColor <- colorRamp2(breaks = seq(-0.6,0.6, by = 0.1), colors = viridis::viridis(13))#"#3878C1"
# UpColor <- colorRamp2(breaks = c(-1,0,1), colors = c("blue","#FFFADD","red"))#"#AB221F"
# DnColor <- colorRamp2(breaks = c(-1,0,1), colors = c("blue","#FFFADD","red"))#"#3878C1"

# UpColor <- colorRamp2(breaks = c(0, 1), colors = c("#FFFADD","red"))#"#AB221F"
# DnColor <- colorRamp2(breaks = c(0, 1), colors = c("#FFFADD","blue"))#"#3878C1"

## 使用up数据集来产生数据
Heatmap(up,
        column_title = "Copy number variation across cancer types", ## 列的标题
        rect_gp = gpar(type = "none"),  #绘制空的数据框
        show_heatmap_legend = F, ##是否显示基本的注释说明
        cluster_rows = T, cluster_columns = T, ## 是否对行列进行聚类
)



row_an <-  HeatmapAnnotation(type = c(rep("R", 10), rep("W", 8), rep("E", 2)), ##注释信息的内容。
                             show_annotation_name = F, ## 是否显示注释的标题
                             col = list(type = c("R" = "#5AC9FA", "W" = "#FAC67A", "E" = "#51B743")), ## 注释信息的颜色
                             show_legend = T,  ## 是否显示注释信息的说明
                             annotation_legend_param = list(title = "m6A group",nrow = 1), ## 注释信息图例的个性化说明，nrow表示把所有分类的图例放到一行。
                             which = "row" #对行或者列进行注释
)

DiagFunc <- function(up, down){
  function(j, i, x, y, width, height, fill){
    grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width),
                 unit.c(y - 0.5*height, y + 0.5*height, y + 0.5*height),
                 gp = gpar(fill = DnColor(down[i, j]), col = "grey"))
    grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width),
                 unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height),
                 gp = gpar(fill = UpColor(up[i, j]), col = "grey"))
  }
}
uppvalue=up
downpavlue=dn
DiagFunc <- function(up,uppvalue,down,downpavlue){
  function(j, i, x, y, width, height, fill){
    grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width),
                 unit.c(y - 0.5*height, y + 0.5*height, y + 0.5*height),
                 gp = gpar(fill = DnColor(down[i, j]), col = "grey"))
    if (downpavlue[i, j] < 0.001) {
      grid.text("***", x, y,vjust = -0.2,rot = 0)
    }
    if (downpavlue[i, j]< 0.01&downpavlue[i, j]>=0.001) {
      grid.text("**", x, y,vjust = -0.2,rot = 0)
    }
    if (downpavlue[i, j]< 0.05&downpavlue[i, j]>=0.01) {
      grid.text("*", x, y,vjust = -0.2,rot = 0)
    }
    grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width),
                 unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height),
                 gp = gpar(fill = UpColor(up[i, j]), col = "grey"))
    if (uppvalue[i, j] < 0.001) {
      grid.text("***", x, y,vjust = 1.5,rot = 0)
    }
    if (uppvalue[i, j]< 0.01&uppvalue[i, j]>=0.001) {
      grid.text("**", x, y,vjust = 1.5,rot = 0)
    }
    if (uppvalue[i, j]< 0.05&uppvalue[i, j]>=0.01) {
      grid.text("*", x, y,vjust = 1.5,rot = 0)
    }
  }

}

p1 <- Heatmap(as.matrix(up), column_title = "Copy number variation across cancer types",
              rect_gp = gpar(type = "none"),
              show_heatmap_legend = F,
              cluster_rows = F,
              cluster_columns = F, ##绘制空的热图框
              left_annotation = row_an, ##添加左侧注释信息
              cell_fun = DiagFunc(up = up,uppvalue=uppvalue, down = dn,downpavlue=downpavlue) ## 绘制表格内的内容
); p1

lgd <- list(Legend(title = "CNV Gain", ## 注释的标题
                   col_fun = UpColor, ## 注释的颜色
                   at = c(0,0.5,1), ## 注释刻度的分组
                   direction = "horizontal" ## 注释的方向
)
,Legend(title = "CNV Loss", ## 注释的标题
        col_fun = DnColor, ## 注释的颜色
        at = c(0,0.5,1), ## 注释刻度的分组
        direction = "horizontal" ## 注释的方向
))

lgd <- list(Legend(title = "CNV Gain", ## 注释的标题
                   col_fun = UpColor, ## 注释的颜色
                   at = c(-0.6,-0.3,0,0.3,0.6), ## 注释刻度的分组
                   direction = "horizontal" ## 注释的方向
))

pdf(file =paste0("对角线热图.pdf"),height = 10,width = 15)
draw(p1, annotation_legend_list = lgd,
     annotation_legend_side = "bottom",
     heatmap_legend_side = "bottom",
     merge_legend = TRUE)
dev.off()

