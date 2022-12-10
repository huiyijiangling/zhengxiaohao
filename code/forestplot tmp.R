library(forestplot)


rs_forest <- read.csv('demo.csv',header = FALSE)
rs_forest
#              V1           V2       V3              V4   V5   V6   V7
# 1         Group Case/Control        p       OR(95%CI)   NA   NA   NA
# 2           Age                                         NA   NA   NA
# 3           <50      127/211 2.43E-11 4.24(2.84~6.66) 4.24 2.84 6.66
# 4         50-59      115/244 4.14E-12   5.48(3.5~9.2) 5.48 3.50 9.20
# 5          >=60      132/237 1.25E-08 2.74(1.96~3.93) 2.74 1.96 3.93
# 6           Sex                                         NA   NA   NA
# 7          Male      323/318 1.00E-15  3.46(2.71~4.5) 3.46 2.71 4.50
# 8        Female       51/374 3.14E-09 4.21(2.69~7.01) 4.21 2.69 7.01
# 9       Smoking                                         NA   NA   NA
# 10          Yes      122/171 2.06E-07 2.56(1.82~3.71) 2.56 1.82 3.71
# 11           No      252/521 1.00E-15 4.16(3.17~5.59) 4.16 3.17 5.59
# 12 Hypertension                                         NA   NA   NA
# 13          Yes      193/207 3.31E-12  3.14(2.31~4.4) 3.14 2.31 4.40
# 14           No      181/485 1.00E-15 3.96(2.95~5.45) 3.96 2.95 5.45
# 15     Drinking                                         NA   NA   NA
# 16          Yes       83/111 2.44E-04 2.42(1.55~4.01) 2.42 1.55 4.01
# 17           No      291/581 1.00E-15 3.83(3.01~4.96) 3.83 3.01 4.96
演示数据集demo的格式
forestplot(labeltext = as.matrix(rs_forest[,c(1:4)]),
           mean = rs_forest$V5, 
           lower = rs_forest$V6, 
           upper = rs_forest$V7, 
           is.summary=c(T,T,F,F,F,T,F,F,T,F,F,T,F,F,T,F,F),
           zero = 1, 
           boxsize = 0.5, 
           lineheight = unit(6,'mm'),
           graphwidth = unit(35,'mm'),
           colgap = unit(2,'mm'),
           lwd.zero = 2,
           lwd.ci = 2,
           col=fpColors(box='#ef233c',summary="#9e2a2b",lines = 'blue1',zero = 'grey'),
           xlab="Odds Ratio",
           lwd.xaxis=2,
           lty.ci = "solid", ci.vertices.height = 0.15,clip=c(1,13),xticks.digits = 6,
           graph.pos = 4)

forestplot(labeltext = as.matrix(rs_forest[,c(1:4)]),
           mean = rs_forest$V5, 
           lower = rs_forest$V6, 
           upper = rs_forest$V7, 
           is.summary=c(T,T,F,F,F,T,F,F,T,F,F,T,F,F,T,F,F),
           graph.pos=4, #为Pvalue箱线图所在的位置
           #定义标题
           title="Hazard Ratio Plot",
           ##定义x轴
           xlab="    <---PCI Better---   ---Medical Therapy Better--->",
           ##根据亚组的位置，设置线型，宽度造成“区块感”
           #hrzl_lines=list("2" = gpar(lwd=1, lineend="butt", columns=c(1:5),col="#99999922"),
           #                "18" = gpar(lwd=1, lineend="butt", columns=c(1:5),col="#99999922")),
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1.25),
                          ticks=gpar(cex=1.1),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2)),
           ##fpColors函数设置颜色
           col=fpColors(box="#1c61b6", lines="#1c61b6", zero = "gray50"),
           #箱线图中基准线的位置
           zero=1,
           cex=0.9, lineheight = "auto",
           colgap=unit(8,"mm"),
           #箱子大小，线的宽度
           lwd.ci=2, boxsize=0.5,
           #箱线图两端添加小竖线，高度
           ci.vertices=TRUE, ci.vertices.height = 0.4)

3. 森林图结果解读
（1）森林图中横短线与中线相交表示无统计学意义；
（2）95% CI上下限均>1，即在森林图中，其95% CI横线不与无效竖线相交，且该横线落在无效线右侧时，说明该指标大于竖线代表的结局；
（3）95% CI上下限均＜1，即在森林图中，其95% CI横线不与无效竖线相交，且该横线落在无效线左侧时，说明该指标小于于竖线代表的结局。