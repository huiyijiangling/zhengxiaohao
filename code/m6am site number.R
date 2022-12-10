#m6am site number
source(r"(C:\Users\zxh\Desktop\R\rsetup\theme_big_simple1.R)")
m6a50=readxl::read_excel(r"(C:\Users\zxh\Desktop\杜永星老师 提交\fto\fto文献\scDART-seq reveals distinct m 6 A signatures and mRNA methylation heterogeneity in single cells (sc m6a exp writer)\1-s2.0-S1097276521011436-mmc5.xlsx)",sheet = 6,skip = 1)
m6a50=head(m6a50,10)
# ggstatsplot::ggbarstats(m6a50, x = Gene, y = nSites,stat = 'identity')
library(ggplot2)
m6a50$Gene=factor(m6a50$Gene,levels = m6a50$Gene)
ggplot(m6a50,aes(x=Gene,y=nSites,fill=Gene))+geom_bar(stat = 'identity')+
ggsci::scale_fill_d3()+
  # scale_fill_manual(values = cell.col) +                                                  # 设置不同组对应的颜色

  # theme_classic() +                                                                       # 去除不必要的元素
  xlab("Gene symbol") +                                                                     # 修改x和y轴列名
  ylab("m6A sites number") +  
  theme_big_simple1()+                                       # 调整坐标轴标题字号
  theme(strip.background = element_rect(fill="grey", color = "white", size = 1),          # 设置顶部seurat_cluster显示为白框灰底
        strip.text = element_text(size = 12, colour="black"),                             # 设置顶部seurat_cluster字号和文字颜色
        legend.title = element_blank(),      
        axis.title = element_text(size = 15))
