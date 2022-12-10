library(readxl)
wuming=readxl::read_excel("1.xls",sheet = 1,na = "")
wuming$TAP=as.numeric(wuming$TAP)
wuming$`PGⅠ（胃蛋白酶原I）`=as.numeric(wuming$`PGⅠ（胃蛋白酶原I）`)
wuming$`PGII（胃蛋白酶原Ⅱ）`=as.numeric(wuming$`PGII（胃蛋白酶原Ⅱ）`)
wuming$`PGⅠ/PGⅡ`=as.numeric(wuming$`PGⅠ/PGⅡ`)
wuming$`CA19-9`=as.numeric(wuming$`CA19-9`)
wuming$CA125=as.numeric(wuming$CA125)
wuming$`CA72-4`=as.numeric(wuming$`CA72-4`)
wuming$CEA=as.numeric(wuming$CEA)
wuming$`CYFRA21-1`=as.numeric(wuming$`CYFRA21-1`)
wuming=as.data.frame(wuming)
box=wuming
# box <- gather(box, key = "variable",value ="value", 1:(nrow(ewid.data_gene_3)) )
# box=subset(box,!is.na(box$value))
library(ggpubr)

# google search : ggpubr boxplot add p-value
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
# ,outlier.shape = NA,order=
box$liange=box$良恶性


###############1
p <- ggboxplot(box, x = "liange", y = "CYFRA21-1",  fill = "liange",outlier.shape = NA,#fill="variable",
               font.label = list(size = 0.4 )#,palette=c(1,2,3),order=c("No changes","Shallow Deletion","Deep Deletion")#"#00AFBB","#E7B800",
)+
  #add = "jitter"
  # scale_color_manual(values=c(1,2,3,4))  +#"#E69F00", "#56B4E9"
  coord_cartesian(ylim = c(0, 20))+
  # theme(legend.position="left")+
  labs(title="",x="Tumor vs. Normal", y = "CYFRA21-1")+#rotate()
  rotate_x_text(90)
#  Add p-value
p + stat_compare_means(,method = "wilcox.test", label.y =20)#aes(group=cnv3), label.y = 
# ggpar(p,legend = "none")
ggsave("CYFRA21-1.pdf",width = 9.6,height=6,dpi = 600)

###################

library(tidyverse)
library(corrplot)
# p.mat=res1$p,

# cor=t(cor)
# cor=as.data.frame(cor)
# colnames(cor)=c("FTO",rownames(TCGA_immune))
cor=box[,c(22,6:14)]
cor=na.omit(cor)
colnames(cor)=paste0("V",1:10)
box
# cor=as.numercoric(cor)
res1 <- cor.mtest(cor, conf.level = .95, method = "spearman")
cor(cor) %>%   corrplot(method = "square",
                        type = "lower",tl.srt = 45,tl.col = "black")

pdf("CoR.pdf",width = 10,height=10)
pp1<- cor(cor, method ="spearman") %>%  corrplot(type="lower",method = "square",
                                                 insig="label_sig",p.mat = res1$p,outline = "white",
                                                 sig.level=c(0.01,.05), tl.cex=0.9 ,#.001,.01,
                                                 pch.cex=.9,pch.col="green",tl.srt=45,tl.col="black")
dev.off()








#+#tag="ENSG00000000003"+
# ??geom_dotplot
dev.off()
p+geom_point(aes( color=source2),shape=1,binaxis='y', stackdir='center', stackratio=0.1, dotsize=0.5,show.legend = T)+  guides(color = guide_legend(title = 'Source',nrow = 1, byrow = TRUE,legend.position = "top", direction = "vertical"),fill="none")+
  theme(legend.position = c(0.8,0.1), legend.background = element_blank())

ggsave('gene_up_all_GO_dotplot.pdf',width = 9.6,height=6,dpi = 600)




# median(wuming$`PGII（胃蛋白酶原Ⅱ）`,na.rm = T)
wuming_cancer=subset(wuming,wuming$良恶性==1)
apply(wuming_cancer,2,sum)
apply(wuming_cancer,2,function(x) table(is.na(x)))
apply(as.data.frame(wuming_cancer[,6:14]),2,function(x) median(as.numeric(x),na.rm = T))
x=as.data.frame(wuming_cancer[,6:14])
colnames(x)=paste0("V",1:9)
apply(x,2,function(x) class(x))
apply(x,2,function(x) median(x,na.rm = T))
wuming_cancer$`PGⅠ（胃蛋白酶原I）`
median(x,na.rm = T)
