
data1<-read_excel("C:\\Users\\zxh\\Desktop\\v11.xlsx",sheet=1)

ggplot(data1,aes(x=ss,y=sqdj,fill=ss))+

geom_violin(alpha=0.99,width=1,adjust=1)+ #É¾³ıÁËscale='count'

geom_boxplot(width=0.05, fill="black",outlier.colour = NA) +
  
stat_summary(fun.y = median, geom="point", fill="white", shape=21, size=1.5)+

labs(x="Treatment",y="Index diameter after chemotherapy (mm)",title="C")+

theme_bw()+

theme(panel.background = element_rect(fill = "black"))+

scale_y_continuous(limit=c(5,25),breaks=c(5,8,10,15,20,25))+

theme(axis.title=element_text(size=5, family = "myFont",color="black"))+

theme(plot.title = element_text(size=8, family = "myFont",color="black", face = "bold",hjust=0,vjust=-0.5))+

theme(plot.margin=unit(c(0.01,0.01,0.01,0.01),'in'))+

theme(axis.text.x = element_text(size=5, family = "myFont",color="black", face = "plain"))+

theme(axis.text.y = element_text(size=5, family = "myFont",color="black", face = "plain"))+

theme(axis.ticks.length = unit(0.03,"in"))+

theme(axis.ticks=element_line(colour="black", size=.4, linetype=1, lineend=1))+

#theme(panel.spacing.x=unit(0.1, "lines"),panel.spacing.y=unit(1,"lines"))+

theme(legend.position = 'none')+# È¥µôÍ¼Àı

theme(panel.border = element_blank())+# È¥µôÍâ²ã±ß¿ò

theme(panel.background = element_rect(fill = "gray91"))+ #»ÒÉ«±³¾°

theme(panel.grid=element_blank())+

scale_fill_discrete(c=100, l=100)


ggsave("C:\\Users\\zxh\\Desktop\\Figure 2C.eps",device="eps",width = 1.6, height = 2,dpi =600)#±£´æÍ¼Æ¬

#savePlot(filename = "C:\\Users\\zxh\\Desktop\\Figure 2d.eps",
#type = c("eps"),
#device = dev.cur(),
#restoreConsole = TRUE)
