theme_big_simple1 <- function(){ 
  theme_bw(base_size = 16, base_family="") %+replace% 
    theme(
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA),
      legend.title=element_text(vjust = 1,face="bold",family="sans",size = 24),
      legend.text=element_text(vjust = 1,face="bold",family="sans",size = 20),
      legend.justification = "right",
      legend.key.width = unit(40, "pt"),
      axis.line = element_line(color = "black", size = 1, linetype = "solid"),
      axis.ticks = element_line(colour = "black", size = 1),
      panel.background  = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),#panel.grid.major = element_line(color="lightgrey"),
      panel.border = element_blank(),
      # legend.position = ("bottom"),
      plot.title = element_text(face="bold",family="sans",size = 24, hjust = 0.0, vjust = 1.75),
      axis.text.x = element_text(face="bold",family="sans",hjust = 0.5,color="black", size = 20, margin = margin(t = 4, r = 0, b = 0, l = 0)),
      axis.text.y = element_text(face="bold",family="sans",hjust = 0.5,color="black", size= 20, margin = margin(t = 0, r = 4, b = 0, l = 0)),
      axis.title.x = element_text(face="bold",family="sans",hjust = 0.5,margin = margin(t = 10, r = 0, b = 0, l = 0), angle = 0, size = 24),
      axis.title.y = element_text(face="bold",family="sans",hjust = 0.5,margin = margin(t = 0, r = 10, b = 0, l = 0), angle = 90, size = 24),
      axis.ticks.length = unit(0.20, "cm"),
      strip.background = element_rect(color="black", size = 1, linetype="solid"),
      strip.text.x = element_text(face="bold",family="sans",size = 20, color = "black"),
      strip.text.y = element_text(face="bold",family="sans",size = 20, color = "black")
    )}




#C:\Users\zxh\Desktop\R\Ballinger\BallingerMack_NYBZase_2022\code\postprocess_RNAseq\useful
#small

theme_small_simple1 <- function(){ 
  theme_bw(base_size = 1.6, base_family="") %+replace% 
    theme(
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA),
      legend.title=element_text(vjust = 1,face="bold",family="sans",size = 2.4),
      legend.text=element_text(vjust = 1,face="bold",family="sans",size = 2),
      legend.justification = "right",
      legend.key.width = unit(4, "pt"),
      axis.line = element_line(color = "black", size = 1, linetype = "solid"),
      axis.ticks = element_line(colour = "black", size = 1),
      panel.background  = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),#panel.grid.major = element_line(color="lightgrey"),
      panel.border = element_blank(),
      # legend.position = ("bottom"),
      plot.title = element_text(face="bold",family="sans",size = 2.4, hjust = 0.0, vjust = 1.75),
      axis.text.x = element_text(face="bold",family="sans",hjust = 0.5,color="black", size = 2, margin = margin(t = 4, r = 0, b = 0, l = 0)),
      axis.text.y = element_text(face="bold",family="sans",hjust = 0.5,color="black", size= 2, margin = margin(t = 0, r = 4, b = 0, l = 0)),
      axis.title.x = element_text(face="bold",family="sans",hjust = 0.5,margin = margin(t = 10, r = 0, b = 0, l = 0), angle = 0, size = 2.4),
      axis.title.y = element_text(face="bold",family="sans",hjust = 0.5,margin = margin(t = 0, r = 10, b = 0, l = 0), angle = 90, size = 2.4),
      axis.ticks.length = unit(0.02, "cm"),
      strip.background = element_rect(color="black", size = 1, linetype="solid"),
      strip.text.x = element_text(face="bold",family="sans",size = 2, color = "black"),
      strip.text.y = element_text(face="bold",family="sans",size = 2, color = "black")
    )}