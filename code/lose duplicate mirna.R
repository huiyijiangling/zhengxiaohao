# mirna 去重
rm(list=ls())
options(stringsAsFactors = F)
load("C:/Users/zxh/Desktop/R/second/STAD_miRNA_rpm.Rdata")

dat=mir_exp
# 去掉蜡块标本
# dat=dat[,substr(colnames(dat),16,17)!="B-"] 
dat=as.data.frame(t(dat))
# dat=log2(dat+1)
dat$name0121=substr(rownames(dat),1,21)
dat$name2225=substr(rownames(dat),22,25)
du=dat[duplicated(dat$name0121),]
du=dat[dat$name0121 %in% du$name0121,]
du=du[du$name2225!="A360",]
dat=dat[,-(ncol(dat)-1):-ncol(dat)]
dat=dat[!rownames(dat) %in% rownames(du),]
mir_exp=as.data.frame(t(dat))
save(mir_exp,file="STAD_miRNA_rpm_477.Rdata")
# write.csv(dat,"dat.csv")#请手动删除.其实我看了数据结构，重做均用了A360
# 
# #####################################################################################
# dat=read.csv("dat477.csv",header = F)
# rownames(dat)=dat[,1]
# dat=dat[,-1]
# # 后面要做pldyr 不能有-.
# colnames(dat)=dat[1,]
# dat=dat[-1,]
# mir_exp=as.data.frame(t(dat))


# mirna 去重
rm(list=ls())
options(stringsAsFactors = F)
load("C:/Users/zxh/Desktop/R/second/STAD_miRNA_count.Rdata")

dat=mir_exp
# 去掉蜡块标本
# dat=dat[,substr(colnames(dat),16,17)!="B-"] 
dat=as.data.frame(t(dat))
# dat=log2(dat+1)
dat$name0121=substr(rownames(dat),1,21)
dat$name2225=substr(rownames(dat),22,25)
du=dat[duplicated(dat$name0121),]
du=dat[dat$name0121 %in% du$name0121,]
du=du[du$name2225!="A360",]
dat=dat[,-(ncol(dat)-1):-ncol(dat)]
dat=dat[!rownames(dat) %in% rownames(du),]
mir_exp=as.data.frame(t(dat))
save(mir_exp,file="STAD_miRNA_count_477.Rdata")

