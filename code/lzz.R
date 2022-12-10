library(readxl)
qwe1=readxl::read_xlsx("工作簿1(已自动还原).xlsx",col_names = F)
qwe1=Reduce(intersect,list(qwe1$...1,qwe1$...2,qwe1$...3))
qwe2=readxl::read_xlsx("肝癌(1).xlsx",col_names = F)
qwe2=as.data.frame(qwe2)
in1=qwe2[qwe2$...1%in%qwe1,]
out1=qwe2[!qwe2$...1%in%qwe1,]
write.csv(in1,file = "in.csv")
write.csv(out1,file = "out.csv")
write.csv(qwe1,file = "common.csv")
