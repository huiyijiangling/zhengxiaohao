#最后没用上
library("meta")
library("mvmeta")
library("forestplot")
col <- seq(from=7,to=122,by=5)
data(Olkin1995)
data(Fleiss93)
data(Fleiss93)

# m1 <- metabin(event.e, n.e, event.c, n.c,
#               
#               data = Fleiss93, studlab = study,
#               
#               sm = "RR", method = "I")
sub_log2=sub_log
sub_log=sub_log2
# sub_log=subset(sub_log2,subgroups=="Age2")
sub_log$level=as.character(sub_log$level)
m1 <- metabin(event_Delayed,
             sample_Delayed,
              event_Normal,
              sample_Normal,
              data = sub_log, studlab =level ,subgroup =subgroups,
              sm = "RR", method = "I")
m1
m2=update(m1,byvar =subgroups ,print.byvar=T)
m2
forest((m1),layout = "RevMan5")
forest(metainf(m1),layout = "RevMan5")

metainf(m1, pooled = "random")#revman5



forest(metainf(m2))
pdf("year.pdf",width = 12,height = 8)
forest((m1),layout = "RevMan5")
# forest(metainf(m1), layout = "RevMan5")
dev.off()
sub_log=subset(sub_log,subgroups=="pstageH")
sub_log$level=as.character(sub_log$level)
m1 <- metabin(event_Delayed,
              sample_Delayed,
              event_Normal,
              sample_Normal,
              data = sub_log, studlab = level,
              sm = "RR", method = "I")
pdf("pstageH.pdf",width = 10,height = 5)
forest(metainf(m1), layout = "RevMan5")
dev.off()


metainf(m1, sortvar = study)

metainf(m1, sortvar = 7:1)



m2 <- update(m1, title = "Fleiss93 meta-analysis",
             
             backtransf = FALSE)

metainf(m2)


# 新建一个matrix放每个类别的分析结果，命名为result。我这里保存了5列，分别是类别名factor、该类别异质性检验的I^2值、p-value和自由度，最后的Peter-P是发表偏倚的peters检验

result <- matrix(c(1:5*length(col)),nrow = length(col), ncol = 5)
result <- as.data.frame(result)
colnames(result) <- c("factor", "I^2", "P-value", "df", "Peter-P")

# 本来是一个循环可以把前面col提出的所有类别一次性完成分析的，但是考虑到每张图里可能要显示不同的内容，所以自己定义了变量i，每次更改i值来做不同类别的分析。下面以i=1为例做第一个类别变量--性别的相关分析。

i <- 1


# 异质性检验

# 利用meta包中的metabin函数进行异质性检验（前面i*5+3到i*5+6分别指定data.raw中性别对应的四列数据[即8、9、10、11列]，"OR"是指定效应量合并的方法，studlab指定下一步森林图中显示在前面的文献名，comb.random = TRUE和comb.fixed = TRUE表示随机效应模型和固定效应模型都做）：
f1 <- metabin(data.raw[,i*5+3], data.raw[,i*5+4], data.raw[,i*5+5], data.raw[,i*5+6], data=data.raw,sm="OR", studlab = paste(data.raw[,3],data.raw[,1]), comb.random = TRUE, comb.fixed = TRUE)


跑出来结果如下：

图片



将需要的信息提出来存到result中对应的位置：

result[i,1] <- colnames(data.raw[col[i]])
result[i,2] <- f1$I2
result[i,3] <- f1$pval.Q
result[i,4] <- f1$df.Q


现在result如下：

图片




森林图

前面已经将meta分析的结果赋给f1这个变量，直接用forest函数对f1做图即可（lab.e和lab.c设置显示在上面的标签，后面两个参数设置是否在图的左下角显示模型效应值）：
jpeg("1性别.jpeg",height=650,width=700)
forest(m1, lab.e = "men", lab.c = "women",test.overall.fixed=TRUE, test.overall.random=TRUE)
dev.off()


图如下：

图片


敏感性分析

直接调用meta包的metainf函数进行敏感性分析，方法有很多种，这里是将每个研究逐一排除做敏感性分析：
inf <- metainf(m1, pooled="fixed")
结果如下：

图片



如果直接用meta包里的forest函数对敏感性分析结果画森林图：

forest(metainf(m1), comb.fixed=TRUE)
结果如下：

图片



可以调的参数比较少（也可能是我没找到吧），上面显示的信息只有study、OR和95%CI，但是如果希望把逐个剔除研究做的敏感性分析得到的异质性I^2值和p值也显示上去，可以用forestplot包。



首先把前面敏感性分析的结果整理出来：

rf <- as.data.frame(matrix(c(1:7*(length(inf$studlab))+1), nrow = length(inf$studlab)+1, ncol = 7))
colnames(rf) <- c("study","I^2", "p.value", "95%CI", "OR", "lower", "upper")
rf[c(2:(1+length(inf$studlab))),1] <- inf$studlab
rf[c(2:(1+length(inf$studlab))),5] <- exp(inf$TE)
rf[c(2:(1+length(inf$studlab))),6] <- exp(inf$lower)
rf[c(2:(1+length(inf$studlab))),7] <- exp(inf$upper)
rf[c(2:(1+length(inf$studlab))),2] <- paste(round(inf$I2,4)*100,"%",sep="")
rf[c(2:(1+length(inf$studlab))),3] <- inf$p.value
rf[c(2:(1+length(inf$studlab))),4] <- paste(paste(paste(round(rf[c(2:(1+length(inf$studlab))),5],3), sep="(", round(rf[c(2:(1+length(inf$studlab))),6],3)), sep="-", round(rf[c(2:(1+length(inf$studlab))),7],3)),sep="",")")
for(m in c(2:(1+length(inf$studlab)))){
  if(is.na(rf[m,7])){
    rf[m,] <-NA
  }
}
rf[1,] <- c("study", "I^2", "Z-test p.value", "OR(95%CI)", NA, NA, NA)


整理出来，格式如下：

图片



用forestplot包作图：

jpeg("1性别敏感性分析.jpeg",height=400,width=650)
forestplot(labeltext = as.matrix(rf[,1:4]),mean = as.numeric(rf$OR), lower = as.numeric(rf$lower), upper=as.numeric(rf$upper), is.summary = c(rep(FALSE,length(rf[,1])-1),TRUE), zero = 1, boxsize = 0.2, graph.pos = 4, title = "性别敏感性分析")
dev.off()


结果如下：

图片




发表偏倚检验

先画漏斗图看一下：

jpeg("1性别漏斗图1.jpeg",height=300,width=400)
funnel(f1)
dev.off()
结果如下：

图片

用metabias进行发表偏倚的检验，使用Peters检验（需要研究数≥10），并把检验结果的p值保存到result里：

bias <- metabias(f1, method.bias="peters")
result[i,5] <- bias$p.value


结果如下：

图片



从上面漏斗图和peters检验如果发现该meta分析存在发表偏倚的可能性，需要进行校正，以trim and filled方法为例：

tf1 <- trimfill(f1, comb.fixed=TRUE)
summary(tf1)


结果如下：

图片



校正后再画漏斗图看一下：

jpeg("1性别漏斗图2.jpeg",height=300,width=400)
funnel(tf1)
dev.off()
如下：

图片



最后可以把result结果导出来：

write.csv(resultI,"result.csv", row.names = TRUE)
