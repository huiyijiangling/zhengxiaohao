anova(myfit,myfit_reduced,test="Chisq")
稳健Logistic回归：当拟合Logistic回归模型数据出现离群点和强影响点时，robust包中的glmRob()函数可用来拟合稳健的广义线性模型

多项分布回归：当响应变量包含两个以上的无序类别（例如已婚、寡居、离婚）时，可使用mlogit包中的mlogit()函数拟合多项Logistic回归

序数Logistic回归：当响应变量是一组有序的类别（例如信用为好/良/差）时，可使用rms包中的lrm()函数拟合序数Logistic回归。

# 360doc.com/content/20/0513/15/52334415_912090926.shtml

当其值大于1非常多时，则认为存在过度离势，同样的，也可以对其进行检验，
pchisq(summary(fit.od)$dispersion * fit$df.residual, fit$df.residual, lower = F)

当出现过度离势时，可使用 family=quasibinomial()对family = binomial()的部分进行替换。

###################
library(ipflasso)
pflist<-list(c(1,1),c(1,2),c(1,4),c(1,8),c(1,16),c(2,1),c(4,1),c(8,1),c(16,1),c(1,32),c(1,64))

set.seed(1)
mod.ipf<-cvr2.ipflasso(Y=y,X=x,family="binomial",standardize=T,
                       blocks=list(block1=1:3,block2=(3+1):ncol(x)),
                       pflist=pflist,nfolds=3,ncv=1,type.measure='class')
coeff <- list()
coeff$ipf<-mod.ipf$a[[mod.ipf$ind.bestpf]]$coeff[,mod.ipf$ind.bestlambda]
names(coeff$ipf)
compAUC<-function(coeff,newdata,y)
{
  linpred<-newdata%*%coeff
  my.auc(linpred,y)
}
lapply(coeff,compAUC,newdata=x,y=y)
# auc
lapply(coeff,compAUC,newdata=as.matrix(cbind(rep(1,length(y.new)),x.clin.new,x.gene.new)),y=y.new)
# dim 
lapply(coeff,function(x) sum(x!=0))
