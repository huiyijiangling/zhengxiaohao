#注意计算时我直接使用RFI置换了dfi
library(readxl)
# phe_gdc=read.csv(file="PAAD_ductal_clinical_gdc.csv)
phe_cellPAAD=read_excel("leisile.xlsx",sheet = 2,col_types="numeric")
phe_cellPAAD_RFI01=phe_cellPAAD[phe_cellPAAD$RFI %in% c(0,1),]
dim(phe_cellPAAD_RFI01)
# #注意一定是用rnaExpr(已经配平的)，而不是用raw_count,raw count 没有校正。
# phe=metaMatrix.RNA[substr(rownames(metaMatrix.RNA),14,15)=="01",]
# phe$event=ifelse(phe$vital_status=='Alive',0,1)
# phe$time=as.numeric(ifelse(phe$vital_status=='Alive',phe$days_to_last_follow_up,phe$days_to_death))
# #4 删除随访NA 
# dim(phe)
# phe=subset(phe,time!="NA")
# #4 删除随访太少30，但是另一个研究直接写的1，怎么办滤过不掉啊
# dim(phe)
# phe$time=phe$time/30
# phe=phe[phe$time>=1,]
# dim(phe)
# phe=phe[-which(phe$time<36 & phe$event=="0"),]#不要随便删病例，对制造模型没有好处
phe=phe_cellPAAD_RFI01
#################create the input Expr tumor only 

# group_list=ifelse(substr(colnames(rnaCounts_Pancreas),1,4)=="GTEX",'normal',
#                   ifelse(substr(colnames(rnaCounts_Pancreas),1,4)=="TCGA" & substr(colnames(rnaCounts_Pancreas),14,15)=="11","normal",
#                          ifelse(substr(colnames(rnaCounts_Pancreas),1,4)=="TCGA" & substr(colnames(rnaCounts_Pancreas),14,15)=="01","tumor","wrong")))
# table(group_list)
# exprSet=rnaCounts_Pancreas[,group_list=='tumor']

# exprSet_rownames <- rownames(exprSet)
# exprSet_colnames <- colnames(exprSet)
# exprSet_quant <- preprocessCore::normalize.quantiles(
#   as.matrix(exprSet))
# rownames(exprSet_quant) <- exprSet_rownames
# colnames(exprSet_quant) <- exprSet_colnames
# exprSet_quant=log2(exprSet_quant+1)

# if(T){
#   load("ans.Rdata")
#   library(stringr)
#   lnchighall=unique(str_split(ansT3,"_",simplify = T)[,1])#ansALL
#   mhighall=unique(str_split(ansT3,"_",simplify = T)[,2])#ansALL
#   # ansT
#   # ansN
# }
#中心化前去掉不分析的样本,col_types="numeric"
exprSet_quant=read_excel("leisile.xlsx",sheet = 3,col_types=c("text",rep("numeric",247)))
rn=exprSet_quant[,1]
exprSet_quant=exprSet_quant[,-1]
exprSet_quant=as.data.frame(exprSet_quant)
rownames(exprSet_quant)=rn$bah

exprSet_quant_with_clinical=exprSet_quant[,which(colnames(exprSet_quant) %in% phe$bah)]
phe=phe[which(phe$bah %in% colnames(exprSet_quant_with_clinical)),]
phe=phe[order(phe$bah,decreasing = F),]
exprSet_quant_with_clinical=exprSet_quant_with_clinical[,order(colnames(exprSet_quant_with_clinical),decreasing = F)]

dat$event=as.numeric(dat$RFI)#dat$event=as.numeric(dat$event)
dat$time=as.numeric(dat$RFI.time)#dat$time=as.numeric(dat$time)
dat$sample_type="PrimaryTumor"
library(GDCRNATools)
survOutput <- gdcSurvivalAnalysis(gene     = rownames(exprSet_quant_with_clinical), 
                                  method   = 'coxph', 
                                  rna.expr = exprSet_quant_with_clinical, 
                                  metadata = dat)
survOutput <- gdcSurvivalAnalysis(gene     = mrna_2list,
                                  method   = 'KM',
                                  rna.expr = rnaExpr,
                                  metadata = phe,
                                  sep      = 'median')
XXXXXXXXXXXXXXXXx=rownames(survOutput[survOutput$pValue<0.05,])

#建立训练集和验证集
# library(caret)
# set.seed(1)
# sam<- createDataPartition(phe$event, p = 0.5,list = FALSE)
# dim(sam)
# t_exp=exprSet_quant_with_clinical[,sam] ##
# v_exp=exprSet_quant_with_clinical[,-sam]## 
# 
# table(ifelse(substr(colnames(t_exp),14,15)=='01','tumor','normal'))
# table(ifelse(substr(colnames(v_exp),14,15)=='01','tumor','normal'))
# 
# t_tumor=t_exp[,substr(colnames(t_exp),14,15)=='01']
# v_tumor=v_exp[,substr(colnames(v_exp),14,15)=='01']
# 
# t_phe=phe[match(substr(colnames(t_tumor),1,12),phe$patient),]
# v_phe=phe[match(substr(colnames(v_tumor),1,12),phe$patient),]
# #查看两组一些临床参数切割比例
# table(t_phe$tumor_stage)
# table(v_phe$tumor_stage)

#这一步是中心化 要不要删除等待验证
if(F){
# library(DescTools)
# tcga_only_expr_mad=t(RobScale(t(exprSet_quant_with_clinical)))
# tcga_only_expr_mean=t(scale(t(exprSet_quant_with_clinical)))
# save(tcga_only_expr_mean,tcga_only_expr_mad,exprSet_quant,file="tcga_only_expr_after.Rdata")
}
# normalize（此步骤可省略，因为glmnet默认会标准化后建模，再返回变换后的真实系数）
# pp = preProcess(x,method = c("center", "scale"))
# x <- predict(pp, x)
## 必须保证生存资料和表达矩阵，两者一致

library(lars) 
library(glmnet) 
library(survival)
library(glmnet)

table(is.na(phe))
# all(substring(colnames(exprSet_quant_with_clinical),1,12)==phe$patient)
x=as.matrix(t(exprSet_quant_with_clinical))
y=phe$event
y <- Surv(phe$RFI.time, phe$RFI)#Surv(phe$time, phe$event)
# fit the model
set.seed(4)
model_lasso <- glmnet(x, y, alpha=1,family ='cox', nlambda=1000)#type.measure="auc",注意auc只能用来算二分类logistic#nlambda=10倍sample
print(model_lasso)
# 列%Dev代表了由模型解释的残差的比例，对于线性模型来说就是模型拟合的R^2(R-squred)。
# 它在0和1之间，越接近1说明模型的表现越好，
# 如果是0，说明模型的预测结果还不如直接把因变量的均值作为预测值来的有效。
# 使用area under the ROC curve, CV 选择压缩参数lambda
# 再设置一次set.seed
head(coef(model_lasso, s=c(model_lasso$lambda[282],0.24690)))#??????
#plot.glmnet(model_lasso, xvar = "norm", label = TRUE)
plot(model_lasso, xvar="lambda", label=TRUE)#???????
set.seed(85)
cv_fit <- cv.glmnet(x, y, alpha=1,family = 'cox',nfolds =50,type.measure="C")#默认的cox的dev使用的是likelihood，但可能c能有更好的结果,type.measure="C",type.measure="C", 
plot(cv_fit)
coef.min = coef(cv_fit, s = "lambda.min")
active.min = which(coef.min != 0)
index.min = coef.min[active.min]
index.min
# 两条虚线分别指示了两个特殊的λ值:
c(log(cv_fit$lambda.min),log(cv_fit$lambda.1se))
c(cv_fit$lambda.min,cv_fit$lambda.1se)

# cv_fit$cvm
fit <- glmnet(x=x, y=y, alpha = 1,family ='cox', lambda=cv_fit$lambda.min)#lambda.1se
head(fit$beta)
#一倍SE内的更简洁的模型,是22个miRNAmin
#fit <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
#head(fit$beta)# 这里是40个miRNA
choRFIe_gene=rownames(fit$beta)[as.numeric(fit$beta)!=0]
length(choRFIe_gene)
myexpr=x[,choRFIe_gene]
mysurv=phe[,c("RFI.time","RFI")]

save(choRFIe_gene,myexpr,mysurv,file="RFIdadada.Rdata")

jjju=DEGAll[choRFIe_gene,]
paste(jjju$symbol,collapse = ", ")
paste(choRFIe_gene,collapse = "+")
paste(names(exprSet4),collapse = " ")


########################################################################################################

model_lasso <- glmnet(x=x, y=y, alpha = 1,family = 'cox', lambda=cv_fit$lambda.1min)
lasso.prob <- predict(cv_fit, newx=x,s=c(cv_fit$lambda.min,cv_fit$lambda.1se))#lasso.prob <- predict(cv_fit, newx=t(v_tumor),s=c(cv_fit$lambda.min,cv_fit$lambda.1se))
re=cbind(y ,lasso.prob)#re=cbind(v_phe$event ,lasso.prob)
dat_0=as.data.frame(re[,1:2])
colnames(dat_0)=c('event','prob')
dat_0$event=as.factor(dat_0$event)
# #一个很好的boxplot
# library(ggpubr) 
# p <- ggboxplot(dat_0, x = "event", y = "prob",
#                color = "event", palette = "jco",
#                add = "jitter")
# #  Add p-value
# p + stat_compare_means()
fit <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
head(fit$beta)
## https://vip.biotrainee.com/d/812-
library(ROCR)
library(glmnet)
library(caret)
# calculate probabilities for TPR/FPR for predictions
pred <- prediction(re[,2], re[,1])
perf <- performance(pred,"tpr","fpr")
performance(pred,"auc") # shows calculated AUC for model
plot(perf,colorize=FALSE, col="black") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

if(F){
  model_lasso <- glmnet(x=x, y=y, alpha = 1,family ='cox', lambda=cv_fit$lambda.min)#lambda.1se
  head(model_lasso$beta)
  lasso.prob <- predict(cv_fit, newx=x , s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
  re=cbind(y ,lasso.prob)
  dat_0=as.data.frame(re[,1:2])
  colnames(dat_0)=c('event','prob')
  dat_0$event=as.factor(dat_0$event)
}
