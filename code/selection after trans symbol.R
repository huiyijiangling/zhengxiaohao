#第3步，转换为正常的genesymbol
set.seed(519)
cv_fit <- glmnet::cv.glmnet(x, y, alpha=1,type.measure="C",standardize = FALSE,family = 'cox',nfolds =3)
fit <- glmnet::glmnet(x=x, y=y, alpha = 1,family ='cox', lambda=cv_fit$lambda.min)#lambda.1se
choose_gene=rownames(fit$beta)[as.numeric(fit$beta)!=0]
library(plotmo)
pdf(file='lasso_total.pdf',height = 8,width = 8)
model_lasso <- glmnet(x, y, alpha=1,family ='cox')#
plot_glmnet(model_lasso,xvar = c("rlambda"), label=F,col=1:6) 
dev.off()
pdf(file='cross validation glmnet.pdf',height = 8,width = 8)
plot(cv_fit)
dev.off()
# coefsort=coef(cv_fit, s = 0)
# coefsort=as.data.frame(as.matrix(coefsort))
# colnames(coefsort)="coefsort"
# coefsort$coefsort=as.numeric(coefsort$coefsort)
# coefsort$ensg=rownames(coefsort)
# coefsort=coefsort[order(coefsort$coefsort,decreasing = T),]
# 
# coefsort$num=1:nrow(coefsort)
# col=rep(0,ncol(x))
# col[which(rownames(coefsort) %in% choose_gene)]=1:table(as.numeric(fit$beta)!=0)[2]

library(org.Hs.eg.db)
df=biomaRt::select(org.Hs.eg.db, keys = colnames(x), columns=c("ENSEMBL","SYMBOL"), keytype="ENSEMBL")
df[duplicated(df$ENSEMBL),]
df=subset(df,SYMBOL!="C9orf47")
df=df[df$ENSEMBL %in% choose_gene,]
if(F){
xx=x
colnames(xx)=df$SYMBOL
set.seed(1)
model_lasso <- glmnet(xx, y, alpha=1,family ='cox')#
pdf(file='lasso_select.pdf',height = 8,width = 8)
plot_glmnet(model_lasso,xvar = c("rlambda"),s=cv_fit$lambda.min,label=F,col=1:20) #颜色不对
dev.off()

xxx=xx[,which(rownames(coefsort) %in% choose_gene)]
}


#dplyr,plyrrename刚好相反
unlist(paste(df$SYMBOL,choose_gene,sep="=", collapse = ","))
# 


v0d=dplyr::rename(v0,ERRFI1=ENSG00000116285,TANC2=ENSG00000170921,ZNF189=ENSG00000136870,MLYCD=ENSG00000103150,TMEM80=ENSG00000177042,EIF4EBP1=ENSG00000187840,ITGA3=ENSG00000005884)
v2d=dplyr::rename(v2,ERRFI1=ENSG00000116285,TANC2=ENSG00000170921,ZNF189=ENSG00000136870,MLYCD=ENSG00000103150,TMEM80=ENSG00000177042,EIF4EBP1=ENSG00000187840,ITGA3=ENSG00000005884)
v3d=dplyr::rename(v3,ERRFI1=ENSG00000116285,TANC2=ENSG00000170921,ZNF189=ENSG00000136870,MLYCD=ENSG00000103150,TMEM80=ENSG00000177042,EIF4EBP1=ENSG00000187840,ITGA3=ENSG00000005884)
v4d=dplyr::rename(v4,ERRFI1=ENSG00000116285,TANC2=ENSG00000170921,ZNF189=ENSG00000136870,MLYCD=ENSG00000103150,TMEM80=ENSG00000177042,EIF4EBP1=ENSG00000187840,ITGA3=ENSG00000005884)
v6d=dplyr::rename(v6,ERRFI1=ENSG00000116285,TANC2=ENSG00000170921,ZNF189=ENSG00000136870,MLYCD=ENSG00000103150,TMEM80=ENSG00000177042,EIF4EBP1=ENSG00000187840,ITGA3=ENSG00000005884)
v7d=dplyr::rename(v7,ERRFI1=ENSG00000116285,TANC2=ENSG00000170921,ZNF189=ENSG00000136870,MLYCD=ENSG00000103150,TMEM80=ENSG00000177042,EIF4EBP1=ENSG00000187840,ITGA3=ENSG00000005884)
datd=dplyr::rename(dat,ERRFI1=ENSG00000116285,TANC2=ENSG00000170921,ZNF189=ENSG00000136870,MLYCD=ENSG00000103150,TMEM80=ENSG00000177042,EIF4EBP1=ENSG00000187840,ITGA3=ENSG00000005884)
save(df,v0d,v2d,v3d,v4d,v6d,datd,file="vd.Rdata")
load(file="vd.Rdata")

