#1.构建模型。【节点数为3.4.5均进行了构建和比较】
library(rms)
fit1<- cph(Surv(time,event)~rcs(LODDS,3),data=dat);anova(fit1);fit1
fit2<- cph(Surv(time,event)~rcs(LODDS,4),data=dat);anova(fit2);fit2
fit3<- cph(Surv(time,event)~rcs(LODDS,5),data=dat);anova(fit3);fit3
# 根据Harrell建议，本数据为n=350，节点选择应为5。
# 
# anova( )函数中，P-Nonlinear<0.05为存在非线性关系； 
# 
# R²和Dxy值越大模型越优
# 
# 结果显示：节点选择为3~5时，拟合后的函数P-Nonlinear均<0.05，但节点取5，R²和Dxy最大，拟合的模型更优，与上述建议一致
# 
# 图片
# 
# 图片
# 

nomo<-datadist(dat)
options(datadist='nomo')
#2.计算风险比HR值与自变量(age)变化关系
HR<-Predict(fit3,LODDS,fun=exp,ref.zero=TRUE);HR
#3.可视化上述关系
P<-ggplot(HR);P



结果解释：基于图表结果

年龄<48岁， 复发风险随年龄增加而减低（当age=44.5时，HR≈1）；

年龄>48岁后，复发风险随年龄变化不是很明显。

综上，可以发现该病理类型乳腺癌患者的复发风险并不是随年龄呈线性变化的，48岁随后年龄对复发的影响趋于稳定且影响度很小。这也是符合实际意义。

#4.图形美化
P2<-ggplot()+
  geom_line(data=HR, #数据来源
            aes(age,yhat),#xy轴的数据
            linetype="solid",#曲线加粗
            size=1,
            alpha=0.7,
            colour="purple")+
  geom_ribbon(data=HR, #加入置信区间
              aes(age, 
                  ymin=lower, 
                  ymax=upper),
              alpha=0.1,
              fill="purple")+
  theme_classic()+
  geom_hline(yintercept=1,linetype=2,size=0.75)+ #y=1的水平线 
  labs(title ="复发风险随年龄的变化曲线：基于限制性立方样条(RCS)",
       x="年龄 (即连续变量)", 
       y="HR (95%CI)"
  );P2

fit4<- cph(Surv(time,event)~rcs(age,4)+ssfsy+NerveF,data=dat);anova(fit4);fit4
HR<-Predict(fit4,age=c(10,20,30,40,50,60,70),fun=exp,ref.zero=TRUE);HR
summary(fit4,data=dat)
with(new_dat,rcorr.cens(fp,Surv(time, event)  ))

?pol(glasgowComaScore, 2)
pol(pmin(dat$age,10),3)
