library(ROSE)
data(hacide)
str(hacide.train)
table(hacide.train$cls)
prop.table(table(hacide.train$cls))           

library(rpart)
treeimb <- rpart(cls ~ ., data = hacide.train)
pred.treeimb <- predict(treeimb, newdata = hacide.test)
accuracy.meas(hacide.test$cls, pred.treeimb[,2])

roc.curve(hacide.test$cls, pred.treeimb[,2], plotit = F)

# 过采样
data_balanced_over <- ovun.sample(cls ~ ., data = hacide.train, method = "over",N = 1960)$data

table(data_balanced_over$cls)
#欠采样
data_balanced_under <- ovun.sample(cls ~ ., data = hacide.train, method = "under", N = 40, seed = 1)$data
table(data_balanced_under$cls)
# 双

data_balanced_both <- ovun.sample(cls ~ ., data = hacide.train, method = "both", p=0.5, N=1000, seed = 1)$data

table(data_balanced_both$cls)
# 人工
data.rose <- ROSE(cls ~ ., data = hacide.train, seed = 42)$data

table(data.rose$cls)

# 训练决策树
tree.rose <- rpart(cls ~ ., data = data.rose)
tree.over <- rpart(cls ~ ., data = data_balanced_over)
tree.under <- rpart(cls ~ ., data = data_balanced_under)
tree.both <- rpart(cls ~ ., data = data_balanced_both)

# 在测试集上做预测
pred.tree.rose <- predict(tree.rose, newdata = hacide.test)
pred.tree.over <- predict(tree.over, newdata = hacide.test)
pred.tree.under <- predict(tree.under, newdata = hacide.test)
pred.tree.both <- predict(tree.both, newdata = hacide.test)

# 人工数据合成AUC值
roc.curve(hacide.test$cls, pred.tree.rose[,2])
Area under the curve (AUC): 0.989
# 过采样AUC值
roc.curve(hacide.test$cls, pred.tree.over[,2])
Area under the curve (AUC): 0.798
# 欠采样AUC值
roc.curve(hacide.test$cls, pred.tree.under[,2])
Area under the curve (AUC): 0.867
# 双采样AUC值
roc.curve(hacide.test$cls, pred.tree.both[,2])
Area under the curve (AUC): 0.798

#这个包为我们提供了一些基于holdout和bagging的模型评估方法，这有助于我们判断预测结果是否有太大的方差。


ROSE.holdout <- ROSE.eval(cls ~ ., data = hacide.train, learner = rpart, method.assess = "holdout", extr.pred = function(obj)obj[,2], seed = 1)

ROSE.holdout
# 
# bagging: 随机森林
# boosting: AdaBoost, GBDT, XGBoost,LightGBM, CatBoost

