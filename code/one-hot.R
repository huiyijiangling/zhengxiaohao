# ONE-hot
set.seed(555)
data <- data.frame(
  Outcome = seq(1,100,by=1),
  Variable = sample(c("Red","Green","Blue"), 100, replace = TRUE)
)
#1
library(mltools)
library(data.table)
data$Variable <- as.factor(data$Variable)
newdata <- one_hot(as.data.table(data))
#2
library(caret)
dummy <- dummyVars(" ~ .", data=data)
newdata <- data.frame(predict(dummy, newdata = data)) 
#
library(reshape2)

newdata <- dcast(data = data, Outcome ~ Variable, length)

  
set.seed(555)
data <- data.frame(ID = seq(1,100,by=1),
                   Colour = sample(c("Red","Green","Blue"), 100, replace = TRUE),
                   Quality = sample(c("Poor","Average","Good"), 100, replace = TRUE)
)
newdata <- one_hot(as.data.table(data))
#
dummy <- dummyVars(" ~ .", data=data)
newdata <- data.frame(predict(dummy, newdata = data))
#
newdata <- dcast(data = melt(data, id.vars = "ID"), ID ~ variable + value, length)
