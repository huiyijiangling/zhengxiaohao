test=iris
x = ifelse(iris$Sepal.Width>3,"big","small")
table(x)
test$jjjd=ifelse(iris$Sepal.Width>3,"big","small")
library(dplyr)
test2 = mutate(test,x)
