library(tidyverse)


head(iris)


iris_setosa <- iris %>%  filter(Species == "setosa")

ggplot(data=iris_setosa, mapping=aes(x=Sepal.Length,y=Sepal.Width))+
  geom_point()+
  stat_smooth(method="lm")+#confidence interval:se =TRUE
  xlim(4,4.5)+
  ylim(1,4.5)
?stat_smooth

#In the context of the R lm (linear model) function, ~1 represents a model with only an intercept term and no predictor variables. Let's break down what this means:
ml<-lm(data=iris_setosa, formula=Sepal.Length ~1 +Sepal.Width)

summary(ml)
