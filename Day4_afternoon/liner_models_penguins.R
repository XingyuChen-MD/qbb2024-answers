library(ggplot2)
library(ggsci)
library(palmerpenguins)
library(tidyverse)
library(ggthemes)
library(broom)

head(penguins)
glimpse(penguins)

ggplot(data=penguins, mapping=aes(x=bill_length_mm,y=bill_depth_mm, color=species))+
  geom_point()+
  stat_smooth(method="lm")

ggplot(data=penguins, mapping=aes(x=bill_length_mm,y=bill_depth_mm, color=sex))+
  geom_point()+
  stat_smooth(method="lm")

penguins1<-penguins %>% filter(!is.na(sex))

lm(data=penguins, formula= bill_length_mm ~1, + bill_depth_mm) %>% summary()

# Ensure species is a factor
penguins$species <- as.factor(penguins$species)

lm(data=penguins, formula= bill_length_mm ~1 + bill_depth_mm + species) %>% summary()

lm(data=penguins, formula= bill_length_mm ~1 + bill_depth_mm * species) %>% summary()

##sex
lm(data=penguins, formula= bill_length_mm ~1 + bill_depth_mm + species +sex) %>% summary()

##Does spcies matter for bill length

full_model <-lm(data=penguins, formula= bill_length_mm ~1 + bill_depth_mm+ species+ sex) 
summary(full_model)

reduced_model <- lm(data=penguins, formula= bill_length_mm ~1 + bill_depth_mm + species) 

anova(full_model,reduced_model)

13.2164 + 13.4033 + 1.3940 *17.5

#Gentoo/male/Penguin with a bill depth of 17.5 mm, what is the bill depth

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept)       27.6224     2.6718  10.339  < 2e-16 ***
  bill_depth_mm      0.5317     0.1513   3.514 0.000504 ***
  speciesChinstrap   9.9709     0.3357  29.699  < 2e-16 ***
  speciesGentoo     10.4890     0.5828  17.999  < 2e-16 ***
  sexmale            2.8938     0.3385   8.548 4.82e-16 ***
  ---
  
27.6224 + 0.5317*17.5 +10.4890 + 2.8938

new_penguin<- tibble(species="Gentoo",bill_depth_mm=17.5,sex="male")
new_penguin
## Using predict function
predict(full_model,new_penguin)
  