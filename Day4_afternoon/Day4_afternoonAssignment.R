
# Step 1.1
library(ggplot2)
library(readr)
library(dplyr)

# Load the data into a tibble
dnm_data <- read_csv("aau1043_dnm.csv")
colnames(dnm_data)
head(dnm_data)

# Step 1.2 Use group_by() and summarize() to tabulate the number of paternally and maternally inherited DNMs in each proband. Note that the maternal versus paternal origin of the mutations are recorded in the column titled Phased_combined.
# Tabulate the number of paternally and maternally inherited DNMs in each proband
dnm_summary <- dnm_data %>%
  group_by(Proband_id, Phase_combined) %>%
  summarize(Number_of_DNMs = n(), .groups = 'drop') %>%
  pivot_wider(names_from = Phase_combined, values_from = Number_of_DNMs, names_prefix = "DNMs_")

# View the summary
print(dnm_summary)

#Step 1.3 Now, load the data from aau1043_parental_age.csv into a new dataframe.

aau1043_parental_age <- read_csv("aau1043_parental_age.csv")

# Step 1.4 You now have two dataframes with complementary information. It would be nice to have all of this in one dataframe. Use the left_join() function to combine your dataframe from step 2 with the dataframe you just created in step 3 based on the shared column Proband_id.
 
merged_data <- left_join(dnm_summary,aau1043_parental_age, by="Proband_id") 
merged_data


## Exercise 2: Fit and interpret linear regression models ###

#Step 2.1
#First, you’re interested in exploring if there’s a relationship between the number of DNMs and parental age. Use ggplot2 to plot the following. All plots should be clearly labelled and easily interpretable.
#1. the count of maternal de novo mutations vs. maternal age

ggplot(data=merged_data, mapping=aes(x=DNMs_mother,y=Mother_age))+
  geom_point()+
  stat_smooth(method="lm")


#2. the count of paternal de novo mutations vs. paternal age
ggplot(data=merged_data, mapping=aes(x=DNMs_father,y=Father_age))+
  geom_point()+
  stat_smooth(method="lm")


#Step 2.2 Now that you’ve visualized these relationships, you’re curious whether they’re statistically significant. Fit a linear regression model to the data using the lm() function.

### maternal
lm(data=merged_data, formula= DNMs_mother ~1, + Mother_age) %>% summary()

### Result:
Call:  lm(formula = DNMs_mother ~ 1, data = merged_data, subset = +Mother_age)
Residuals:
  Min      1Q  Median      3Q     Max 
-10.677  -4.677  -2.677   5.323  16.323 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept)  14.6768     0.3767   38.97   <2e-16 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 7.496 on 395 degrees of freedom

### Parental
lm(data=merged_data, formula= DNMs_father ~1, + Father_age) %>% summary()

### results
Call:
  lm(formula = DNMs_father ~ 1, data = merged_data, subset = +Father_age)

Residuals:
  Min      1Q  Median      3Q     Max 
-29.308 -12.308  -2.308  11.692  31.692 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept)  54.3081     0.7205   75.38   <2e-16 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 14.34 on 395 degrees of freedom

# Q: What is the “size” of this relationship? In your own words, what does this mean? Does this match what you observed in your plots in step 2.1?
Intercepts: The average number of maternally inherited DNMs is significantly lower than the average number of paternally inherited DNMs. Specifically, the number of maternally inherited DNMs is around 14.68, while the number of paternally inherited DNMs is around 54.31.
Residuals: The range and spread of residuals are larger for paternally inherited DNMs compared to maternally inherited DNMs, indicating more variability in the number of paternally inherited DNMs.

## Does this match what you observed in your plots in step 2.1?
A higher average number of paternally inherited DNMs compared to maternally inherited DNMs.

# Q:Is this relationship significant? How do you know? In your own words, what does this mean?
Yes, In both models, the intercepts have very low p-values (less than 2e-16), which indicates that the average number of maternally and paternally inherited DNMs is significantly different from zero. This is a strong indicator that the observed differences are not due to random chance.
  

##Step 2.4 Using your results from step 2.3, predict the number of paternal DNMs for a proband with a father who was 50.5 years old at the proband’s time of birth. Record your answer and your work (i.e. how you got to that answer).
maternal_model=lm(data=merged_data, formula= DNMs_mother ~1, + Mother_age) 
Parental_model=lm(data=merged_data, formula= DNMs_father ~1, + Father_age)

new_Parental<- tibble(Father_age= 50.5)
new_Parental
## Using predict function
predict(Parental_model,new_Parental)
54.30808 
## Re: the  predicted the number of paternal DNMs is 54.30808 


##Step 2.5 Next, you’re curious whether the number of paternally inherited DNMs match the number of maternally inherited DNMs. Plot the distribution of maternal DNMs per proband (as a histogram). In the same panel (i.e. the same set of axes) plot the distribution of paternal DNMs per proband. Make sure to make the histograms semi-transparent so you can see both distribution
ggplot(merged_data) +
  geom_histogram(aes(x = DNMs_mother, y = ..density..), binwidth = 1, fill = "blue", alpha = 0.5, position = "identity") +
  geom_histogram(aes(x = DNMs_father, y = ..density..), binwidth = 1, fill = "red", alpha = 0.5, position = "identity") +
  labs(title = "Distribution of Maternal and Paternal DNMs per Proband",
       x = "Number of DNMs",
       y = "Density") +
  theme_minimal()


#Step 2.6 Now that you’ve visualized this relationship, you want to test whether there is a significant difference between the number of maternally vs. paternally inherited DNMs per proband. What would be an appropriate statistical model to test this relationship? Fit this model to the data.
t_test_result <- t.test(merged_data$DNMs_mother, merged_data$DNMs_father, paired = TRUE)
# Display the result
print(t_test_result)
###
data:  merged_data$DNMs_mother and merged_data$DNMs_father
t = -61.609, df = 395, p-value < 2.2e-16
alternative hypothesis: true mean difference is not equal to 0
95 percent confidence interval:
  -40.48685 -37.98284
sample estimates:
  mean difference 
-39.23485 

# After performing your test, answer the following questions:
# What statistical test did you choose? Why?
Re: I chose t-test, because t-test is a statistical tool used to determine if_ there is a statistically significant difference between the means of two groups or between a group''s mean and a standard value.

#  Was your test result statistically significant? Interpret your result as it relates to the number of paternally and maternally inherited DNMs.
Re: Yes, it is significant because the p-value < 2.2e-16
