library(ggplot2)
library(ggsci)
library(palmerpenguins)
library(tidyverse)
library(ggthemes)

setwd("~/Data/GTEx")

#Q1. Load the tidyverse package, and use the function read_delim() to read in sample-level metadata that was obtained from the GTEx Portal (GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt). In addition, open the data dictionary GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx in Excel, which provides a description of each column in the .txt file.
data<-read_delim("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

#Q2. View the first rows of the tibble by simply entering the variable name in which you stored it. Notice that some of the columns were cut off due to the limits of the display. Use the glimpse() function to examine the data types and first entries of all of the columns.
library(dplyr)
data
data[1,]
glimpse(data)

#Q3. Use the filter() function to subset the dataset to only the RNA-seq data by selecting rows for which the SMGEBTCHT column contains the value "TruSeq.v1". TruSeq is a library preparation kit from Illumina.
# Filter the data to include only rows with SMGEBTCHT equal to "TruSeq.v1"
filtered_data <- data %>%
  filter(SMGEBTCHT == "TruSeq.v1")
# View the first few rows of the filtered dataset
head(filtered_data)


#Q4. Plot the number of samples from each tissue (SMTSD) as a barplot. (Hint: if you do not specify a y-axis, ggplot will use stat = count as the default, so the y-axis will represent the number of occurrences of each value of your x variable). See this webpage for a code snippet for rotating axis labels, which will be relevant throughout this exercise. Always be sure to label your axes with informative names!
library(ggplot2)
# Create the bar plot
ggplot(filtered_data, aes(x = SMTSD)) +
  geom_bar() +
  labs(x = "Tissue", y = "Number of Samples", title = "Number of Samples from Each Tissue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#Q5. The RNA integrity number is a measurement of the degree of RNA degradation based on characteristics of an electropherogram trace. It ranges from 1 to 10, with 10 being the least degraded. Plot the distribution of RNA integrity numbers across your samples. What type of plot is best for visualizing a single continuous distribution? Take a look at this “cheat sheet” for hints.
colnames(filtered_data)
# Histogram
ggplot(filtered_data, aes(x = SMRIN)) +
  geom_histogram(binwidth = 0.5, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(x = "RNA Integrity Number (RIN)", y = "Frequency", title = "Distribution of RNA Integrity Numbers") +
  theme_minimal(base_size = 12) +
  theme(
    title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

#Density Plot
ggplot(filtered_data, aes(x = SMRIN)) +
  geom_density(fill = "steelblue", color = "black", alpha = 0.5) +
  labs(x = "RNA Integrity Number (RIN)", y = "Density", title = "Density Plot of RNA Integrity Numbers") +
  theme_minimal(base_size = 12) +
  theme(
    title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
#What is the shape of the distribution? Is it unimodal?

#Q6. Copy your code from above, but now plot the distribution of RIN, stratified by tissue. Consider what type of plot is best for contrasting continuous distributions across multiple groups.

# Faceted Density Plot
ggplot(filtered_data, aes(x = SMRIN, fill = SMTSD)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ SMTSD) +
  labs(x = "RNA Integrity Number (RIN)", y = "Density", title = "Distribution of RNA Integrity Numbers by Tissue") +
  theme_minimal(base_size = 12) +
  theme(
    title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggplot(filtered_data, aes(x = SMTSD, y = SMRIN)) +
  labs(x = "SMTSD", y = "SMRIN", title = "Distribution of RNA Integrity Numbers by Tissue")+
  geom_boxplot()  +
  labs(x = "Tissue", y = "RNA Integrity Number (RIN)", title = "Number of RNA Integrity Number (RIN) from Each Tissue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Do you notice any differences across tissues? Are certain tissues outliers? What are your hypotheses to explain these observations?
#Re:Yes, there are difference


#Q7. Visualize the number of genes detected per sample, stratifying by tissue. Again consider what type of plot is best for contrasting continuous distributions across multiple groups.

if (!require(devtools)) {
  install.packages('devtools')
}
devtools::install_github('erocoar/gghalves')

library(gghalves)
ggplot(filtered_data, aes(x = SMTSD, y = SMGNSDTC)) +
  geom_half_boxplot() + 
  geom_dotplot(binaxis = "y", method="histodot", stackdir="up",dotsize = 0.01) +
  labs(x = "Tissue", y = "SMGNSDTC", title = "Number of SMGNSDTC from Each Tissue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
?geom_dotplot

#ANOTHER WAY
devtools::install_github("eclarke/ggbeeswarm")
library(ggbeeswarm)
ggplot(filtered_data,aes(SMTSD, SMGNSDTC)) + geom_quasirandom(size = 0.05)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("SMTSD + SMGNSDTC")

#Q8. Plot the relationship between ischemic time and RIN. Consider what type of plot is best for visualizing the relationship between two continuous variables. 
#Create sub-panels that stratify the data by tissue using facet_wrap(). Resize points to size = 0.5 and set the opacity to alpha = 0.5. Add linear trend lines to your plot with `geom_smooth(method = “lm”).
ggplot(data=filtered_data)+
  geom_point(mapping=aes(x = SMTSISCH, 
                         y = SMRIN, 
                         color =  SMTSD#, shape=SMTSD
  )) +
  scale_fill_jama() + 
  geom_smooth(methods=lm, mapping=aes(x = SMTSISCH, y = SMRIN))

ggplot(data=filtered_data)+
  geom_point(mapping=aes(x = SMTSISCH, 
                         y = SMRIN, 
                         size = 0.5,alpha = 0.5,#Resize points to size = 0.5 and set the opacity to alpha = 0.5. 
                         color =  SMTSD#, shape=SMTSD
                         )) +
  scale_fill_jama() + 
  facet_wrap(vars(SMTSD))+
  geom_smooth(methods=lm, mapping=aes(x = SMTSISCH, y = SMRIN))+#Add linear trend lines to your plot with `geom_smooth(method = “lm”).
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#What relationships do you notice? Does the relationship depend on tissue?
#Re:most of them are nagetively correlated in several tissues,  the relationship does not depend on tissue

#Q9. Copy your answer from question 6 above, but modify it to color your points by autolysis score (SMATSSCR). Note that if we place aes(color = SMATSSCR) within the ggplot() portion of the code, it will attempt to apply this mapping to all geom_s, including geom_smooth. To avoid this, place aes(color = SMATSSCR) within the geom_point() portion of the code.
#What relationships do you notice? Does the relationship depend on tissue?
#Re: We found a positive correlation between SMTSISCH and SMATSSCR and indeed saw that this relationship depended on the organization
ggplot(data=filtered_data)+
  geom_point(mapping=aes(x = SMTSISCH, 
                         y = SMRIN, color = SMATSSCR,
                         size = 0.5,alpha = 0.5,#Resize points to size = 0.5 and set the opacity to alpha = 0.5. 
  )) +
  scale_fill_jama() + 
  facet_wrap(vars(SMTSD))+
  geom_smooth(methods=lm, mapping=aes(x = SMTSISCH, y = SMRIN))+#Add linear trend lines to your plot with `geom_smooth(method = “lm”).
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Q10. If you finished early, make some more plots! What else can you learn about these data? 
#Consider which type of plot visualizes a particular relationship most effectively. 
#Keep it simple. Ideally, each figure should convey one main point. The purpose of a figure is to convey that point to the audience as clearly and concisely as possible.

# Install library(ggpubr)
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library(ggpubr)
p <- ggboxplot(filtered_data, x = "SMTSD", y = "SMRIN",
               color = "SMTSD",
               add = "jitter", shape = "SMTSD")  + scale_fill_jama()+#Add linear trend lines to your plot with `geom_smooth(method = “lm”).
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

ggviolin(filtered_data, x = "SMTSD", y = "SMRIN", fill = "SMTSD",
         add = "boxplot", add.params = list(fill = "white"))+ scale_fill_jama()+#Add linear trend lines to your plot with `geom_smooth(method = “lm”).
  theme(axis.text.x = element_text(angle = 90, hjust = 1))