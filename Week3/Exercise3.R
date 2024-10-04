library(tidyverse)

#Load the data
AF <- read_tsv("/Users/cmdb/qbb2024-answers/week3-all/AF.txt")
DP <- read_tsv("/Users/cmdb/qbb2024-answers/week3-all/DP.txt")

#Create a histogram with the ggplot 
AF_plot <- ggplot(AF, aes(x = Allele_Frequency)) + 
  geom_histogram(bins = 11) + 
  labs(x = "Allele Frequency", y = "Number of Variation", title = "Allele Frequency")
ggsave("/Users/cmdb/qbb2024-answers/week3-all/AF_plot.png", plot = AF_plot)

#Create a histogram with the ggplot
DP_plot <- ggplot(DP, mapping = aes(x = Read_Depth)) + 
  geom_histogram(bins = 21) + 
  xlim(0, 20) + 
  labs(x = "Read Depth", y = "Number of Variation", title = "Read Depth")
ggsave("/Users/cmdb/qbb2024-answers/week3-all/DP_plot.png", plot = DP_plot)
