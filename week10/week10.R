# Load libraries
library(ggplot2)

# Load data
data <- read.csv("Nuclei_data.csv")

# Plot nascent RNA signals
p1 <- ggplot(data, aes(x = gene, y = nascentRNA, fill = gene)) +
  geom_violin() +
  theme_minimal() +
  ggtitle("Nascent RNA Signal by Gene") +
  xlab("Gene") +
  ylab("Nascent RNA Signal")

# Plot PCNA signals
p2 <- ggplot(data, aes(x = gene, y = PCNA, fill = gene)) +
  geom_violin() +
  theme_minimal() +
  ggtitle("PCNA Signal by Gene") +
  xlab("Gene") +
  ylab("PCNA Signal")

# Plot log2 ratio of nascent RNA to PCNA
p3 <- ggplot(data, aes(x = gene, y = log2_ratio, fill = gene)) +
  geom_violin() +
  theme_minimal() +
  ggtitle("Log2 Ratio (Nascent RNA / PCNA) by Gene") +
  xlab("Gene") +
  ylab("Log2 Ratio")

# Arrange plots
library(gridExtra)
grid.arrange(p1, p2, p3, nrow = 3)
