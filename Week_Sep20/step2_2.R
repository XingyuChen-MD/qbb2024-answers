library(tidyverse)
library(ggthemes)

file <- read.delim("~/qbb2024-answers/9_20/snp_counts.txt") 

# Log-transform enrichment values
snpenrich <- file %>% mutate(logEnrichment = log2(Enrichment)) 

# Filter data by feature type
exons <- snpenrich %>% filter(Feature == "exons")
introns <- snpenrich %>% filter(Feature == "introns")
cCREs <- snpenrich %>% filter(Feature == "cCREs")
other <- snpenrich %>% filter(Feature == "other")

nature_palette <- c("Exons" = "#E69F00", "Introns" = "#56B4E9", "cCREs" = "#009E73", "Other" = "#F0E442")

p <- ggplot() +
  geom_line(data = exons, aes(MAF, logEnrichment, color = "Exons"), size = 1.2) +
  geom_line(data = introns, aes(MAF, logEnrichment, color = "Introns"), size = 1.2) +
  geom_line(data = cCREs, aes(MAF, logEnrichment, color = "cCREs"), size = 1.2) +
  geom_line(data = other, aes(MAF, logEnrichment, color = "Other"), size = 1.2) +
  
  # Add axis labels and title
  labs(
    color = "Genomic Feature",
    x = "Minor Allele Frequency (MAF)",
    y = expression(Log[2]~"SNP Enrichment"),
    title = "SNP Enrichment across Genomic Features by MAF"
  ) +
  
  # Customize the color palette
  scale_color_manual(values = nature_palette) +
  
  # Use a clean, publication-ready theme
  theme_classic() +
  
  # Refine the plot text sizes and style
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "top"
  )

ggsave(filename = "~/qbb2024-answers/9_20/snp_enrichments.pdf", plot = p, width = 8, height = 6)
