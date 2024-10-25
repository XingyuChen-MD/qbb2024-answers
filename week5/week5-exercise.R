#BiocManager::install("vsn")

library(DESeq2)
library(vsn)
library(matrixStats)
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)

# Load the RNA-seq data file
RNAseq_data <- read_tsv("salmon.merged.gene_counts.tsv")

# Convert gene names into row names and remove the 'gene_id' column
RNAseq_data <- RNAseq_data %>%
  column_to_rownames(var = "gene_name") %>%
  select(-gene_id)

# Convert all numeric columns to integers
RNAseq_data <- RNAseq_data %>%
  mutate(across(where(is.numeric), as.integer))

# Filter out low coverage genes (keep genes with more than 100 total counts)
RNAseq_data <- RNAseq_data[rowSums(RNAseq_data) > 100, ]

# Select a subset of columns for analysis ('narrow' dataset)
narrow <- RNAseq_data %>% select("A1_Rep1":"P2-4_Rep3")

# Create metadata tibble for DESeq2 analysis
narrow_metadata <- tibble(
  tissue = factor(c("A1", "A1", "A1", "A2-3", "A2-3", "A2-3", "Cu", "Cu", "Cu", 
                    "LFC-Fe", "LFC-Fe", "Fe", "LFC-Fe", "Fe", "Fe", "P1", "P1", 
                    "P1", "P2-4", "P2-4", "P2-4")),
  rep = factor(rep(1:3, length.out = 21))
)

# Create DESeq2 data object
narrow_data <- DESeqDataSetFromMatrix(countData = as.matrix(narrow), 
                                      colData = narrow_metadata, 
                                      design = ~ tissue)

# Perform variance stabilizing transformation (VST) to correct batch effects
narrow_vst_data <- vst(narrow_data)

# Visualize the mean-variance relationship
meanSdPlot(assay(narrow_vst_data))

# Perform PCA and save PCA plot
narrow_pca_data <- plotPCA(narrow_vst_data, intgroup = c("rep", "tissue"), returnData = TRUE)
# Create the PCA plot
pca_plot <- ggplot(narrow_pca_data, aes(PC1, PC2, color = tissue, shape = rep)) +
  geom_point(size = 5) +
  labs(title = "PCA Plot of Midgut Tissues") +
  theme_minimal()

# Display the plot (optional, for visual inspection)
print(pca_plot)

# Save the plot to a file
ggsave("PCAPlot.png", plot = pca_plot, width = 8, height = 6, dpi = 300)

# Convert VST data into a matrix for downstream analysis
mat_narrow <- assay(narrow_vst_data)
# Calculate mean of replicates for each gene
replicate_means <- (mat_narrow[, seq(1, 21, 3)] +
                      mat_narrow[, seq(2, 21, 3)] +
                      mat_narrow[, seq(3, 21, 3)]) / 3

# Filter out low variance genes directly from the matrix
filtered_genes <- rowSds(replicate_means) > 1
mat_narrow_filtered <- mat_narrow[filtered_genes, ]
mat_narrow <- mat_narrow[filtered_genes, ]

# Perform k-means clustering (set seed for reproducibility)
set.seed(42)
k_clusters <- kmeans(mat_narrow, centers = 12)$cluster

# Order matrix by cluster
ordering <- order(k_clusters)
k_clusters <- k_clusters[ordering]
mat_narrow <- mat_narrow[ordering, ]

# Plot heatmap with clusters and save as an image
heatmap_colors <- brewer.pal(12, "Paired")[k_clusters]
png("heatmap.png", width = 800, height = 600)
heatmap(mat_narrow, Rowv = NA, Colv = NA, RowSideColors = heatmap_colors, scale = "row")
dev.off()

# Extract genes from the first cluster and save to a file
cluster_1_genes <- rownames(mat_narrow[k_clusters == 1, ])
write.table(cluster_1_genes, "cluster1.txt", sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)

cluster_1_genes <- rownames(mat_narrow[k_clusters == 1, ])
write.table(cluster_1_genes, "cluster1_genes.txt", sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
