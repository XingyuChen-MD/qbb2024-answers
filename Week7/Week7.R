
# 1.1.1: Load DESeq2 and tidyverse libraries, read data
library(DESeq2)
library(tidyverse)

# Load gene expression data and metadata
gene_data <- read_delim("gtex_whole_blood_counts_downsample.txt", delim = ",")
metadata <- read_delim("gtex_metadata_downsample.txt", delim = ",")

# 1.1.2: Move gene names to rownames in gene expression data
gene_data <- gene_data %>%
  column_to_rownames(var = "GENE_NAME")

# 1.1.3: Move subject IDs to rownames in metadata
metadata <- metadata %>%
  column_to_rownames(var = "SUBJECT_ID")

# 1.1.4: Peek at the first few rows of gene_data and metadata
head(gene_data)
head(metadata)

# Step 1.2: Create a DESeq2 object
# 1.2.1: Verify column order of count matrix and metadata row order
all(rownames(metadata) == colnames(gene_data)) # should return TRUE

# 1.2.2: Create DESeq2 object with formula including SEX, DTHHRDY, AGE
dds <- DESeqDataSetFromMatrix(countData = gene_data,
                              colData = metadata,
                              design = ~ SEX + DTHHRDY + AGE)

# Step 1.3: Normalization and PCA
# 1.3.1: Apply variance stabilizing transformation
vsd <- vst(dds)

# 1.3.2: Plot PCA with different intragroup values, save plots as PNG
# Generate PCA plot using top 500 features by variance, colored by "SEX"
# Generate PCA plot for "SEX" and save it
pca_plot_sex <- plotPCA(vsd, intgroup = "SEX", ntop = 500)
ggsave("pca_sex.png", plot = pca_plot_sex, width = 7.49, height = 8.17)

# Generate PCA plot for "DTHHRDY" and save it
pca_plot_dthhrdy <- plotPCA(vsd, intgroup = "DTHHRDY", ntop = 500)
ggsave("pca_dthhrdy.png", plot = pca_plot_dthhrdy, width = 7.49, height = 8.17)

# Generate PCA plot for "AGE" and save it
pca_plot_age <- plotPCA(vsd, intgroup = "AGE", ntop = 500)
ggsave("pca_age.png", plot = pca_plot_age, width = 7.49, height = 8.17)



# 1.3.3: PCA interpretation in comments
# PC1 and PC2 explain major variance. Observed that SEX aligns more with PC1, 
# suggesting sex-based gene expression differences. Age and DTHHRDY correlate 
# with PC2, indicating these might influence gene expression as well.

# Exercise 2: Perform differential expression analysis

# Step 2.1: Homemade test for differential expression
# Convert VST expression matrix to tibble with metadata
vsd_df <- assay(vsd) %>%
  t() %>%
  as_tibble() %>%
  bind_cols(metadata)

# 2.1.1: Test WASH7P gene for differential expression by SEX
# Load the broom package
library(broom)

# Run the linear model and tidy the output
m1 <- lm(WASH7P ~ DTHHRDY + AGE + SEX, data = vsd_df) %>%
  summary() %>%
  tidy()
m1

## A tibble: 4 × 5
term                   estimate std.error statistic  p.value
<chr>                     <dbl>     <dbl>     <dbl>    <dbl>
  1 (Intercept)              6.57      0.165     39.9   4.44e-64
2 DTHHRDYventilator_case  -0.0424    0.114     -0.372 7.11e- 1
3 AGE                      0.0587    0.0395     1.49  1.40e- 1
4 SEXmale                  0.119     0.109      1.09  2.79e- 1
# Interpretation:
# WASH7P shows not evidence of differential expression by SEX beacuse p-value is 0.279 > 0.05.
# the estimate for SEXMale > 0, then upregulated in males.

# 2.1.2: Test SLC25A47 gene for differential expression by SEX
m2 <- lm(SLC25A47 ~ DTHHRDY + AGE + SEX, data = vsd_df) %>%
  summary() %>%
  tidy()

m2
# A tibble: 4 × 5
term                   estimate std.error statistic  p.value
<chr>                     <dbl>     <dbl>     <dbl>    <dbl>
  1 (Intercept)              3.53      0.345      10.2  2.38e-17
2 DTHHRDYventilator_case  -0.760     0.239      -3.17 1.98e- 3
3 AGE                      0.0943    0.0827      1.14 2.57e- 1
4 SEXmale                  0.518     0.229       2.26 2.57e- 2
# Interpretation:
#Yes, it seems have significant effect on expression when the sample is male
# the estimate for SEXMale is 0.518 > 0, then upregulated in males.
#The p-value is less than p=0.05, making it statistically significant 


# Step 2.2: DE analysis with DESeq2
# 2.2.1: Run DESeq function on raw counts data (dds object)
dds <- DESeq(dds)

# Step 2.3: Extract and interpret sex differential expression results
# 2.3.1: Extract results for SEX
res_sex <- results(dds, name = "SEX_male_vs_female") %>%
  as_tibble(rownames = "GENE_NAME")
# 2.3.2: Filter for significant genes at 10% FDR
sig_genes_sex <- res_sex %>%
  filter(padj < 0.1)

# Count significant genes
num_sig_genes_sex <- nrow(sig_genes_sex)
num_sig_genes_sex # Display the count of significant genes

# 2.3.3: Load gene-to-chromosome mapping and merge with results
chromosome_map <- read_delim("gene_locations.txt")
# Check column names in chromosome_map
colnames(chromosome_map)

merged_res <- left_join(res_sex, chromosome_map, by = "GENE_NAME") %>% 
  arrange(padj)

# Interpretation:
# Observed X and Y chromosome associations for sex-differentiated genes, with Y-linked genes more upregulated in males and some X-linked in females.
#The top differential expression genes were exclusively found on sex chromosomes
#It matches with our prediction as we focused on the sex difference 

# 2.3.4: Check consistency with previous results for WASH7P and SLC25A47
res_sex %>% filter(GENE_NAME == "WASH7P" | GENE_NAME == "SLC25A47")
# A tibble: 2 × 7
GENE_NAME baseMean log2FoldChange lfcSE  stat        pvalue        padj
<chr>        <dbl>          <dbl> <dbl> <dbl>         <dbl>       <dbl>
  1 WASH7P       106.          0.0893 0.121 0.739 0.460         0.899      
2 SLC25A47      18.5         3.06   0.501 6.10  0.00000000105 0.000000832
# Exercise 2.4: Differential expression by death classification
# 2.4.1: Extract results for death classification
res_death <- results(dds, name = "DTHHRDY_ventilator_case_vs_fast_death_of_natural_causes") %>%
  as_tibble(rownames = "GENE_NAME")

# Filter significant genes for death classification
sig_genes_death <- res_death %>%
  filter(padj < 0.1)

# Count significant genes
num_sig_genes_death <- nrow(sig_genes_death)
num_sig_genes_death # Display the count of significant genes

# 2.4.2: Interpretation
# Death type may affect a wider gene set than SEX due to stress responses,
# affecting expression patterns linked to trauma or prolonged illness.

# Exercise 3: Visualization

# 3.1: Volcano plot for differential expression by SEX
res_sex <- res_sex %>%
  mutate(log_padj = -log10(padj))

# Define significance criteria
volcano_plot <- ggplot(res_sex, aes(x = log2FoldChange, y = log_padj)) +
  geom_point(aes(color = (padj < 0.1 & abs(log2FoldChange) > 1))) +
  theme_minimal() +
  labs(title = "Volcano Plot for Sex Differential Expression", 
       x = "log2 Fold Change", 
       y = "-log10 Adjusted P-value")

# Save the plot as PNG
ggsave("volcano_plot_sex.png", plot = volcano_plot)
