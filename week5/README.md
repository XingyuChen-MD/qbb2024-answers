Exercise 1: Checking FASTQ Quality
# Step 1.1: Initial Quality Assessment with FastQC
Observations: No significant issues were observed across most quality metrics; however, the per-base sequence quality and per-sequence GC content metrics showed irregularities at the beginning of the sequences.
Explanation: The poor quality at the beginning of the sequences may be due to the limitations of sequencing technology. Often, the initial bases can have lower quality because the signal at the start of sequencing can be unstable, leading to errors in base calling. This is a common issue and can also arise if adapter sequences or other technical artifacts are present.
# Step 1.2: Overrepresented Sequences
Most Overrepresented Sequence: The most overrepresented sequence was identified as coding for serine proteases in Drosophila.
Interpretation: This finding makes sense because serine proteases play a crucial role in the digestive system, breaking down proteins during digestion. Since the samples were taken from Drosophila midgut tissues, it is logical to see a high abundance of these enzymes, reflecting their biological function in the gut.
# Exercise 2: Using MultiQC to Check Processed Data Quality
Sample Rejection Criteria: Based on the requirement to retain samples with unique read percentages ≥ 45%, 25 samples were found to be acceptable for downstream analysis.
Replicate Consistency: The DESeq2 sample-to-sample distance heatmap showed clear clustering of triplicates, indicating a high level of consistency between the replicates.
Adjusting the Minimum Slider: Increasing the minimum threshold slider on the heatmap can improve the clarity of the clusters, further confirming the consistency. This suggests that the experimental replicates were prepared and processed under consistent conditions, which is essential for reliable RNA-seq analysis.
# Exercise 3: Normalization and Clustering
## Step 3.3: PCA Analysis
PCA Plot Observations: There was a notable discrepancy, where A1 clustering with A2-3,  with Fe clustering with LFC-Fe samples and LFC-Fe clustering with Fe samples
Possible Solutions:
Revisit Metadata: Double-check the sample labels in the metadata to ensure they match the actual biological samples. Correcting any labeling errors could resolve the issue.
Reordering Using Matrix Indexing: If the problem persists after verifying the metadata, manually adjust the data matrix to reorder these samples and rerun the PCA.
## Step 3.6: Gene Ontology Enrichment Analysis
GO Enrichment Results:
Molecular Functions (MF):
Peptidase activity, catalytic activity, metallopeptidase activity, endopeptidase activity, hydrolase activity
Monooxygenase activity, serine-type endopeptidase activity, heme binding, iron ion binding
Oxidoreductase activity, serine hydrolase activity, acid sphingomyelin phosphodiesterase activity
Biological Processes (BP):
Proteolysis, lipid metabolic process
Interpretation: These enrichment categories are consistent with the known biological function of midgut tissues. The high representation of terms related to metabolic pathways and proteolysis aligns with the gut's role in digestion and nutrient absorption. For instance:
Peptidase and proteolysis activities reflect the breakdown of dietary proteins into smaller peptides and amino acids.
Lipid metabolic processes are essential for absorbing and processing dietary fats.
Enrichments in heme binding and iron ion binding can indicate the involvement of enzymes related to redox reactions, important for metabolic functions within the gut cells.