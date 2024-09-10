#!/usr/bin/env python3

# Question1
import sys
import numpy

### gunzip GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
### less GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct


# import sys
# # open File
# fs = open(sys.argv[1],mode = 'r')
# # skip 2 lines
# fs.readline()
# fs.readline()
# # split colume by tabs and skips two entries
# line = fs.readline()
# fields = line.strip('\n').split('\t')
# # create way to hold gene names
# tissues = fields[2:]
# gene_names = []
# gene_IDs = []
# expression = []
# # for each line
# for line in fs:
#     fields = line.strip('\n').split('\t')     
#     gene_IDs.append(fields[0])                 
#     gene_names.append(fields[1])               
#     expression.append(fields[2:])             
# fs.close()
# print(expression)

file_path = '/Users/Shared/Data/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'

# Initialize lists to store the data
tissue_names = []
gene_ids = []
gene_names = []
expression_data = []

# Open the file
with open(file_path, 'r') as fs:
    # Skip the first 2   lines of the file (header lines)
    for _ in range(2):
        fs.readline()
    
    # Read the rest of the file line by line
    for line in fs:
        # Split the line by tab character
        parts = line.strip().split('\t')
        
        # Check if the line has the expected number of parts
        if len(parts) < 3:
            continue
        
        # The first part is gene ID, second is gene name, the rest are expression levels
        gene_id = parts[0]
        gene_name = parts[1]
        expressions = parts[2:]
        
        # Append gene ID and name to their respective lists
        gene_ids.append(gene_id)
        gene_names.append(gene_name)
        
        # Process expression levels
        expression_data.append(expressions)
    
    # Extract tissue names from the first line of the expression data
    tissue_names = expression_data[0][1:]  # Skip the first column which is for gene IDs

# Print the results (for debugging purposes)
print("Tissue Names:", tissue_names)
# Tissue Names: ['Adipose - Visceral (Omentum)', 'Adrenal Gland', 'Artery - Aorta', 'Artery - Coronary', 'Artery - Tibial', 'Bladder', 'Brain - Amygdala', 'Brain - Anterior cingulate cortex (BA24)', 'Brain - Caudate (basal ganglia)', 'Brain - Cerebellar Hemisphere', 'Brain - Cerebellum', 'Brain - Cortex', 'Brain - Frontal Cortex (BA9)', 'Brain - Hippocampus', 'Brain - Hypothalamus', 'Brain - Nucleus accumbens (basal ganglia)', 'Brain - Putamen (basal ganglia)', 'Brain - Spinal cord (cervical c-1)', 'Brain - Substantia nigra', 'Breast - Mammary Tissue', 'Cells - Cultured fibroblasts', 'Cells - EBV-transformed lymphocytes', 'Cervix - Ectocervix', 'Cervix - Endocervix', 'Colon - Sigmoid', 'Colon - Transverse', 'Esophagus - Gastroesophageal Junction', 'Esophagus - Mucosa', 'Esophagus - Muscularis', 'Fallopian Tube', 'Heart - Atrial Appendage', 'Heart - Left Ventricle', 'Kidney - Cortex', 'Kidney - Medulla', 'Liver', 'Lung', 'Minor Salivary Gland', 'Muscle - Skeletal', 'Nerve - Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin - Not Sun Exposed (Suprapubic)', 'Skin - Sun Exposed (Lower leg)', 'Small Intestine - Terminal Ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole Blood']
print("Gene IDs:", gene_ids[1:10])
Gene IDs: ['ENSG00000223972.5', 'ENSG00000227232.5', 'ENSG00000278267.1', 'ENSG00000243485.5', 'ENSG00000237613.2', 'ENSG00000268020.3', 'ENSG00000240361.1', 'ENSG00000186092.4', 'ENSG00000238009.6']print("Gene Names:", gene_names[1:10])
print("Gene Names:", gene_names[1:10])
#Gene Names: ['DDX11L1', 'WASH7P', 'MIR6859-1', 'MIR1302-2HG', 'FAM138A', 'OR4G4P', 'OR4G11P', 'OR4F5', 'RP11-34P13.7']
print("Expression Data:", expression_data[1:10,1:5])


# Question2
import numpy as np
# Convert lists to numpy arrays
tissue_names = np.array(tissue_names)
gene_ids = np.array(gene_ids)
gene_names = np.array(gene_names)
expression_data = np.array(expression_data, dtype=float) 
# Print the results (for debugging purposes)
print("Tissue Names (numpy array):", tissue_names)
print("Gene IDs (numpy array):", gene_ids)
print("Gene Names (numpy array):", gene_names)
print("Expression Data (numpy array, first 5 rows):", expression_data[:5])



# Question 3
import numpy as np
# Number of genes to process
num_genes = 10
# Initialize a list to store mean expression values for the first 10 genes
mean_expression_values = []
# Iterate over the first 10 genes (rows)
for i in range(num_genes):
    # Initialize a variable to accumulate the sum of expression values for the current gene
    total_expression = 0.0
    # Iterate over all the tissues (columns) for the current gene
    for j in range(expression_data.shape[1]):
        total_expression += expression_data[i, j]
    
    # Calculate the mean expression value for the current gene
    mean_expression = total_expression / expression_data.shape[1]
    
    # Append the mean expression value to the list
    mean_expression_values.append(mean_expression)

# Print the mean expression values for the first 10 genes
for i, mean_expression in enumerate(mean_expression_values):
    print(f"Mean expression for gene {i+1}: {mean_expression:.2f}")
## Result:
#Mean expression for gene 1: 3.76
#Mean expression for gene 2: 0.00
#Mean expression for gene 3: 0.01
#Mean expression for gene 4: 0.00
#Mean expression for gene 5: 0.02
#Mean expression for gene 6: 0.04
#Mean expression for gene 7: 0.05
#Mean expression for gene 8: 0.03
#Mean expression for gene 9: 0.06
#Mean expression for gene 10: 8.11



## Question 4
import numpy as np
mean_expression_values = np.mean(expression_data[:10, :], axis=1)
# Print the mean expression values for the first 10 genes
for i, mean_expression in enumerate(mean_expression_values):
    print(f"Mean expression for gene {i+1}: {mean_expression:.2f}")
## Result:
# Mean expression for gene 1: 3.76
# Mean expression for gene 2: 0.00
# Mean expression for gene 3: 0.01
# Mean expression for gene 4: 0.00
# Mean expression for gene 5: 0.02
# Mean expression for gene 6: 0.04
# Mean expression for gene 7: 0.05
# Mean expression for gene 8: 0.03
# Mean expression for gene 9: 0.06
# Mean expression for gene 10: 8.11
#Q:Do they match the means you found with the nested for loop approach?
#Re: Yes!


## Question 5
import numpy as np
mean_expression = np.mean(expression_data)
print(f"Mean expression value of the entire dataset: {mean_expression:}")
### Mean expression value of the entire dataset: 16.56

# Calculate the median expression value of the entire dataset
median_expression = np.median(expression_data)

print(f"Median expression value of the entire dataset: {median_expression:.2f}")
### Median expression value of the entire dataset: 0.03

#Q:What can you infer from the difference between these statistics?
#Re: it suggests that there might be a significant number of very high gene expression across these tissues values skewing the mean upwards.


## Question 6
# Add a pseudo-count of 1 to avoid taking the log of zero
pseudo_count = 1
expression_data_with_pseudo = expression_data + pseudo_count

# Apply log2 transformation
log_expression_data = np.log2(expression_data_with_pseudo)

# Calculate the mean and median for the transformed data
mean_log_expression = np.mean(log_expression_data)
median_log_expression = np.median(log_expression_data)

# Print the results for the transformed data
print(f"Mean log2-transformed expression value of the entire dataset: {mean_log_expression:.2f}")
## Result: Mean log2-transformed expression value of the entire dataset: 1.12
print(f"Median log2-transformed expression value of the entire dataset: {median_log_expression:.2f}")
## Result: Median log2-transformed expression value of the entire dataset: 0.04



# Question 7
import numpy as np
#To do this, you will need to create a copy of the expression array since you will be altering the array to find this information. In order to make a copy of an array, you can use the function numpy.copy.
expression_data_copy = np.copy(log_expression_data)
# Sort the expression values along the tissue axis (axis=1)
sorted_expression_data = np.sort(expression_data_copy, axis=1)
# Calculate the difference between the highest and the second highest expression values
# For each gene (row), the last value is the highest and the second to last value is the second highest
diff_array = sorted_expression_data[:, -1] - sorted_expression_data[:, -2]
# array([9.25200e-02, 0.00000e+00, 2.35543e-02, ..., 1.08270e+03,2.60095e+00, 3.80328e+00])

#Print the results
print("Expression gaps (highest - second highest) for each gene:")
for i, diff in enumerate(diff_array):
    print(f"Gene {i+1}: {diff:.2f}")



# Question 8
specificity_threshold = 10
# Count the number of genes where the difference is greater than or equal to 10
high_specificity_count = np.sum(diff_array >= specificity_threshold)
print(f"Number of genes with high single-tissue specificity (difference >= {specificity_threshold}): {high_specificity_count}")
#Number of genes with high single-tissue specificity (difference >= 10): 33


#Advanced exercises
# Question 9
# Create a zero matrix of the same shape as expression_data
max_expression_matrix = np.zeros_like(log_expression_data)
# Find the indices of the highest expression values for each tissue (column)
# axis=0 gives us the index of the maximum value for each column (tissue) across all rows (genes)
max_indices = np.argmax(log_expression_data, axis=0)
# Create an index array for the rows, which will be the same length as the number of tissues
rows = np.arange(log_expression_data.shape[0])
# Set the highest expression values in the zero matrix to 1
max_expression_matrix[max_indices, np.arange(log_expression_data.shape[1])] = 1
high_tissue=log_expression_data[max_indices,]



# Question 10
# Define the threshold for high tissue specificity
specificity_threshold = 10
#Q8: diff_array = sorted_expression_data[:, -1] - sorted_expression_data[:, -2]
mask = diff_array > specificity_threshold
# Reshape the boolean array to match the dimensions of high_tissue
# Reshape to (num_genes, 1) to align with the second axis of high_tissue
mask_reshaped = mask.reshape(-1, 1)
# Multiply the high_tissue array by the reshaped boolean array
filtered_high_tissue = high_tissue * mask_reshaped
# Print the results
print("Filtered high_tissue array with specific genes retained:")
print(filtered_high_tissue)



# Question 11
specificity_threshold = 10
high_specificity_mask = diff_array > specificity_threshold
# Create a mask for the high_tissue array by broadcasting
mask_reshaped = high_specificity_mask.reshape(-1, 1)
filtered_high_tissue = high_tissue * mask_reshaped

# Find the indices of gene-tissue pairs where the value is 1
gene_indices, tissue_indices = np.where(filtered_high_tissue == 1)

print("Gene-Tissue pairs with high expression and specificity:")
for gene_idx, tissue_idx in zip(gene_indices, tissue_indices):
    gene_name = gene_names[gene_idx]
    tissue_name = tissue_names[tissue_idx]
    print(f"Gene: {gene_name}, Tissue: {tissue_name}")

# Optional: Save the results to a file
with open('high_specificity_genes_simplified.txt', 'w') as f:
    for gene_idx, tissue_idx in zip(gene_indices, tissue_indices):
        gene_name = gene_names[gene_idx]
        tissue_name = tissue_names[tissue_idx]
        f.write(f"Gene: {gene_name}, Tissue: {tissue_name}\n")

print("Results saved to 'high_specificity_genes_simplified.txt'")



