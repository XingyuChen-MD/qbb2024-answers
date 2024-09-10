#!/usr/bin/env python3

import sys
import numpy as np
 # get file name

#Q1

filename = sys.argv[1]  # Retrieve the filename for the gene-tissue pairs from command-line arguments
fs = open(filename, mode='r')  # Open the file in read mode

# Create a dictionary to hold gene-tissue pairs
relevant_samples = {}

# Iterate through each line in the file
for line in fs:
    fields = line.rstrip("\n").split("\t")  # Split the line into fields using tab as the delimiter
    key = (fields[0], fields[2])  # Create a key from geneID (fields[0]) and tissue type (fields[2])
    relevant_samples[key] = []  # Initialize an empty list for each gene-tissue pair
fs.close()  # Close the file after reading all lines

print("Relevant Samples:", relevant_samples)  # Print the dictionary to verify the contents
#Explanation:Opening the File: The open(filename, mode='r') function opens the file for reading.
#Dictionary Initialization: relevant_samples is a dictionary where keys are tuples of (geneID, tissue type) and values are empty lists, initialized for future data.
#Processing Each Line: Each line is split into fields, and a tuple of (geneID, tissue type) is used as the key in the dictionary.
#How to run:
##chmod +X Day4_lunch.py   
##python Day4_lunch.py GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

#Q2
filename = sys.argv[2]  # Retrieve the filename for the sample attributes from command-line arguments
fs = open(filename, mode='r')  # Open the file in read mode

fs.readline()  # Skip the header line

# Create a dictionary to hold tissue samples
tissue_samples = {}

# Iterate through each line in the file
for line in fs:
    fields = line.rstrip("\n").split("\t")  # Split the line into fields using tab as the delimiter
    key = fields[6]  # Tissue type is in column 6
    value = fields[0]  # Sample ID is in column 0
    tissue_samples.setdefault(key, [])  # Ensure a list exists for the tissue type
    tissue_samples[key].append(value)  # Append the sample ID to the list for the tissue type
fs.close()  # Close the file after reading all lines

print("Tissue Samples:", tissue_samples)  # Print the dictionary to verify the contents


# Q3 You now need to get the list of sampleIDs that are present in the gene expression file GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.
# RNA-seq data
filename=sys.argv[3] #less -s GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
# sampleId is in coloum1
fs = open(filename, mode='r') # Open the file in read mode
#Skip line
fs.readline()# Skip the second header line
fs.readline()
# creat dict to hold samoles for gene-tissue pairs
header = fs.readline().rstrip("\n").split("\t")
header = header[2:]# Get column indices for sample IDs (excluding first two columns)

print(header)


#python Day4_lunch.py gene_tissue.tsv GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct

sample_id_to_column = {sample_id: idx for idx, sample_id in enumerate(header)}

# Q4: Map each tissue to the relevant columns in the expression data
tissue_columns = {}  # Initialize a dictionary to hold tissue-to-column mappings

for tissue, samples in tissue_samples.items():
    tissue_columns.setdefault(tissue, [])  # Ensure a list exists for each tissue type
    for sample in samples:
        if sample in sample_id_to_column:  # Check if sample ID is in the expression data
            position = sample_id_to_column[sample]  # Get the column index of the sample
            tissue_columns[tissue].append(position)  # Append the column index to the list for the tissue type

print("Tissue Columns:", tissue_columns)  # Print the dictionary to verify the contents


# Usage: python Day4_lunch.py gene_tissue.tsv GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct

# Q5: Analyze the number of samples for each tissue type
tissue_sample_counts = {tissue: len(indices) for tissue, indices in tissue_columns.items()}

# Identify tissues with the most and fewest samples
max_samples_tissue = max(tissue_sample_counts, key=tissue_sample_counts.get)
min_samples_tissue = min(tissue_sample_counts, key=tissue_sample_counts.get)

print("Number of samples per tissue type:", tissue_sample_counts)
print(f"Tissue with the most samples: {max_samples_tissue} ({tissue_sample_counts[max_samples_tissue]} samples)")
print(f"Tissue with the fewest samples: {min_samples_tissue} ({tissue_sample_counts[min_samples_tissue]} samples)")
#### Results: 
#Number of samples per tissue type: {'Whole Blood': 755, 'Brain - Frontal Cortex (BA9)': 209, 'Adipose - Subcutaneous': 663, 'Muscle - Skeletal': 803, 'Artery - Tibial': 663, 'Artery - Coronary': 240, 'Heart - Atrial Appendage': 429, 'Adipose - Visceral (Omentum)': 541, 'Ovary': 180, 'Uterus': 142, 'Vagina': 156, 'Breast - Mammary Tissue': 459, 'Skin - Not Sun Exposed (Suprapubic)': 604, 'Minor Salivary Gland': 162, 'Brain - Cortex': 255, 'Adrenal Gland': 258, 'Thyroid': 653, 'Lung': 578, 'Spleen': 241, 'Pancreas': 328, 'Esophagus - Muscularis': 515, 'Esophagus - Mucosa': 555, 'Esophagus - Gastroesophageal Junction': 375, 'Stomach': 359, 'Colon - Sigmoid': 373, 'Small Intestine - Terminal Ileum': 187, 'Colon - Transverse': 406, 'Prostate': 245, 'Testis': 361, 'Skin - Sun Exposed (Lower leg)': 701, 'Nerve - Tibial': 619, 'Heart - Left Ventricle': 432, 'Pituitary': 283, 'Brain - Cerebellum': 241, 'Cells - Cultured fibroblasts': 504, 'Artery - Aorta': 432, 'Cells - EBV-transformed lymphocytes': 174, 'Brain - Cerebellar Hemisphere': 215, 'Brain - Caudate (basal ganglia)': 246, 'Brain - Nucleus accumbens (basal ganglia)': 246, 'Brain - Putamen (basal ganglia)': 205, 'Brain - Hypothalamus': 202, 'Brain - Spinal cord (cervical c-1)': 159, 'Liver': 226, 'Brain - Hippocampus': 197, 'Brain - Anterior cingulate cortex (BA24)': 176, 'Brain - Substantia nigra': 139, 'Kidney - Cortex': 85, 'Brain - Amygdala': 152, 'Cervix - Ectocervix': 9, 'Fallopian Tube': 9, 'Cervix - Endocervix': 10, 'Bladder': 21, 'Kidney - Medulla': 4, 'Cells - Leukemia cell line (CML)': 0}
#Tissue with the most samples: Muscle - Skeletal (803 samples)
#Tissue with the fewest samples: Cells - Leukemia cell line (CML) (0 samples)

# Usage: python Day4_lunch.py gene_tissue.tsv GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct


# Q5: 6. Go through the code from steps 1 through 5 and comment each line or block of code (if multiple lines are working together and it doesn’t make sense to comment the individually) and explain what the purpose of each part of the code is doing.
#Your comments should demonstrate that you understand not only the purpose of a piece of code, buy also how it is accomplishing that purpose
###Re： Already did this.

#7. (formerly 8) Finally, you can visualize how variable each gene’s expression is. Download the data, load it into R, and create a violin plot of expression levels broken down by gene (ggplot2’s geom_violin()).
#For categories, create a combination of tissue names and gene IDs (dplyr::mutate(Tissue_Gene=paste0(Tissue, " ", GeneID)))
#You will need to log-transform your data with a psuedo-count of one (you can use dplyr::mutate for this step as well)
#Switch the axes for the violin plot so the categories are on the y-axis (coord_flip())
#Make sure to label your axes

#Given the tissue specificity and high expression level of these genes, are you surprised by the results?
#What tissue-specific differences do you see in expression variability? Speculate on why certain tissues show low variability while others show much higher expression variability.

