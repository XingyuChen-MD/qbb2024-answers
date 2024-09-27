#!/bin/bash

# Step 2: Count feature SNPs and determine enrichment

# Create the results file with a header
echo -e "MAF\tFeature\tEnrichment" > snp_counts.txt

# Loop through each possible MAF value
for MAF in 0.1 0.2 0.3 0.4 0.5
do 
    MAF_file=chr1_snps_${MAF}.bed # File name for the SNP MAF file

    # Find the SNP coverage of the whole chromosome
    bedtools coverage -a genome_chr1.bed -b ${MAF_file} > coverage_${MAF}.txt 

    # Sum SNPs from coverage (assuming the 4th column is SNP counts)
    num_SNP=$(awk '{s+=$4} END {print s}' coverage_${MAF}.txt)
    # Sum total bases from coverage (assuming the 6th column is the base count)
    num_bases=$(awk '{s+=$6} END {print s}' coverage_${MAF}.txt)

    # Calculate the background SNP density (SNPs per base in the whole genome)
    density=$(echo "scale=6; ${num_SNP}/${num_bases}" | bc -l)

    # Loop through each feature
    for feature in exons introns cCREs other
    do
        # Get the file name for the feature
        feature_file=${feature}_chr1.bed

        # Find the SNP coverage of the current feature
        bedtools coverage -a ${feature_file} -b ${MAF_file} > ${MAF}_${feature}.txt

        # Sum SNPs from feature coverage (4th column)
        feature_num_SNP=$(awk '{s+=$4} END {print s}' ${MAF}_${feature}.txt)
        # Sum total bases from feature coverage (6th column)
        feature_num_bases=$(awk '{s+=$6} END {print s}' ${MAF}_${feature}.txt)

        # Calculate the SNP density for the current feature
        feature_density=$(echo "scale=6; ${feature_num_SNP}/${feature_num_bases}" | bc -l)
        
        # Calculate the enrichment as the ratio of feature density to background density
        enrichment=$(echo "scale=6; ${feature_density}/${density}" | bc -l)

        # Save the result to the output file
        echo -e "${MAF}\t${feature}\t${enrichment}" >> snp_counts.txt
    done
done
    