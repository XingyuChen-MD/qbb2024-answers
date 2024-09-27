#!/bin/bash

# Step 1: Obtain the data

# Step 1.4: Sort and Merge Feature Files (exons, genes, cCREs)
for FEATURE in exons genes cCREs  # Loops through the three features
do
    # Sort the BED file by chromosome and coordinate
    bedtools sort -i ${FEATURE}.bed > ${FEATURE}_sort.bed
    # Merge overlapping regions to eliminate redundancy
    bedtools merge -i ${FEATURE}_sort.bed > ${FEATURE}_chr1.bed
done

# Step 1.5: Identify Introns by subtracting exon regions from gene regions
bedtools subtract -a genes_chr1.bed -b exons_chr1.bed > introns_chr1.bed

# Step 1.6: Identify "other" regions (neither exons, cCREs, nor introns)
bedtools subtract -a genome_chr1.bed -b exons_chr1.bed -b cCREs_chr1.bed -b introns_chr1.bed > other_chr1.bed
