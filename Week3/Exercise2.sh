#!/usr/bin/env bash

# Task 2
## 2.1
wget https://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
gunzip sacCer3.fa.gz
# Download and decompress the reference genome for Saccharomyces cerevisiae.
bwa index sacCer3.fa
# Creates a BWA index for efficient alignment of sequencing reads to the reference genome.

### 2.1
grep "^>" sacCer3.fa
# Extracts chromosome names from the reference file by searching for lines beginning with '>'.
grep -c "^>" sacCer3.fa
# Counts the number of chromosomes in the file, with or without the mitochondrial genome, resulting in either 16 or 17 chromosomes.

## 2.2
for sample in A01_*.fastq
do 
    sample_name=`basename ${sample} .fastq`
    # Loop through each sample file and extract the base name.
    bwa mem -t 3 -R "@RG\tID:${sample_name}\tSM:${sample_name}" sacCer3.fa ${sample}.fastq > ${sample}.sam
    # Align reads to the reference genome using BWA and output in SAM format.
    samtools sort -@ 4 -O bam -o ${sample_name}.bam ${sample_name}.sam
    # Convert the SAM file to a sorted BAM file.
    samtools index ${sample_name}.bam
    # Index the BAM file for fast access to specific regions.
done

## 2.3
### 2.2
grep "^HWI-ST387" A01_09.sam | wc -l
# Counts the number of read alignments in the SAM file, resulting in 669,548 alignments.

### 2.3
grep -w "chrIII" A01_09.sam | wc -l
# Filters alignments on chromosome III, yielding a total of 18,196 reads mapped to this chromosome.

### 2.4
# Alignments and BAM file preparation are covered in previous steps.

### 2.5
# While average depth can be approximated, read coverage varies significantly across the genome, with some regions having no coverage and others exceeding the expected depth.

### 2.6
# Three single nucleotide polymorphisms (SNPs) are identified. Two of them are well-supported by four reads, while one is supported by only two reads, raising concerns about sequencing accuracy.
# The SNP is located at position chrIII:825834, outside any gene but between SCC2 and SAS4.
