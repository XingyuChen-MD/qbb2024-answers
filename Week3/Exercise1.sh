#!/usr/bin/env bash

# Task 1
## 1.1 
awk 'NR%4 == 2 {print length($0)}' A01_09.fastq
# This command isolates every second line in each group of four, corresponding to the sequencing read, and outputs its length.
# Sequence read lengths are found to be 76 nucleotides.

## 1.2 
wc -l A01_09.fastq | awk '{print $1 / 4}'
# Divides the total number of lines in the FASTQ file by 4 to estimate the total number of reads.
# The resulting count is 669,548 reads.

## 1.3 
echo "scale=3; 76 * 669548 / 12200000" | bc
# Computes the sequencing depth given a reference genome size of 12.2 million base pairs.
# The depth of coverage is approximately 4.171x.

## 1.4 
du -h *.fastq
# This gives a human-readable summary of file sizes for all FASTQ files in the directory.
# File sizes vary, with the largest file being 149MB and the smallest 110MB.

## 1.5 
fastqc *.fastq
# Runs FastQC for quality control, generating HTML reports to assess base quality, GC content, and more.
# Median base quality hovers around 36, implying an error rate of roughly 1 in 3981 bases.
