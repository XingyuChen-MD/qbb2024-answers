#!/usr/bin/env python3

#Import library
import sys
import numpy as np 
import pandas as pd 

#Create empty list for AF and DP
AF = []
DP = []

#Q3
##3.1
#Open the VCF file 
for line in open(sys.argv[1]):
    #Remove the lines start with #
    if line.startswith("#"): 
        continue 
    else: 
        #Divide the columns with tab
        fields = line.rstrip('\n').split('\t')
        #Select info column
        info_f = fields[7]
        #Divide the info column by semi-colon
        columns = info_f.rstrip('\n').split(';')
        #Select column 4 for allele frequency 
        AF.append(columns[3])
        #Create for loop for all the variant columns
        for variant in range(9, 19): 
            column = fields[variant]
            #Divide variant data columns by colon 
            data = column.rstrip('\n').split(':')
            #Select column 3 for the read depth
            DP.append(data[2])
        

#Write a file containing header and allele frequency data (only select the number)
with open(f"AF.txt", "w") as file:
    file.write("Allele_Frequency\n")
    for freq in AF:
        file.write(f"{freq[3:]}\n")

#Write a file containing header and read depth data
with open(f"DP.txt", "w") as file:
    file.write("Read_Depth\n")
    for depth in DP:
        file.write(f"{depth}\n")

###Q3.1
#We expected to see such distribution since the allele frequencies should favor heterozygous as to having only one or the other allele
#The distribution is a binomial distribution 

###Q3.2
#The distribution is expected as we have read depth averaging on the lower ends but still has regions that have a substantial amount of reads
#The distribution is a poisson distribution

        