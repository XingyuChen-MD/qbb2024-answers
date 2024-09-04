#!/usr/bin/env python3

import sys

my_file  = open(sys.argv[1])

for line in my_file:
    if "##" in line:  ## skip the description lines 
        continue
    i = line.split("\t")
    chr = i[0]
    start = i[3]
    stop = i[4]
    info = i[8] # the last colume has id, name, type... info
    info_list = info.split(";") # split all info make it list
    
    x = 0
    gene_name = ""
    for a in info_list:
        if "gene_name" in a:
            a = a.lstrip('gene_name "').rstrip('"')
            gene_name = a
            print(chr + "\t" + start + "\t" + stop + "\t"+ gene_name)
            x = x + 1

my_file.close()