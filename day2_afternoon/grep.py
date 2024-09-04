#!/usr/bin/env python3

import sys

my_file  = open(sys.argv[2])
gene = sys.argv[1]

for i in my_file:
    i = i.rstrip("\n")
    if gene in i:
        print(i)

my_file.close()