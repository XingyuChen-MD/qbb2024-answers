#!/usr/bin/env python3
import sys

my_file= opne ( sys.argv[1])
 for  line in my_file:
    line =line.restirp("\n")
    print(line)
my_file.close()


i=0
for line in my_file:
    if  i=10:
        break
    line = line.rstrip("\n")
    print( line )
    i=i+1
my_file.close()

import sys
if len(sys.argv) != 2:
    print("Usage: python script.py filename")
    sys.exit(1)

filename = sys.argv[1]

try:
    with open(filename, 'r') as my_file:
        for line in my_file:
            print(line, end='')  
except FileNotFoundError:
    print(f"Error: The file '{filename}' does not exist.")
except IOError as e:
    print(f"Error: An IOError occurred. {e}")



