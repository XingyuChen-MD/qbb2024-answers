#!/bin/bash

echo "Hello, World!"

#unit stores command line argument in $0, $1
echo "0th: " $0
echo "1st: " $1

cut -f 1 $q >column.txt
sort column.txt | uniq -c