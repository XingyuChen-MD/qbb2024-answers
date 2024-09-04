#!/usr/bin/env python3

import sys

print(sys.argv)

num_arg = len(sys.argv)

print(num_arg)
print("Number of arguments: " + str(num_arg))

i = 0
for my_arg in sys.argv:
    print(str(my_arg) + " the argument is " + sys.argv[i])
    i = i + 1

#(qb24) cmdb@QuantBio-26 qbb2024-answers % python3 03_arguments.py apple banana coconut
