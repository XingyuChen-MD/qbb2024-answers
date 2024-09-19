#!/usr/bin/env python3

# Step 1.5

import numpy as np

genome_size = 1_000_000  # 1Mbp
read_length = 100  # 100bp

coverage_10x = 10
num_reads_10x = int((coverage_10x * genome_size) / read_length)

genome_coverage_10x = np.zeros(genome_size, dtype=int)

for _ in range(num_reads_10x):
    start_pos = np.random.randint(0, genome_size - read_length + 1)
    genome_coverage_10x[start_pos:start_pos + read_length] += 1

output_file_path_10x = 'genome_coverage_10x.txt'
np.savetxt(output_file_path_10x, genome_coverage_10x, fmt='%d')
