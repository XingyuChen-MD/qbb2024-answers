#!/usr/bin/env python3

# Step 1.6
import numpy as np
genome_size = 1_000_000  # 1Mbp
read_length = 100  # 100bp
coverage_30x = 30  # 30x coverage

num_reads_30x = int((coverage_30x * genome_size) / read_length)

genome_coverage_30x = np.zeros(genome_size, dtype=int)
for _ in range(num_reads_30x):
    start_pos = np.random.randint(0, genome_size - read_length + 1)
    genome_coverage_30x[start_pos:start_pos + read_length] += 1

output_file_path_30x = 'genome_coverage_30x.txt'
np.savetxt(output_file_path_30x, genome_coverage_30x, fmt='%d')
