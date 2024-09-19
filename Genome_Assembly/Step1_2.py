#!/usr/bin/env python3

# Step 1.2
import numpy as np

# Parameters
genome_size = 1_000_000  # 1Mbp
read_length = 100  # 100bp
coverage = 3  # 3x coverage

# Calculate number of reads needed
num_reads = int((coverage * genome_size) / read_length)

# Initialize an array to keep track of coverage at each position in the genome
genome_coverage = np.zeros(genome_size, dtype=int)

# Generate random start positions for each read
start_positions = np.random.randint(0, genome_size - read_length + 1, size=num_reads)

# Efficiently update the genome coverage using np.add.at
for start_pos in start_positions:
    np.add.at(genome_coverage, np.arange(start_pos, start_pos + read_length), 1)

# Save the genome coverage array to a text file
output_file_path = 'genome_coverage.txt'
np.savetxt(output_file_path, genome_coverage, fmt='%d')
