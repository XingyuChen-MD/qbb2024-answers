
#!/usr/bin/env python3

# Step 2.1

import numpy as np

reads = ['ATTCA', 'ATTGA', 'CATTG', 'CTTAT', 'GATTG', 'TATTT', 'TCATT', 'TCTTA', 'TGATT', 'TTATT', 'TTCAT', 'TTCTT', 'TTGAT']

edges = set()
k = 3

for read in reads:
    for i in range(len(read) - k):
        kmer1 = read[i:i+k]
        kmer2 = read[i+1:i+1+k]
        edge = f"{kmer1} -> {kmer2}"
        edges.add(edge)


with open('debruijn_edges.txt', 'w') as file:
    for edge in edges:
        file.write(edge + '\n')

print("Edges written to debruijn_edges.txt")

def write_dot_file(edges, output_file):
    with open(output_file, 'w') as f:
        f.write("digraph G {\n")
        
        for edge in edges:
            kmer1, kmer2 = edge.split(' -> ')
            f.write(f'    "{kmer1}" -> "{kmer2}";\n')
        
        f.write("}\n")

write_dot_file(edges, 'debruijn_graph.dot')

print("Graph written to debruijn_graph.dot")