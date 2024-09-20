

# Step 1.1 In your README.md for this assignment, answer the following question (show your work): How many 100bp reads are needed to sequence a 1Mbp genome to 3x coverage?

Coverage (C) = (Number of reads × Length of each read) / Genome size
Coverage (C) = 3x
Genome size = 1,000,000 bp (1Mbp)
Length of each read = 100 bp
So, 30,000 reads of 100 bp are needed to sequence a 1Mbp genome to 3x coverage.


# Step 1.4 Using your results from Step 1.3, answer the following questions in your README.md:

# In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
# Rcode：
 zero_coverage_count <- sum(genome_coverage == 0)
> total_bases <- length(genome_coverage)
> percentage_unsequenced <- (zero_coverage_count / total_bases) * 100
> print(percentage_unsequenced)
[1] 5.0997
Proportion of Genome with 0x Coverage: 5.0997

Around 4.98% of the genome positions are expected to have 0x coverage under a Poisson distribution with λ = 3.
Comparison with Poisson Expectations: The observed proportion of genome positions with 0x coverage was 5.0997%, which closely matches the Poisson expectation of approximately 4.98%.

# How well does this match Poisson expectations? How well does the normal distribution fit the data?
The normal distribution provides a reasonable fit around the mean but does not accurately capture the probabilities for very low coverage values (e.g., 0 or 1). The Poisson distribution, on the other hand, aligns more closely with the observed histogram for lower coverage values. Overall, the Poisson distribution is a better fit for this data, particularly at low coverages.

# Step 1.5 In your README.md, answer the following questions:
## In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
> zero_coverage_count_10x <- sum(genome_coverage_10x == 0)
> total_bases_10x <- length(genome_coverage_10x)
> percentage_unsequenced_10x <- (zero_coverage_count_10x / total_bases_10x) * 100
> print(percentage_unsequenced_10x)
[1] 0.0045
Proportion of Genome with 0x Coverage: 0.0045%
## How well does this match Poisson expectations? How well does the normal distribution fit the data?
Comparison with Poisson Expectations: The observed percentage of 0x coverage (0.0045%) aligns perfectly with the Poisson expectation of approximately 0.0045%, suggesting that the simulation accurately reflects a Poisson distribution for 10x coverage.
Fit of the Normal Distribution: The normal distribution provides a reasonable fit around the mean and captures the overall shape of the histogram fairly well. However, the Poisson distribution still more accurately models the discrete nature of coverage counts, particularly at lower values.

# Step 1.6
# In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
zero_coverage_count_30x <- sum(genome_coverage_30x == 0)
> total_bases_30x <- length(genome_coverage_30x)
> percentage_unsequenced_30x <- (zero_coverage_count_30x / total_bases_30x) * 100
> print(percentage_unsequenced_30x)
[1] 2e-04
Proportion of Genome with 0x Coverage: Approximately 0.0002%.

# How well does this match Poisson expectations? How well does the normal distribution fit the data?
Comparison with Poisson Expectations: The observed percentage of 0x coverage is extremely close to the theoretical Poisson expectation of nearly 0%, indicating a strong match. This is consistent with the expected result that very high coverage (30x) leaves almost no part of the genome unsequenced.
Fit of the Normal Distribution: The normal distribution provides a very good approximation of the coverage data at 30x, closely matching the histogram's symmetric shape around the mean. Although the Poisson distribution remains more accurate for discrete counts, the normal distribution’s fit improves significantly as coverage increases.

# Step 2.2
conda create -n graphviz -c conda-forge graphviz
conda activate graphviz

# Step 2.3
Graphviz’s command line tool is called dot. Read more about how to use dot and the file format it’s expecting here.

NOTE: We’re going to want to produce a directed graph.
Based on what you’ve read, modify your code from Step 2.1 to output the edges in a format dot can use.

# Step 2.4
(qb24) cmdb@QuantBio-26 Desktop % conda activate graphviz
(graphviz) cmdb@QuantBio-26 Desktop % dot -Tpng debruijn_graph.dot -o ex2_digraph.png

# Step 2.5
One possible genome sequence: ATTCATTCATTGA

# Step 2.6: Requirements for Accurate Genome Reconstruction
To correctly piece together the genome sequence:

1.  Make sure that every part of the genome is included in the reads.
2. Pick a suitable length for k-mers. Longer k-mers help distinguish repeated sections in the genome but require more data to work correctly.
3. Handling Repeats: Use smart algorithms to deal with repeated sequences in the genome, such as finding paths in the de Bruijn graph.
4. High Sequencing Depth: Collect enough data (reads) to reduce mistakes and make sure overlaps between sequences are clear.
