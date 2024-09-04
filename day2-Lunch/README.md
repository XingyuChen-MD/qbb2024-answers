# (1) Instructions
### Complete the following exercises using a combination of cut, sort, uniq, and grep. Document each command you use, along with any short answers, in a single README.md file stored in ~/qbb2024-answers/day2-lunch. Use Markdown headings, lists, and code blocks to organize your answers into sections and format content for proper display e.g.

# 1. Extracting Specific Columns
## Extract specific columns from hg38-gene-metadata-feature.tsv and hg38-gene-metadata-go.tsv.
### for hg38-gene-metadata-feature.tsv: Extract columns 1 and 3. for `hg38-gene-metadata-go.tsv`: Extract columns 2 and 4. 
```
cut -f1,3 hg38-gene-metadata-feature.tsv > feature_columns.txt
cut -f2,4 hg38-gene-metadata-go.tsv > go_columns.txt
head feature_columns.txt
```
```
(qb24) cmdb@QuantBio-26 Day2Morning % head feature_columns.txt
ensembl_gene_id	chromosome_name
ENSG00000228037	1
ENSG00000142611	1
ENSG00000284616	1
ENSG00000157911	1
ENSG00000260972	1
ENSG00000224340	1
ENSG00000226374	1
ENSG00000229280	1
ENSG00000142655	1 
```

# 2. Sorting and Removing Duplicates
## Sort the contents of length.txt and remove duplicate lines.
```
sort length.txt | uniq > unique_sorted_length.txt
head unique_sorted_length.txt
```
```
(qb24) cmdb@QuantBio-26 Day2Morning % head unique_sorted_length.txt
   61633 hg38-gene-metadata-feature.tsv
   69601 hg38-gene-metadata-homologs.tsv
  293343 hg38-gene-metadata-go.tsv
  424577 total
  ```
## sort organizes the lines alphabetically, and uniq removes any repeated lines that are adjacent.

# 3. Filtering Lines with Specific Patterns
## Task: Filter lines containing the pattern "gene" from gencode.v46.basic.annotation.gtf.
```
grep "gene" gencode.v46.basic.annotation.gtf > genes_only.txt
head genes_only.txt
```

# 4. Combining Multiple Commands
#Task: Extract the first column from hg38-gene-metadata-homologs.tsv, sort the results, and then remove duplicates.
```
cut -f1 hg38-gene-metadata-homologs.tsv | sort | uniq > unique_homologs.txt
```
#cut extracts the first column, sort organizes the lines, and uniq removes duplicates.

# 5. Verifying File Contents
#Task: Check the contents of sorted_length.txt to ensure it's properly sorted.

`cat sorted_length.txt`

# 6 Answer 1
```
(qb24) cmdb@QuantBio-26 Day2Morning % wc -l hg38-gene-metadata-feature.tsv
   61633 hg38-gene-metadata-feature.tsv
```


# (2) Exercises
## Tally the number of each gene_biotype in hg38-gene-metadata-feature.tsv. How many protein_coding genes are there? Pick one biotype you would want to learn more about and explain why.
```
cut -f7 hg38-gene-metadata-feature.tsv | uniq -c > biotype_counts.txt
```

## Which ensembl_gene_id in hg38-gene-metadata-go.tsv has the most go_ids? Create a new file that only contains rows corresponding to that gene_id, sorting the rows according to the name_1006 column. Describe what you think this gene does based on the GO terms.
```
(qb24) cmdb@QuantBio-26 Day2Morning % cut -f1 hg38-gene-metadata-go.tsv | uniq -c | sort | > ensembl_gene_id_counts.txt
(qb24) cmdb@QuantBio-26 Day2Morning % tail ensembl_gene_id_counts.txt
 183 ENSG00000097007
 184 ENSG00000142208
 187 ENSG00000197122
 201 ENSG00000141510
 207 ENSG00000114251
 221 ENSG00000148400
 224 ENSG00000232810
 227 ENSG00000105329
 234 ENSG00000164690
 273 ENSG00000168036
```
### ENSG00000168036 has the most go_ids, 273 go_ids

```
(qb24) cmdb@QuantBio-26 Day2Morning % grep -w "ENSG00000168036" hg38-gene-metadata-go.tsv | sort -u -k3 > ENSG00000168036_only.txt 
(qb24) cmdb@QuantBio-26 Day2Morning % head ENSG00000168036_only.txt
ENSG00000168036	GO:0005515	protein binding
ENSG00000168036	GO:0007155	cell adhesion
ENSG00000168036	GO:0045296	cadherin binding
ENSG00000168036	GO:0005737	cytoplasm
ENSG00000168036	GO:0005634	nucleus
ENSG00000168036	GO:0016020	membrane
ENSG00000168036	GO:0042995	cell projection
ENSG00000168036	GO:0005886	plasma membrane
ENSG00000168036	GO:0005856	cytoskeleton
ENSG00000168036	GO:0016055	Wnt signaling pathway
```
#### This is gene function may related to Cellular Processes and Functions, Cell Morphology and Adhesion, and Cellular Components and Structures, and et al


## Immunoglobin (Ig) genes are present in over 200 copies throughout the human genome. How many IG genes (not pseudogenes) are present on each chromosome? You can use a dot (.) in a regular expression pattern to match any single character. How does this compare with the distribution of IG pseudogenes?

`grep -w -e "IG...gene" -e "IG....gene" gene.gtf | cut -f1 | sort | uniq -c`
```
sort | uniq -c
  91 chr14
  16 chr15
   6 chr16
  52 chr2
   1 chr21
  48 chr22
```
```
`grep -w -e "IG.pseudogene" -e "IG...pseudogene" gene.gtf | cut -f1 | sort | uniq -c`
   1 chr1
   1 chr10
  84 chr14
   6 chr15
   8 chr16
   1 chr18
  45 chr2
  48 chr22
   1 chr8
   5 chr9
````


## Why is grep pseudogene gene.gtf not an effective way to identify lines where the gene_type key-value pair is a pseudogene (hint: look for overlaps_pseudogene)? What would be a better pattern? Describe it in words if you are having trouble with the regular expression.

Re: The term pseudogene might appear in multiple contexts in the GTF file. It could be part of other fields or annotations not related to the gene_type specifically, leading to incorrect or irrelevant results.
For example: chr1	HAVANA	gene	12010	13670	.	+	.	gene_id "ENSG00000223972.6"; gene_type "transcribed_unprocessed_***pseudogene***"; gene_name "DDX11L1"; level 2; hgnc_id "HGNC:37102"; havana_gene "OTTHUMG00000000961.2";
To find lines where gene_type is exactly pseudogene, use a regular expression that looks for gene_type "pseudogene" within the attributes field.

`grep -E 'gene_type "pseudogene"' gene.gtf`
`grep -w "pseudogene*" gene.gtf`

## Convert the annotation from .gtf format to .bed format. Specifically, print out just the chromosome, start, stop, and gene_name. As cut splits lines into fields based on the tab character, first use sed to create a new file where spaces are replaced with tabs.
```
sed "s/ /\t/g" gene.gtf > gene-tabs.gtf
less -S gene-tabs.gtf 
# we found that the chromosome, start, stop, and gene_name are located in  column 1,4,5,14
cut -f 1,4,5,14 > new_gene_tabs.bed
```
