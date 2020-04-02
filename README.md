# gwaR
A collection of R functions to conduct downstream GWAS analysis for Arabidopsis thaliana.
This package provides functions for fast analysis of GWAS tables. This package interacts with the 1001genomes API for variant tables by default, but can also deal with local SNPmatrix (needs to be fst format).
It provides a generic plotting function for manhattan plots, if desired peaks can be annotated with closest, or overlapping genomic features. Annotations are fetched from biomart upon package installation, should be TAIR10.
Overlapping or nearest features can be retrieved for the desired number of gwas hits. The pvalue table is sorted by -log10(pv), and the top rank is the one with the highest value (== lowest p-value).
If provided with a phenotype table, gwaR provides functions to split the phenotype by the presence of a SNP in the gwas table. It also provides a function to do this for the expression level of any gene, making again use of the 1001genomes API.
Finally, it provides a function to compute and plot linkage between a SNP from the gwas table, and all known SNPs around it. This can again be done using either 1001genomes information, or a local SNPmatrix. SNPIMPACTS are in any case retrieved from 1001genomes.org.

Niklas Schandry

