# gwaR

A collection of R functions to conduct downstream GWAS analysis for Arabidopsis thaliana.

# About

This package provides functions for fast analysis of GWAS tables. This package interacts with the 1001genomes API for variant tables by default, but can also deal with local SNPmatrix (needs to be fst format). Annotations are fetched from biomart upon package installation, and should be TAIR10. The overall idea is to be able to have a package that integrates tasks that are usually carried out by visiting different websites and trying to collect and curate information.
This package provides no new classes, but the table containing the SNP information (chromosome, position, pvalue) needs to conform to the expectations. There are wrappers that import csv files to match the desired format.

# Functions

This pacakge provides a generic plotting function for manhattan plots, if desired peaks can be annotated with closest, or overlapping genomic features. Overlapping, or nearest, genomic features can also be retrieved in table format. 

gwaR is aimed at making it easy to understand the effect of SNPs on the phenotype of interest. For this, gwaR provides functions that first retrieve which accessions contain a SNP, and then split the phenotype table by SNP presence. In addition, gwaR can make use of 1001genomes.org to retrieve expression data, and integrate this information with SNP information, for example to investigate if SNPs in promoters or UTRs have an impact on expression.

Finally, it provides a function to compute and plot linkage between a SNP from the gwas table, and all known SNPs around it. This can again be done using either 1001genomes information (however, this often fails), or a local SNPmatrix. SNPIMPACTS (low, mid, high) are in any case retrieved from 1001genomes.org

As of recently, this package is also able to query thalemine, to retreive curated annotations and a list of publications referencing the gene of interest.

# Other
This package was developed and used to analyze data in [aradeepopsis: From images to phenotypic traits using deep transfer learning] (https://www.biorxiv.org/content/10.1101/2020.04.01.018192v2) and works well with [aradeepopsis pipeline](https://github.com/Gregor-Mendel-Institute/aradeepopsis)
