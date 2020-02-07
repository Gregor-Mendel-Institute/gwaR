% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GWAS_functions.R
\name{polymorph_impact}
\alias{polymorph_impact}
\title{Retrieve impacts of SNPs within a certain range around a GWAS peak. Uses snpmatrix for genotype lookup}
\usage{
polymorph_impact(
  gwas_table,
  rank,
  nuc_range,
  snpmatrix_path,
  use_phenotype_table = NULL,
  acc_col
)
}
\arguments{
\item{gwas_table}{Object returned from \code{\link{read_gwas}} function}

\item{rank}{The (-log10(p)) rank of the SNP of interest}

\item{nuc_range}{Window around the SNP of interest to retrieve information for.}

\item{snpmatrix_path}{The path to the snpmatrix to use for identifying accessions that carry the relevant SNP.}

\item{use_phenotype_table}{If supplied: Genotypes listed here will be used for linkage analysis. Otherwise, all accessions that carry this SNP will be included.}

\item{acc_col}{the column that contains accession identifiers.}
}
\description{
Retrieve impacts of SNPs within a certain range around a GWAS peak. Uses snpmatrix for genotype lookup
}
\details{
Based on a Phenotype table, the name of the phenotype and a GWAS table and a rank, produces a boxplot of that phenotype, grouped by presence of that SNP.
}
\seealso{
\code{\link{get_phenotype}}

\code{\link{read_gwas}}

\code{\link{intersect_phenotype_snp}}
}