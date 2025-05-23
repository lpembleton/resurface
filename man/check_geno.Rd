% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{check_geno}
\alias{check_geno}
\title{Check Genotype Validity}
\usage{
check_geno(geno)
}
\arguments{
\item{geno}{Matrix of allele frequency genotype data.
Samples across rows, Loci down columns.}
}
\value{
Genotype matrix with any homozygous loci prefilled.
}
\description{
Function to check the validity of input geno and also fill any homozygous loci.
Note: as there is no variance for homozygous loci to correlate with, imputation
isnt possible and therefore the whole loci is recoded as as homozygous
}
