# resurface

An R package for the imputation of genotypic allele frequencies using the LD-kNNi algorithm.

<img src="man/figures/resurface_high_200ms.gif" height="200"/>

## Overview

Many software packages exist to impute biallelic (AA, AB, BB) genotype calls.
However, allele frequencies from pooled population genotyping or autopolyploid genotypes are often overlooked.
resurface aims to fill this gap by providing a method to impute missing data in allele frequency datasets without requiring ordered genotypes (genetic or physical loci positions).

resurface is inspired by the original implementation of LD-kNNi for genotypic class imputation in the software [`LinkImpute`](http://www.cultivatingdiversity.org/software.html) by Daniel Money.
This methodology was extended to work on continuous allele frequencies rather than categorical genotype calls and implemented in R.
With sufficient loci in linkage disequilibrium (LD) in the dataset, resurface can achieve very high imputation accuracies.
Functions are included to empirically test the expected accuracy on your datasets.

## Installation

Installing is simple, you will just need to make sure you have the [remotes](https://github.com/r-lib/remotes) package from [cran](https://cran.r-project.org/web/packages/remotes/index.html) installed first.

``` r
install.packages("remotes") # if not already installed
remotes::install_github("lpembleton/resurface")
```

## Basic Usage

resurface accepts input in two formats:

-   A numeric matrix of allele frequencies (between 0 and 1) with sample IDs as row names, or
-   A 2-column tibble with sample IDs in column 1 and a matrix-column of allele frequencies in column 2. Samples are across rows, and loci are down columns.

> [!N
> OTE] If you are unfamiliar with matrix-columns in tibbles, check out [Eric Scott's](https://github.com/Aariq) [blog post](https://ericrscott.com/posts/2020-12-11-matrix-columns/) 👈.
> Matrix-columns are an efficient way to include large numeric matrices within tibble-like structures, allowing you to perform all the dplyr operations on other columns that might contain metadata such as population or location.

Missing data should be coded as `NA`.

There are only two functions that you need care about, `impute()` and `imputation_accuracy()` and they do what the name implies.

**`impute()`** will impute missing frequency data with the following parameters:

`k_neighbours` - Number of k-neighbours (samples) to use in the KNN approach (default: 10).

`l_loci` - Number of loci, ordered based on correlation to the SNP to be impute, to use when calculating nearest neighbours (default: 10).

`cpus` - Number of CPUs to use for parallel processing (default: 1).

`method` - Imputation method to use.
Currently the options are "ld_knni" or "ld_knni_fast" (default: "ld_knni_fast").

**`imputation_accuracy()`** simulates a user defined level of missing data and then imputes it in order to predict the expected level of imputation accuracy.
This is good for determining whether resurface is suitable for your data and also to optimise the `k_neighbours` and `l_loci` parameters.
The statistic that is currently return is Raw Bias (RB) which is calculate as the difference between the imputed allele frequency and the true allele frequency, averaged across all loci.
This function take two additional parameters:

`prop_na` The proportion of data to mask as missing in the dataset for accuracy testing (default: 0.2).

`repl` Number of replicate runs for imputation accuracy calculation (default: 1).

> [!W
> ARNING] Calculating the correlation matrix between all pairwise loci is time-consuming.
> Faster correlation methods exist but require no missing data.
> If the `ld_knni_fast` method is used missing data will be temporarily replaced with the mean loci allele frequency for correlation calculations.
> This can reduce imputation accuracy, especially in datasets with high missing data, while the impact is smaller in datasets with lower levels of missing data.

If you are unsure about the `k_neighbours` and `l_loci` values to use, test different combinations with the `imputation_accuracy()` function.

## Empircal Demonstration

Wondering how effective resurface is.
Below is an example (default parameters, completely un-optimised) utilising a publicly available dataset, so you can replicated it on your own setup.

### GBS allele frequency data from 1012 Alfalfa accessions

This dataset comes from [Pegard et al. 2023](https://doi.org/10.3389/fpls.2023.1196134) and is available to download from <https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/L0FLJD>

Simulated level of missing data: 40%

Imputation accuracy (1 - RB): \~ 0.9496

Run time: \~ 4 hours 48 minutes

``` r
Mt_data <- readr::read_table("Genotyping_pegard_etal_Alfalfa.txt")
#For demonstration purposes and to reduce computation load only chromsome 1 is imputed
idx <- grepl("chr1_", Mt_data$Chr_Pos)
Mt_geno_chr1 <- Mt_data[idx,-c(1:3)]
Mt_geno_chr1 <- t(Mt_geno_chr1)
colnames(Mt_geno_chr1) <- Mt_data$Chr_Pos[idx]

imputation_accuracy(geno = Mt_geno_chr1, prop_na = 0.4, repl = 1, k_neighbours = 20, l_loci = 10, cpus = 6, mem = 20, method = "ld_knni_fast")

# A tibble: 1 × 6
  Method       Prop_NA k_Neighbours l_Loci Run_ID     RB
  <chr>          <dbl>        <dbl>  <dbl>  <dbl>  <dbl>
1 ld_knni_fast     0.4           10     10      1 0.0514
```

If you want to speed up this demonstration you can subset chromosome 1 further, for example the first 10k loci took around 79 minutes on the same setup.

## A Note About Speed

resurface is by no means fast, there are other programming languages for that, but R is my jam.
For larger genotype datasets you would be wise to submit it as a long running job, and try and use multiple CPUs where possible.
That said you are unlikely to be imputing large genotypic datasets every second day, so resurface should hopefully get you out of a pickle (that's if you are in one).

Think of it like the old farm ute, its not fast but it gets the job done 👍

⚠️ However if you do have suggestions on speeding up the code, especially the costly correlation calculation please shoutout.
I would be excited to have contributions to this package.

## Citation

Until resurface is on CRAN and has its own separate DOI, please cite the following paper and this github repo.

> *Exploitation of data from breeding programs supports rapid implementation of genomic selection for key agronomic traits in perennial ryegrass.*
>
> **Luke W. Pembleton**, Courtney Inch, Rebecca C. Baillie, Michelle C. Drayton, Preeti Thakur, Yvonne O. Ogaji, German C. Spangenberg, John W. Forster, Hans D. Daetwyler, Noel O. I. Cogan
>
> Theoretical and Applied Genetics, 131, 1891-1902, September 2018.
> doi: [10.1007/s00122-018-3121-7](https://doi.org/10.1007/s00122-018-3121-7)

Luke W. Pembleton, Courtney Inch, Rebecca C. Baillie, Michelle C. Drayton, Preeti Thakur, Yvonne O. Ogaji, German C. Spangenberg, John W. Forster, Hans D. Daetwyler, Noel O. I. Cogan (2018).
Exploitation of data from breeding programs supports rapid implementation of genomic selection for key agronomic traits in perennial ryegrass.
Theoretical and Applied Genetics, 131(9).

## Acknowledgments

Inspiration for resurface came from the original implementation of LD-kNNi for the imputation of genotypic classes in the software [`LinkImpute`](http://www.cultivatingdiversity.org/software.html) by Daniel Money

> *LinkImpute: Fast and Accurate Genotype Imputation for Nonmodel Organisms*
>
> Daniel Money, Kyle Gardner, Zoë Migicovsky, Heidi Schwaninger, Gan-Yuan Zhong, Sean Myles
>
> G3 Genes\|Genomes\|Genetics, Volume 5, Issue 11, 1 November 2015, Pages 2383–2390, doi: [10.1534/g3.115.021667](https://doi.org/10.1534/g3.115.021667)
