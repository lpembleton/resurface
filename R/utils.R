#' Weighted Variance Calculation
#'
#' @param x Vector of allele frequencies.
#' @param w Genotype distance value ^-1.
#' @return Weighted variance
#' @export
weighted_var <- function(x, w) {
  sum_w <- sum(w)
  sum_w2 <- sum(w^2)
  mean_w <- sum(x * w) / sum(w)
  (sum_w / (sum_w^2 - sum_w2)) * sum(w * (x - mean_w)^2, na.rm = TRUE)
}

#' Parse Genotype Matrix
#'
#' Parse genotype matrix out of input format and store sample IDs if present.
#'
#' @param geno Either a matrix of allele frequencies with sample IDs as row names,
#'             or a 2-column tibble with samples IDs in Col1 and a matrix-column
#'             of allele frequencies in Col2. Samples across rows, Loci down columns.
#' @return A list containing a vector of sample names as the first element and a
#'         genotype matrix as the second element.
#' @import dplyr tibble
#' @export
parse_geno <- function(geno) {
  if (is.matrix(geno)) {
    return(list(IDs = row.names(geno), geno = geno, intype = "matrix"))
  }
  if (tibble::is_tibble(geno)) {
    return(list(IDs = geno[, 1], geno = dplyr::pull(geno[, 2]), intype = "tibble"))
  }
}

#' Check Genotype Validity
#'
#' Function to check the validity of input geno and also fill any homozygous loci.
#' Note: as there is no variance for homozygous loci to correlate with, imputation
#' isnt possible and therefore the whole loci is recoded as as homozygous
#'
#' @param geno Matrix of allele frequency genotype data.
#'             Samples across rows, Loci down columns.
#' @return Genotype matrix with any homozygous loci prefilled.
#' @import cli
#' @export
check_geno <- function(geno) {
  if (!is.matrix(geno)) {
    cli::cli_abort("Input geno must be a matrix")
  }
  if (!is.numeric(geno)) {
    cli::cli_abort("Input geno must be numeric")
  }
  if (max(geno, na.rm = T) > 1 | min(geno, na.rm = T) < 0) {
    cli::cli_abort("Allele frequencies must be between 0 and 1")
  }

  v <- apply(geno, 2, var, na.rm = T)
  idx <- which(v == 0)
  if (length(idx) != 0) {
    for (i in idx) {
      mono_AF <- unique(geno[, i] |> (\(x) x[!is.na(x)])())
      geno[, i] <- mono_AF
    }
  }
  geno
}


#' Mean Allele Frequency Impute
#'
#' Function to impute missing data with just the mean allele frequency at the locus
#'
#' @param geno Matrix of allele frequency genotype data.
#'             Samples across rows, Loci down columns.
#' @return Genotype matrix imputed using the basic mean allele frequency method.
#' @import cli
#' @export
mean_freq_impute <- function(geno){

  incomplete_loci <- which(colSums(is.na(geno)) > 0) # index of loci that contain missing values
  for(i in incomplete_loci){
    na_idx <- which(is.na(geno[,i]))
    mean_freq <- mean(geno[,i], na.rm = T)
    geno[na_idx,i] <- mean_freq
  }
  geno

}

