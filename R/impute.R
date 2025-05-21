#' Imputation of a Single Missing Allele Frequency
#'
#' This is an internal function, please use impute() for imputation.
#' This function uses knni logic to impute an allele frequency for a single
#' datapoint.
#'
#' @param k Number of nearest (samples) to use.
#' @param l Number of closest loci to use.
#' @param geno Matrix of allele frequency genotypes, SNPs in columns,
#'             samples across rows.
#' @param locus Column index value for the locus that need imputing.
#' @param ld_order Vector containing a ordered index from pairwise LD/
#'                 correlation values between locus and the geno matrix.
#' @param ind Row index val for the sample needing imputing.
#' @param prop_l due to missing data using l loci might not always be possible. prop_l is the the min prop of l loci required, or max non na loci will be used.
#' @return An imputed genotype frequency.
#' @import cli
#' @importFrom stats weighted.mean rnorm
#' @export
singlepoint_knni <- function(k = 20, l = 10, locus, geno, ld_order, ind, prop_l = 0.5) {
  l_loci_matrix <- geno[, ld_order[1:l]] # gather the l closest loci based on LD

  # TODO: add variable l selection, for example loci with cor values less then x are excluded even if this means the num of loci is <l

  while ((sum(!is.na(l_loci_matrix[ind, ])) / l) < prop_l & l < ncol(geno)) {
    l <- l + 1 # increment l loci by 1 until min prop_l is met
    l_loci_matrix <- geno[, ld_order[1:l]]
  }

  if (sum(!is.na(l_loci_matrix[ind, ])) == 0) {
    cli::cli_abort("Cannot continue, insufficient data for individual ", ind, " locus ", locus, " to able to imputed")
  }

  ind_repl_matrix <- t(replicate(nrow(geno), l_loci_matrix[ind, ])) # replicated matrix of ind geno for l loci
  k_dist <- rowMeans(abs(l_loci_matrix - ind_repl_matrix), na.rm = TRUE) # distance or difference between individuals geno at k loci
  valid_idx <- which(!is.na(geno[, locus]))

  if (length(valid_idx) < 0) {
    cli::cli_abort("Cannot continue, insufficient neigbours with non missing data to perform imputation")
  } else {
    if (length(valid_idx) < k) {
      k_int <- length(valid_idx)
    } else {
      k_int <- k
    }
  }
  k_idx <- order(k_dist) |> (\(x) x[x %in% valid_idx][1:k_int])()

  d <- k_dist[k_idx]
  d[d == 0] <- 0.0001 # prevent weights becoming Inf in later calculations
  d[is.na(d)] <- 0.0001 # those that had no pairwise comparisions of non missing genotype should have their NA dist recoded to 0.0001 to prevent downstream error 

  weighted_mu <- weighted.mean(x = geno[k_idx, locus], w = d^-1, na.rm = TRUE)
  weighted_sd <- sqrt(weighted_var(x = geno[k_idx, locus], w = d^-1))
  i_af <- rnorm(1, mean = weighted_mu, sd = weighted_sd) # imputed allele freq
  # Note: Because any rnorm value beyond the theoretical limits of 0 or 1 provides support for the limit, rnorm doesn't need to be rerun until it is between the upper or lower bound. Instead it can just be forced to the boundary value.
  if (i_af > 1) {
    i_af <- 1
  } else {
    if (i_af < 0) {
      i_af <- 0
    }
  }
  i_af
}



#' Imputation of a Matrix of Missing Allele Frequencies
#'
#' This is an internal function, please use impute() for imputation.
#' This function uses knni logic to impute all missing allele frequency
#' datapoints withing a matrix.
#'
#' @param geno Matrix of allele frequency genotype data to be imputed.
#'             Samples across rows, Loci down columns.
#' @param k Number of nearest neighbours (samples) to use.
#' @param l Number of closet loci to use.
#' @param cpus Number of cpus or threads to use for parallel processing.
#' @param mem Memory availability in GiB. This is divided by cpus for per worker allocation (default: 1).
#' @param fast TRUE/FALSE use a faster method but with an accuracy reduction.
#' @return Matrix of allele frequency genotype data containing imputed
#'         genotype frequencies.
#' @import foreach future doFuture dplyr zoo progressr
#' @importFrom stats cor
#' @export
matrix_knni <- function(geno, k = 10, l = 6, cpus = 1, mem = 1, fast = FALSE) {
  future::plan(future::multisession, workers = cpus)

  progressr::handlers(global = TRUE)
  progressr::handlers("cli")

  oopts <- options(future.globals.maxSize = 500 * 1024^2) ## 500 MB
  on.exit(options(oopts)) # when exiting reset future maxSize to 500MB

  options(future.globals.maxSize = ((mem / cpus) * (1000 * 1024^2)))

  if (fast == TRUE) {
    cor_mat <- abs(cor(zoo::na.aggregate(geno)))

    incomplete_loci <- which(colSums(is.na(geno)) > 0) # index of loci that contain missing values
    p <- progressr::progressor(along = incomplete_loci) # setup progress bar reporting
    indexed_iAF <- foreach::foreach(locus_i = incomplete_loci, .options.future = list(seed = TRUE)) %dofuture% {
      na_idx <- which(is.na(geno[, locus_i]))
      cor_order <- order(
        cor_mat[, locus_i],
        decreasing = TRUE
      )
      imputed_AF <- tibble::tibble(
        row = na_idx,
        column = locus_i,
        iAF = as.numeric(NA)
      )
      for (s in na_idx) {
        imputed_AF$iAF[imputed_AF$row == s] <- singlepoint_knni(
          k = k,
          l = l,
          locus = locus_i,
          geno = geno,
          ld_order = cor_order,
          ind = s
        )
      }
      p(sprintf("locus_i=%g", locus_i))
      imputed_AF
    }
  } else {
    incomplete_loci <- which(colSums(is.na(geno)) > 0) # index of loci that contain missing values
    p <- progressr::progressor(along = incomplete_loci) # setup progress bar reporting
    indexed_iAF <- foreach::foreach(locus_i = incomplete_loci, .options.future = list(seed = TRUE)) %dofuture% {
      na_idx <- which(is.na(geno[, locus_i]))
      cor_order <- suppressWarnings(cor(geno[, locus_i], geno, use = "pairwise.complete")) |>
        abs() |>
        order(na.last = TRUE, decreasing = TRUE)
      imputed_AF <- tibble::tibble(
        row = na_idx,
        column = locus_i,
        iAF = as.numeric(NA)
      )
      for (s in na_idx) {
        imputed_AF$iAF[imputed_AF$row == s] <- singlepoint_knni(
          k = k,
          l = l,
          locus = locus_i,
          geno = geno,
          ld_order = cor_order,
          ind = s
        )
      }
      p(sprintf("locus_i=%g", locus_i))
      imputed_AF
    }
  }

  indexed_iAF <- dplyr::bind_rows(indexed_iAF)
  geno[as.matrix(indexed_iAF[, 1:2])] <- indexed_iAF$iAF
  geno
}


#' Impute Missing Allele Frequencies
#'
#' This function imputes missing allele frequency data using the
#' linkage disequilibrium k-nearest neighbour imputation algorithm
#' for allele frequencies.
#'
#' @param geno Either a matrix of allele frequencies with sample IDs as row names,
#'             or a 2-column tibble with samples IDs in Col1 and a matrix-column
#'             of allele frequencies in Col2. Samples across rows, Loci down columns.
#' @param k_neighbours The number of k-neighbours (samples) to use in the KNN approach (default: 10).
#' @param l_loci Number of loci, ordered based on correlation to the SNP to be impute, to use when calculating nearest neighbours (default: 10).
#' @param cpus The number of CPUs to use for parallel processing (default: 1).
#' @param mem Memory availability in GiB. This is divided by cpus for per worker allocation (default: 1)
#' @param method Imputation method to use. Currently the options are "ld_knni" or "ld_knni_fast" (default: "ld_knni_fast"). When using the ld_knni_fast method missing data will be temporarily replaced with the mean loci allele frequency for correlation calculations (correlation calculation is faster on a complete matrix). This can reduce imputation accuracy, especially in datasets with high missing data, while the impact is smaller in datasets with lower levels of missing data.
#' @return Depending on the input geno format, either a matrix of imputed allele
#'         frequencies with sample IDs as the row names, or a 2-column tibble with
#'         samples IDs in Col1 and a matrix-column of imputed allele frequencies
#'         in Col2.
#' @import dplyr tibble
#' @export
impute <- function(geno, k_neighbours = 10, l_loci = 10, cpus = 1, mem = 1, method = "ld_knni_fast") {

  # Check method
  if( !(method %in% c("ld_knni", "ld_knni_fast")) ){
    cli::cli_abort("method does not exist")
  }

  data <- parse_geno(geno)
  geno <- check_geno(data$geno)

  if( method == "ld_knni" ){
    geno_imp <- matrix_knni(
      geno = geno,
      k = k_neighbours,
      l = l_loci,
      cpus = cpus,
      mem = mem,
      fast = FALSE
    )
  }

  if( method == "ld_knni_fast" ){
    geno_imp <- matrix_knni(
      geno = geno,
      k = k_neighbours,
      l = l_loci,
      cpus = cpus,
      mem = mem,
      fast = TRUE
    )
  }


  if (data$intype == "matrix") {
    cli::cli_alert_info("matrix type")
    row.names(geno_imp) <- data$IDs
    return(geno_imp)
  }
  if (data$intype == "tibble") {
    cli::cli_alert_info("tibble type")
    tibble::tibble(IDs = data$IDs, geno = geno_imp)
    return(geno_imp)
  }
}

