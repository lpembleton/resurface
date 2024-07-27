#' Simulated Missing Data
#'
#' This function masks a proportion of data as missing. Useful for imputation
#' accuracy testing purposes
#'
#' @param geno Matrix of allele frequency genotype data. Samples across
#'             rows, Loci down columns.
#' @param prop_na Proportion of missing data to be introduced, in addition
#'                to any already existing
#' @return Matrix of genotype data with missing data added as NAs.
#' @import dplyr tibble
#' @export
sim_na_geno <- function(geno, prop_na = 0.2) {
  # get all indexes of geno matrix
  geno_idx <- tibble::tibble(
    row = rep(1:nrow(geno), times = ncol(geno)),
    col = rep(1:ncol(geno), each = nrow(geno))
  )

  # get all indexes of already present missing data
  na_idx <- tibble::as_tibble(which(is.na(geno), arr.ind = TRUE))

  # retain only those indexes with no missing data
  non_na_idx <- dplyr::bind_rows(geno_idx, na_idx) |>
    dplyr::group_by_all() |>
    dplyr::count() |>
    dplyr::filter(n == 1) |>
    dplyr::select(-n)

  # sample indexes for simulated missing data at prop_na rate
  sim_na_idx <- non_na_idx[sample(1:nrow(non_na_idx), ceiling(prop_na * nrow(geno_idx))), ]

  # overwrite missing data into genotype matrix and return
  geno[as.matrix(sim_na_idx)] <- NA
  geno
}

#' Calculate Imputation Accuracy
#'
#' This function calculates the imputation accuracy for missing data using
#' k-nearest neighbors (KNN) approach by masking a proportion of data as
#' missing. It can also be used to explore optimal k and l values
#'
#' @param geno Either a matrix of allele frequencies with sample IDs as row names,
#'             or a 2-column tibble with samples IDs in Col1 and a matrix-column
#'             of allele frequencies in Col2. Samples across rows, Loci down columns.
#' @param prop_na The proportion of data to mask as missing in the dataset for
#'                accuracy testing. Default is 0.2.
#' @param repl Number of replicate runs for imputation accuracy calculation.
#'             Default is 1.
#' @param k_neighbours The number of k-neighbours (samples) to use in the KNN approach.
#' @param l_loci Number of loci, ordered based on correlation to the SNP to be imputed, to use when calculating nearest neighbours.
#' @param cpus The number of CPUs to use for parallel processing. Default is 1.
#' @param mem Memory availability in GiB. This is divided by cpus for per worker allocation (default: 1)
#' @param method Imputation method to use. Currently the options are "ld_knni" or "ld_knni_fast" (default: "ld_knni_fast"). When using the ld_knni_fast method missing data will be temporarily replaced with the mean loci allele frequency for correlation calculations (correlation calculation is faster on a complete matrix). This can reduce imputation accuracy, especially in datasets with high missing data, while the impact is smaller in datasets with lower levels of missing data.
#' @return A tibble containing imputation statistics, currently raw bias (RB) for each replicate run.
#' @import dplyr tibble
#' @export
imputation_accuracy <- function(geno, prop_na = 0.2, repl = 1, k_neighbours = 10, l_loci = 10, cpus = 1, mem = 1, method = "ld_knni_fast") {
  # TODO: add option to live report imputation accuracy as it steps through each locus. This might allow you to determine optimal l and k quicker

  # Check method is valid
  if( !(method %in% c("ld_knni", "ld_knni_fast", "mean_freq")) ){
    cli::cli_abort("method does not exist")
  }

  data <- parse_geno(geno)
  geno <- check_geno(data$geno)

  accuracy_boots <- tibble::tibble(
    Method = character(),
    Prop_NA = numeric(),
    k_Neighbours = numeric(),
    l_Loci = numeric(),
    Run_ID = numeric(),
    RB = numeric()
  )

  for (b in 1:repl) {

    geno_na <- sim_na_geno(geno = geno, prop_na = prop_na) #add simulated missing data

    if (method == "ld_knni") { # Impute using ld-knni method
      geno_imp <- matrix_knni(
        geno = geno_na,
        k = k_neighbours,
        l = l_loci,
        cpus = cpus,
        mem = mem,
        fast = FALSE
      )
    }

    if (method == "ld_knni_fast") { # Impute using ld-knni fast method
      geno_imp <- matrix_knni(
        geno = geno_na,
        k = k_neighbours,
        l = l_loci,
        cpus = cpus,
        mem = mem,
        fast = TRUE
      )
    }

    if (method == "mean_freq") { # Impute using locus mean allele freq method (currently this is for internal dev/testing purposes)
      geno_imp <- mean_freq_impute(geno_na)
    }

    # Calculate accuracy statistics

    na_idx <- which(is.na(geno_na) & !is.na(geno), arr.ind = TRUE)
    raw_bias <- mean(abs(geno_imp[na_idx] - geno[na_idx])) #raw bias

    accuracy_boots <- accuracy_boots |>
      dplyr::add_row(
        Method = method,
        Prop_NA = prop_na,
        k_Neighbours = k_neighbours,
        l_Loci = l_loci,
        Run_ID = b,
        RB = raw_bias
      )

  }

  accuracy_boots
}

