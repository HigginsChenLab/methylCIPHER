#' Inherit parameters for functions
#'
#' Template function to inherit parameters for other functions.
#'
#' @param DNAm A matrix of methylation beta values with samples as rows and CpGs as columns. Must include row names (samples) and unique column names (CpGs).
#' @param pheno A data.frame containing sample phenotype data, with samples as rows. Must include the following columns with no missing values:
#' \itemize{
#'   \item \code{Sample_ID}: Character vector of sample IDs to match with `DNAm`.
#' }
#' @param ID Column in `pheno` that indicates the ID of each sample.
#' @param imputation Logical value that will allows you to perform (T)/ skip (F) imputation of mean values for missing CpGs. Warning: when imputation = F if there are missing CpGs, it will automatically ignore these CpGs during calculation, making the clock values less accurate.
#'
#' @return If you added the `pheno` input (preferred) the function appends a column with the clock calculation and returns the data.frame. Otherwise, it will return a data.frame of calculated clock values.
#'
#' @keywords internal
param_template <- function(DNAm, pheno, ID, imputation) {}

#' Inherit parameters for functions. With Female and Age
#'
#' Template function to inherit parameters for other functions requiring phenotype data.
#'
#' @inheritParams param_template
#' @param pheno A data.frame containing sample phenotype data, with samples as rows. Must include the following columns with no missing values:
#' \itemize{
#'   \item \code{Sample_ID}: Character vector of sample IDs to match with `DNAm`.
#'   \item \code{Female}: Numeric vector of 1 (female) or 0 (non-female).
#'   \item \code{Age}: Numeric vector of chronological ages (finite values).
#' }
#'
#' @return A data.frame of calculated clock values appended to `pheno`.
#'
#' @keywords internal
param_template_female_age <- function(DNAm, pheno, ID, imputation) {}

#' @keywords internal
#' @noRd
head1 <- function(x, n = 5, m = 5) {
  x[seq_len(min(n, nrow(x))), seq_len(min(m, ncol(x))), drop = FALSE]
}

#' Align Pheno
#'
#' Filter and sort pheno given a list of ids. ids must be in pheno.
#'
#' @inheritParams param_template
#'
#' @examples
#' \dontrun{
#' align_pheno(examplePheno, c("7786915023_R02C02", "7786915023_R02C02"))
#' }
#'
#' @keywords internal
align_pheno <- function(pheno, Sample_ID, ID) {
  indices <- match(Sample_ID, pheno[[ID]])
  return(pheno[indices, , drop = F])
}

#' Convert to Vector of Zeros with CpG Names
#'
#' Convert a vector of CpGs to a numeric vector of zeros, with the names of the vector set to the CpGs.
#'
#' @inheritParams param_template
#'
#' @return A numeric vector of zeros with names set to the input CpGs.
#'
#' @keywords internal
zero_cpgs <- function(CpGs) {
  z <- numeric(length(CpGs))
  names(z) <- CpGs
  z
}

#' Remove Columns with All NA Values from DNAm
#'
#' Remove columns from the DNAm matrix where all values are NA.
#'
#' @inheritParams param_template
#'
#' @return DNAm with columns removed where all values are NA.
#'
#' @keywords internal
removeNAcol <- function(DNAm) {
  na <- is.na(DNAm)
  keep <- (colSums(na) < nrow(DNAm))
  return(DNAm[, keep])
}

#' trafo
#'
#' @param x A vector of sample ages
#' @param adult.age Age boundary of adulthood--set to 20
#'
#' @return a vector of transformed ages
trafo <- function(x, adult.age = 20) {
  x <- (x + 1) / (1 + adult.age)
  y <- ifelse(x <= 1, log(x), x - 1)
  y
}

#' Developmental Age Transformation
#'
#' @param x A vector of sample ages
#' @param adult.age The age considered to be the cutoff for adulthood
#'
#' @return transformed age prediction
anti.trafo <- function(x, adult.age = 20) {
  ifelse(x < 0, (1 + adult.age) * exp(x) - 1, (1 + adult.age) * x + adult.age)
}
