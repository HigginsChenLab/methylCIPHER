#' Inherit parameters for functions
#'
#' Template function to inherit parameters for other functions.
#'
#' @param DNAm A matrix of methylation beta values with samples as rows and CpGs as columns. Must include row names (samples) and unique column names (CpGs).
#' @param pheno Optional: A data frame or matrix containing sample phenotype data, with samples as rows. The clock results will be appended to this data. Defaults to \code{NULL}.
#' @param imputation Logical value that will allows you to perform (T)/ skip (F) imputation of mean values for missing CpGs. Warning: when imputation = F if there are missing CpGs, it will automatically ignore these CpGs during calculation, making the clock values less accurate.
#' @keywords internal
param_template <- function(DNAm, pheno, imputation) {}

#' @keywords internal
#' @noRd
head1 <- function(x, n = 5, k = 6) {head(x, n = n, k = k)}
