#' Inherit parameters for functions
#'
#' Template function to inherit parameters for other functions.
#'
#' @param DNAm A matrix of methylation beta values with samples as rows and CpGs as columns. Must include row names (samples) and unique column names (CpGs).
#' @param pheno A data frame or matrix containing sample phenotype data, with samples as rows. Must include a column called `Sample_ID` that contains the ids of the samples to be matched with the beta matrix.
#' @param imputation Logical value that will allows you to perform (T)/ skip (F) imputation of mean values for missing CpGs. Warning: when imputation = F if there are missing CpGs, it will automatically ignore these CpGs during calculation, making the clock values less accurate.
#' @keywords internal
param_template <- function(DNAm, pheno, imputation) {}

#' @keywords internal
#' @noRd
head1 <- function(x, n = 5, m = 5) {x[seq_len(n), seq_len(m), drop = F]}
