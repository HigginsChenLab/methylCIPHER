#' impute DNAm matrix
#'
#' @inheritParams param_template
#' @param method which method to impute missing CpGs
#' @param CpGs vector of CpGs expected from the matrix. Vector of values of CpGs to be filled in with the CpG names as names.
#'
#' @returns An imputed matrix with just the CpGs needed for the calculation of the clock
# impute_DNAm <- function(DNAm, method = c("mean"), CpGs) {
#   check_DNAm(DNAm)
#   method <- match.arg(method)
#
#   if(sum(is.na(DNAm)) == 0) {
#     return(DNAm)
#   }
# }
