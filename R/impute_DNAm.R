#' Impute Missing Values and Add Missing CpGs in DNA Methylation Matrix
#'
#' This function imputes missing values in a DNA methylation matrix and adds completely missing CpG sites using values provided in a named vector.
#'
#' @inheritParams param_template
#' @param method The method used to impute missing values within existing CpG sites. Currently, only "mean" is supported, which imputes missing values with the column mean.
#' @param CpGs A named numeric vector where names are CpGs names and values are used to impute completely missing CpG sites across all samples. If `NULL`, no new CpG sites are added, and only missing values in existing CpG sites are imputed.
#' @param subset A boolean. If `TRUE`, then only the CpGs provided in the `CpGs` argument are imputed and returned. Otherwise, all CpGs that have missing values in `DNAm` will be imputed and the full DNAm matrix with missing CpGs are returned.
#'
#' @details
#' The function handles two types of missing data in a DNA methylation beta matrix: 1) completely missing CpGs are filled in with the provided values and 2) missing values of existing CpGs are filled in by supported method (e.g, mean imputation).
#' If `CpGs` is `NULL` only missing values in the existing columns are imputed.
#'
#' @return An imputed numeric matrix with no missing values
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Sample DNAm matrix with missing values
#' DNAm <- matrix(
#'   c(0.1, 0.2, NA, 0.4, 0.5, 0.6),
#'   nrow = 2,
#'   dimnames = list(c("sample1", "sample2"), c("cpg1", "cpg2", "cpg3"))
#' )
#' # Named CpGs vector where names are CpGs and values are used to impute completely missing probes
#' CpGs <- c(cpg1 = 0.3, cpg2 = 0.3, cpg3 = 0.3, cpg4 = 0.7)
#' # Impute the matrix
#' imputed_DNAm <- impute_DNAm(DNAm, method = "mean", CpGs = CpGs)
#' imputed_DNAm
#' }
impute_DNAm <- function(DNAm, method = c("mean"), CpGs = NULL, subset = TRUE) {
  # Input validation
  suppressWarnings(check_DNAm(DNAm))
  method <- match.arg(method)
  checkmate::assert_numeric(CpGs, min.len = 1, names = "unique", null.ok = TRUE)
  checkmate::assert_logical(subset, len = 1, any.missing = FALSE, null.ok = FALSE)

  if(is.null(CpGs)) {
    CpGs <- numeric(length = ncol(DNAm))
    names(CpGs) <- colnames(DNAm)
  }

  # Completely missing CpG
  needed_cpgs <- setdiff(names(CpGs), colnames(DNAm))

  needed_matrix <- matrix(
    CpGs[needed_cpgs],
    ncol = length(needed_cpgs),
    nrow = nrow(DNAm),
    byrow = TRUE,
    dimnames = list(row.names(DNAm), needed_cpgs)
  )

  # Imputed CpG
  to_be_imputed <- if(subset) {
    intersect(names(CpGs), colnames(DNAm))
  } else {
    colnames(DNAm)
  }
  n_miss <- colSums(is.na(DNAm[, to_be_imputed, drop = F]))

  imputed_matrix <- if(sum(n_miss) == 0) {
    DNAm[, to_be_imputed, drop = F]
  } else if(method == "mean") {
    cbind(
      DNAm[, names(which(n_miss == 0)), drop = F],
      mean_impute(DNAm[, names(which(n_miss > 0)), drop = F])
    )
  } else {
    stop("Unsupported Method")
  }

  return(cbind(imputed_matrix, needed_matrix)[, c(to_be_imputed, needed_cpgs), drop = F])
}
