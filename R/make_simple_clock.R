#' Create a simple linear epigenetic clock calculator
#'
#' @param obj List/data.frame with CpGs and coefficients
#' @param intercept Numeric intercept
#' @param clock Column name for the resulting clock
#' @param cpg Name of CpG vector element in `obj`
#' @param coefficients Name of coefficient vector element in `obj`
#'
#' @return Function that adds the clock to a phenotype data frame
#' @keywords internal
make_simple_clock <- function(obj, intercept, clock, cpg, coefficients) {
  stopifnot(all(c(cpg, coefficients) %in% names(obj)))
  stopifnot(length(obj[[cpg]]) == length(obj[[coefficients]]))

  # force evaluation so the closure is self-contained
  cpg_names <- obj[[cpg]]
  coefs <- setNames(obj[[coefficients]], cpg_names)
  intercept_ <- intercept
  clock_ <- clock

  function(DNAm, pheno = NULL) {
    check_DNAm(DNAm)
    if (is.null(pheno)) {
      ids <- rownames(DNAm)
      if (is.null(ids)) ids <- as.character(seq_len(nrow(DNAm)))
      pheno <- data.frame(Sample_ID = ids, stringsAsFactors = FALSE)
    }

    # Imputation
    CpGs <- numeric(length(cpg_names))
    names(CpGs) <- cpg_names
    DNAm <- impute_DNAm(DNAm = DNAm, method = "mean", CpGs = CpGs, subset = TRUE)

    # Realign DNAm columns and coefficients to a common order
    DNAm <- DNAm[, cpg_names, drop = FALSE]
    pheno[[clock_]] <- as.vector(DNAm %*% coefs[cpg_names]) + intercept_
    pheno
  }
}
