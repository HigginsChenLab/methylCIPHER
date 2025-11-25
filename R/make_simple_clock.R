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
  function(DNAm, pheno = NULL, ID = "Sample_ID") {
    # Input validation
    # Check DNAm
    check_DNAm(DNAm)
    # Check Pheno
    if (is.null(pheno)) {
      pheno <- data.frame(Sample_ID = row.names(DNAm))
      ID <- "Sample_ID"
    }
    check_pheno(pheno, ID = ID)
    pheno_cols <- names(pheno)
    # Check Consistent between `pheno` and `DNAm`
    need_align <- !isTRUE(all.equal(row.names(DNAm), pheno[[ID]]))
    if (need_align) {
      samples <- intersect(row.names(DNAm), pheno[[ID]])
      if (length(samples) == 0) {
        stop("DNAm and pheno have no ID in common.")
      }
      DNAm <- DNAm[samples, , drop = FALSE]
      pheno <- align_pheno(pheno, samples, ID = ID)
      stopifnot("`DNAm` and `pheno` samples alignment failed. Check ID of pheno and row.names() of `DNAm`" = isTRUE(all.equal(row.names(DNAm), pheno[[ID]])))
      message("Samples inconsistencies between DNAm and Pheno were detected and corrected.")
    }

    ## Imputation
    CpGs <- numeric(length = length(obj[[cpg]]))
    names(CpGs) <- obj[[cpg]]

    DNAm <- impute_DNAm(
      DNAm = DNAm,
      method = "mean",
      CpGs = CpGs,
      subset = TRUE
    )

    ## Re-align to make sure things lined up with the object
    DNAm <- DNAm[, obj[[cpg]], drop = F]
    pheno[[clock]] <- as.vector(DNAm %*% obj[[coefficients]]) + intercept
    return(pheno)
  }
}
