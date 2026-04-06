#' RetroAgeEPICv2
#'
#' @description A function to calculate RetroAgeEPICv2
#'
#' @inheritParams param_template
#'
#' @inherit param_template return
#'
#' @export
#'
#' @examples
#' \dontrun{
#' calcRetroAgeEPICv2(exampleBetas, examplePheno)
#' }
calcRetroAgeEPICv2 <- function(DNAm, pheno = NULL) {
  # Input validation
  check_DNAm(DNAm)

  # Build imputation reference from bundled data, falling back to 0
  impute_vals <- RetroAgeEPICv2_ref$imputation
  cpg_names <- RetroAgeEPICv2_CpGs$name
  if (is.null(impute_vals)) {
    impute_vals <- zero_cpgs(cpg_names)
  } else {
    missing_from_ref <- setdiff(cpg_names, names(impute_vals))
    if (length(missing_from_ref) > 0) {
      impute_vals <- c(impute_vals, zero_cpgs(missing_from_ref))
    }
  }
  impute_vals <- impute_vals[cpg_names]

  # Imputation
  DNAm <- impute_DNAm(
    DNAm = DNAm,
    method = "mean",
    CpGs = impute_vals,
    subset = TRUE
  )

  # Re-align to make sure things lined up with the object
  DNAm <- DNAm[, cpg_names, drop = FALSE]

  # Calculate the epigenetic age using the dot product of the coefficients and the methylation data
  predictions <- as.numeric(
    as.matrix(DNAm) %*% RetroAgeEPICv2_CpGs$coefficient + RetroAgeEPICv2_ref$intercept
  )

  if (is.null(pheno)) {
    result <- data.frame(RetroAgeEPICv2 = predictions)
    row.names(result) <- row.names(DNAm)
    return(result)
  }
  pheno$RetroAgeEPICv2 <- predictions
  return(pheno)
}
