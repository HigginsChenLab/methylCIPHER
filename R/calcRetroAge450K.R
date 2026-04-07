#' RetroAge450K
#'
#' @description A function to calculate RetroAge450K
#'
#' @inheritParams param_template
#'
#' @inherit param_template return
#'
#' @export
#'
#' @examples
#' \dontrun{
#' calcRetroAge450K(exampleBetas, examplePheno)
#' }
calcRetroAge450K <- function(DNAm, pheno = NULL) {
  # Input validation
  check_DNAm(DNAm)

  # Imputation
  DNAm <- impute_DNAm(
    DNAm = DNAm,
    method = "mean",
    CpGs = zero_cpgs(RetroAge450K_CpGs$name),
    subset = TRUE
  )

  # Re-align to make sure things lined up with the object
  DNAm <- DNAm[, RetroAge450K_CpGs$name, drop = FALSE]

  # Calculate the epigenetic age using the dot product of the coefficients and the methylation data
  predictions <- as.numeric(
    as.matrix(DNAm) %*% RetroAge450K_CpGs$coefficient + 84.2097786
  )

  if (is.null(pheno)) {
    result <- data.frame(RetroAge450K = predictions)
    row.names(result) <- row.names(DNAm)
    return(result)
  }
  pheno$RetroAge450K <- predictions
  return(pheno)
}
