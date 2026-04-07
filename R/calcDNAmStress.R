#' Calculate DNAmStress
#'
#' @description Calculates the DNAmStress (MS Stress) methylation index,
#'   a DNA methylation-based predictor of cumulative stress exposure.
#'
#' @inheritParams param_template
#'
#' @inherit param_template return
#'
#' @references
#' Katrinli, S., et al. (2023). Epigenetic aging and perceived psychological
#' stress in old age. Biological Psychiatry, 93(5), 451-459.
#' \url{https://doi.org/10.1016/j.biopsych.2022.08.012}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' calcDNAmStress(exampleBetas, examplePheno)
#' }
calcDNAmStress <- function(DNAm, pheno = NULL) {
  check_DNAm(DNAm)

  cpg_names <- DNAmStress_CpGs$CpG

  DNAm <- impute_DNAm(
    DNAm = DNAm,
    method = "mean",
    CpGs = zero_cpgs(cpg_names),
    subset = TRUE
  )
  DNAm <- DNAm[, cpg_names, drop = FALSE]

  predictions <- as.numeric(
    as.matrix(DNAm) %*% DNAmStress_CpGs$coefficient + DNAmStress_ref$intercept
  )

  if (is.null(pheno)) {
    result <- data.frame(DNAmStress = predictions)
    row.names(result) <- row.names(DNAm)
    return(result)
  }
  pheno$DNAmStress <- predictions
  return(pheno)
}
