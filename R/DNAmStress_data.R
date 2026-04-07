#' @title DNAmStress CpGs
#'
#' @description CpG probes and coefficients for the DNAmStress (MS Stress)
#'   methylation index.
#'
#' @format A data frame with 211 rows (CpGs) and 2 variables:
#' \describe{
#'   \item{CpG}{CpG probe ID}
#'   \item{coefficient}{Model coefficient}
#' }
#'
#' @source Katrinli et al. (2023) Biological Psychiatry, Table S10.
"DNAmStress_CpGs"

#' @title DNAmStress Reference Data
#'
#' @description Intercept for the DNAmStress clock.
#'
#' @format A list with one element:
#' \describe{
#'   \item{intercept}{Numeric scalar intercept (-4.494069).}
#' }
"DNAmStress_ref"
