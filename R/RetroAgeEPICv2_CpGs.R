#' @title RetroAgeEPICv2 Weights
#'
#' @description The CpGs and coefficients for calculating RetroAgeEPICv2
#'
#' @format A data frame with 1378 rows (CpGs) and 2 variables:
#' \describe{
#'   \item{name}{The Illumina ID of the CpG included}
#'   \item{coefficient}{The coefficient weight for calculating RetroAgeEPICv2}
#' }
"RetroAgeEPICv2_CpGs"

#' @title RetroAgeEPICv2 Reference Data
#'
#' @description Intercept and imputation reference for RetroAgeEPICv2
#'
#' @format A list with two elements:
#' \describe{
#'   \item{intercept}{Numeric scalar intercept for the linear model.}
#'   \item{imputation}{Named numeric vector for imputing missing CpGs.}
#' }
"RetroAgeEPICv2_ref"
