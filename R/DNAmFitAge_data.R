#' @title DNAmFitAge Model Data
#'
#' @description Bundled model coefficients, sex-specific imputation medians,
#'   and CpG lists for calculating DNAmFitAge fitness biomarkers.
#'
#' @format A list containing 14 elements:
#' \describe{
#'   \item{Gait_noAge_Females}{Female gait speed model (age-unadjusted)}
#'   \item{Gait_noAge_Males}{Male gait speed model (age-unadjusted)}
#'   \item{Gait_wAge_Females}{Female gait speed model (age-adjusted)}
#'   \item{Gait_wAge_Males}{Male gait speed model (age-adjusted)}
#'   \item{Grip_noAge_Females}{Female grip strength model (age-unadjusted)}
#'   \item{Grip_noAge_Males}{Male grip strength model (age-unadjusted)}
#'   \item{Grip_wAge_Females}{Female grip strength model (age-adjusted)}
#'   \item{Grip_wAge_Males}{Male grip strength model (age-adjusted)}
#'   \item{FEV1_wAge_Females}{Female FEV1 model (age-adjusted)}
#'   \item{FEV1_wAge_Males}{Male FEV1 model (age-adjusted)}
#'   \item{VO2maxModel}{VO2max model (both sexes)}
#'   \item{Female_Medians_All}{Sex-specific training medians for females}
#'   \item{Male_Medians_All}{Sex-specific training medians for males}
#'   \item{AllCpGs}{Character vector of all 627 required CpG probes}
#' }
#'
#' @source \url{https://github.com/kristenmcgreevy/DNAmFitAge}
"DNAmFitAge_data"

#' @title DNAmFitAge CpGs
#'
#' @description The CpG probes required by the DNAmFitAge clock.
#'
#' @format A data frame with 627 rows and 1 variable:
#' \describe{
#'   \item{CpG}{CpG probe ID}
#' }
#'
#' @source \url{https://github.com/kristenmcgreevy/DNAmFitAge}
"DNAmFitAge_CpGs"
