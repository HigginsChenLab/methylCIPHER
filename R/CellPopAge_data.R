#' @title CellPopAge CpGs
#'
#' @description The 42 CpG probes with non-zero coefficients in the CellPopAge clock,
#'   used by \code{getClockProbes} to report probe coverage.
#'
#' @format A data frame with 42 rows and 2 variables:
#' \describe{
#'   \item{CpG}{CpG probe ID}
#'   \item{coefficient}{Non-zero model coefficient}
#' }
#'
#' @source \url{https://github.com/ucl-medical-genomics/CellPopAge-epigenetic-clock}
"CellPopAge_CpGs"

#' @title CellPopAge Reference Data
#'
#' @description Intercept and full coefficient set for the CellPopAge clock.
#'
#' @format A list with two elements:
#' \describe{
#'   \item{intercept}{Numeric scalar intercept (11.527).}
#'   \item{all_CpGs}{Data frame with all 2,543 CpGs and their coefficients
#'     (including zeros), used for calculation.}
#' }
"CellPopAge_ref"
