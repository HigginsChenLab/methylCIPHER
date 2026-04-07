#' Calculate CellPopAge
#'
#' @description Calculates CellPopAge, an epigenetic clock trained to predict
#'   cumulative population doublings from DNA methylation data.
#'
#' @inheritParams param_template
#'
#' @inherit param_template return
#'
#' @references
#' Lowe, D., et al. (2022). CellPopAge: an epigenetic clock for replicative
#' aging of cell populations.
#' \url{https://github.com/ucl-medical-genomics/CellPopAge-epigenetic-clock}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' calcCellPopAge(exampleBetas, examplePheno)
#' }
calcCellPopAge <- function(DNAm, pheno = NULL) {
  check_DNAm(DNAm)

  all_cpgs <- CellPopAge_ref$all_CpGs
  cpg_names <- all_cpgs$CpG
  intercept <- CellPopAge_ref$intercept

  # Imputation: fill missing CpGs with 0 (matching original behavior)
  DNAm <- impute_DNAm(
    DNAm = DNAm,
    method = "mean",
    CpGs = zero_cpgs(cpg_names),
    subset = TRUE
  )
  DNAm <- DNAm[, cpg_names, drop = FALSE]

  predictions <- as.numeric(
    as.matrix(DNAm) %*% all_cpgs$coefficient + intercept
  )

  if (is.null(pheno)) {
    result <- data.frame(CellPopAge = predictions)
    row.names(result) <- row.names(DNAm)
    return(result)
  }
  pheno$CellPopAge <- predictions
  return(pheno)
}
