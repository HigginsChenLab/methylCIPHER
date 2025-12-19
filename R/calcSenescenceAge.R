#' SenChronoAge
#'
#' A function to calculate SenChronoAge
#'
#' @inheritParams param_template
#' @inherit param_template return
#'
#' @export
calcSenChronoAge <- function(DNAm, pheno = NULL, ID = "Sample_ID") {
  make_simple_clock(
    obj = SenChronoAge_CpGs,
    intercept = -8.005320e+01,
    clock = "SenChronoAge",
    cpg = "CpG",
    coefficients = "Coefficient"
  )(DNAm, pheno, ID)
}

#' SenCultureAge
#'
#' A function to calculate SenCultureAge
#'
#' @inheritParams param_template
#' @inherit param_template return
#' @export
calcSenCultureAge <- function(DNAm, pheno = NULL, ID = "Sample_ID") {
  make_simple_clock(
    obj = SenCultureAge_CpGs,
    intercept = -2.546817e+02,
    clock = "SenCultureAge",
    cpg = "CpG",
    coefficients = "Coefficient"
  )(DNAm, pheno, ID)
}

#' SenMortalityAge
#'
#' A function to calculate SenMortalityAge
#'
#' @inheritParams param_template
#' @inherit param_template return
#' @export
calcSenMortalityAge <- function(DNAm, pheno = NULL, ID = "Sample_ID") {
  make_simple_clock(
    obj = SenMortalityAge_CpGs,
    intercept = 0,
    clock = "SenMortalityAge",
    cpg = "CpG",
    coefficients = "Coefficient"
  )(DNAm, pheno, ID)
}
