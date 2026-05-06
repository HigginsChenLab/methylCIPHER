#' @title Ying Causality-Enriched Epigenetic Clock Weights
#'
#' @description The CpGs and Weights to calculate the Ying causality-enriched
#' epigenetic clocks: `CausAge` (causal aging), `AdaptAge` (adaptive aging),
#' and `DamAge` (damaging aging).
#'
#' @format A data frame with 2 variables:
#' \describe{
#'   \item{CpG}{The CpG names from Illumina Array IDs used in the current clock}
#'   \item{Beta}{The weighted linear regression beta value of fit for each CpG in the clock.}
#' }
#'
#' @references
#' Please cite the source when using these clocks:
#'
#' Ying, K., Liu, H., Tarkhov, A.E., Sadler, M.C., Lu, A.T., Moqri, M., Horvath,
#' S., Kutalik, Z., Shen, X., and Gladyshev, V.N., 2024.
#' Causality-enriched epigenetic age uncouples damage and adaptation.
#' Nature Aging, 4(2), pp.231-246. <https://doi.org/10.1038/s43587-023-00557-0>
#'
#' @name CausalityAge_CpGs
NULL
#' @rdname CausalityAge_CpGs
"CausAge_CpGs"
#' @rdname CausalityAge_CpGs
"AdaptAge_CpGs"
#' @rdname CausalityAge_CpGs
"DamAge_CpGs"

#' Calculate Ying AdaptAge
#'
#' A function to calculate the AdaptAge adaptive-aging epigenetic clock.
#'
#' @inheritParams param_template
#' @inherit param_template return
#' @inherit CausalityAge_CpGs references
#'
#' @export
#'
#' @examples
#' AdaptAge <- calcAdaptAge(
#'   exampleBetas,
#'   examplePheno
#' )
#' AdaptAge
calcAdaptAge <- function(DNAm, pheno = NULL) {
  make_simple_clock(
    obj = AdaptAge_CpGs,
    intercept = -511.9742762,
    clock = "AdaptAge",
    cpg = "CpG",
    coefficients = "Beta"
  )(DNAm, pheno)
}

#' Calculate Ying CausAge
#'
#' A function to calculate the CausAge causality-enriched epigenetic clock.
#'
#' @inheritParams param_template
#' @inherit param_template return
#' @inherit CausalityAge_CpGs references
#'
#' @export
#'
#' @examples
#' CausAge <- calcCausAge(
#'   exampleBetas,
#'   examplePheno
#' )
#' CausAge
calcCausAge <- function(DNAm, pheno = NULL) {
  make_simple_clock(
    obj = CausAge_CpGs,
    intercept = 86.80816381,
    clock = "CausAge",
    cpg = "CpG",
    coefficients = "Beta"
  )(DNAm, pheno)
}

#' Calculate Ying DamAge
#'
#' A function to calculate the DamAge damaging-aging epigenetic clock.
#'
#' @inheritParams param_template
#' @inherit param_template return
#' @inherit CausalityAge_CpGs references
#'
#' @export
#'
#' @examples
#' DamAge <- calcDamAge(
#'   exampleBetas,
#'   examplePheno
#' )
#' DamAge
calcDamAge <- function(DNAm, pheno = NULL) {
  make_simple_clock(
    obj = DamAge_CpGs,
    intercept = 543.4315887,
    clock = "DamAge",
    cpg = "CpG",
    coefficients = "Beta"
  )(DNAm, pheno)
}
