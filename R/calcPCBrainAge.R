#' Calculate PCBrainAge
#'
#' @description Calculates PCBrainAge, a principal component-based epigenetic
#'   clock trained on brain cortical tissue.
#'
#' @inheritParams param_template
#' @param RData
#'   Default to `NULL`, which loads from the default path.
#'   Either a character string specifying the path to the
#'   `PCBrainAge_data.qs2` file, or a list containing the pre-loaded data
#'   from [load_PCBrainAge_data()]. See Details.
#'
#' @details
#' PCBrainAge requires the `PCBrainAge_data.qs2` object, which can be
#' downloaded using [download_methylCIPHER()]. Provide either the path to the
#' downloaded file or the object loaded with [load_PCBrainAge_data()] to the
#' `RData` parameter.
#'
#' @inherit param_template return
#'
#' @references
#' Thrush, K.L., et al. (2022). Aging the brain: multi-region methylation
#' principal component based clock in the context of Alzheimer's disease.
#' Aging, 14(14), 5641-5668.
#' \url{https://doi.org/10.18632/aging.204196}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' download_methylCIPHER(clocks = "PCBrainAge")
#' calcPCBrainAge(exampleBetas, examplePheno)
#' }
calcPCBrainAge <- function(DNAm, pheno = NULL, RData = NULL) {
  # Input validation
  check_DNAm(DNAm)
  checkmate::assert(
    checkmate::check_null(RData),
    checkmate::check_character(RData, len = 1, any.missing = FALSE),
    checkmate::check_list(RData, any.missing = FALSE),
    combine = "or"
  )

  # Load data
  if (is.null(RData)) {
    RData <- load_PCBrainAge_data()
  } else if (is.character(RData)) {
    RData <- load_PCBrainAge_data(RData)
  }

  # Imputation using brain tissue training medians
  DNAm <- impute_DNAm(
    DNAm = DNAm,
    method = "mean",
    CpGs = RData$imputeMissingCpGs,
    subset = TRUE
  )

  # Re-align to match model CpG order
  model_cpgs <- rownames(RData$rotation)
  DNAm <- DNAm[, model_cpgs, drop = FALSE]

  # Calculate: (DNAm - center) %*% rotation %*% coefs + intercept
  predictions <- as.numeric(
    sweep(as.matrix(DNAm), 2, RData$center, "-") %*%
      RData$rotation %*% RData$coefs + RData$intercept
  )

  if (is.null(pheno)) {
    result <- data.frame(PCBrainAge = predictions)
    row.names(result) <- row.names(DNAm)
    return(result)
  }
  pheno$PCBrainAge <- predictions
  return(pheno)
}
