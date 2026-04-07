#' Calculate CellDRIFT
#'
#' @description Calculates CellDRIFT, a DNA methylation-based predictor of
#'   cumulative population doublings using PCA projections from yellow and tan
#'   WGCNA modules combined with elastic net regression.
#'
#' @inheritParams param_template
#'
#' @inherit param_template return
#'
#' @references
#' Minteer, C.J., et al. (2022). More than peeling back the layers:
#' the biological significance of DNA methylation entropy as a measure
#' of cellular proliferation and aging.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' calcCellDRIFT(exampleBetas, examplePheno)
#' }
calcCellDRIFT <- function(DNAm, pheno = NULL) {
  check_DNAm(DNAm)

  models <- CellDRIFT_data
  DNAm <- as.matrix(DNAm)

  # Replace NAs with 0 (matching original behavior)
  if (any(is.na(DNAm))) {
    message("CellDRIFT: NA values set to 0.")
    DNAm[is.na(DNAm)] <- 0
  }

  # --- Yellow module PCA ---
  yellow_cpgs <- names(models$yellow_means)
  missing_y <- setdiff(yellow_cpgs, colnames(DNAm))
  if (length(missing_y) > 0) {
    impute_y <- matrix(
      models$yellow_means[missing_y],
      nrow = nrow(DNAm), ncol = length(missing_y), byrow = TRUE,
      dimnames = list(row.names(DNAm), missing_y)
    )
    DNAm <- cbind(DNAm, impute_y)
  }
  Betas_y <- DNAm[, yellow_cpgs, drop = FALSE]
  PCs_yellow <- predict(models$PCA_yellow, Betas_y)
  colnames(PCs_yellow) <- paste0("PC", 1:31, "_y")

  # --- Tan module PCA ---
  tan_cpgs <- names(models$tan_means)
  missing_t <- setdiff(tan_cpgs, colnames(DNAm))
  if (length(missing_t) > 0) {
    impute_t <- matrix(
      models$tan_means[missing_t],
      nrow = nrow(DNAm), ncol = length(missing_t), byrow = TRUE,
      dimnames = list(row.names(DNAm), missing_t)
    )
    DNAm <- cbind(DNAm, impute_t)
  }
  Betas_t <- DNAm[, tan_cpgs, drop = FALSE]
  PCs_tan <- predict(models$PCA_tan, Betas_t)
  colnames(PCs_tan) <- paste0("PC", 1:31, "_t")

  # --- Combine and predict ---
  PCs_combined <- cbind(PCs_yellow, PCs_tan)
  predictions <- as.numeric(
    glmnet::predict.glmnet(models$Fit, PCs_combined, s = models$CV$lambda.min)
  )

  if (is.null(pheno)) {
    result <- data.frame(CellDRIFT = predictions)
    row.names(result) <- row.names(DNAm)
    return(result)
  }
  pheno$CellDRIFT <- predictions
  return(pheno)
}
