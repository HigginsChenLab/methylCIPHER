#' @title CellDRIFT Model Data
#'
#' @description Bundled PCA models, elastic net model, and module mean values
#'   for calculating CellDRIFT.
#'
#' @format A list with six elements:
#' \describe{
#'   \item{PCA_yellow}{prcomp object for the yellow WGCNA module (2,025 CpGs, 31 PCs)}
#'   \item{PCA_tan}{prcomp object for the tan WGCNA module (297 CpGs, 31 PCs)}
#'   \item{Fit}{Trained glmnet elastic net model}
#'   \item{CV}{cv.glmnet cross-validation object (lambda.min = 0.928)}
#'   \item{yellow_means}{Named numeric vector of training means for yellow module CpGs}
#'   \item{tan_means}{Named numeric vector of training means for tan module CpGs}
#' }
#'
#' @source \url{https://github.com/MorganLevineLab/CellDRIFT}
"CellDRIFT_data"

#' @title CellDRIFT CpGs
#'
#' @description The CpG probes required by CellDRIFT (2,322 CpGs from yellow
#'   and tan WGCNA modules).
#'
#' @format A data frame with 2,322 rows and 1 variable:
#' \describe{
#'   \item{CpG}{CpG probe ID}
#' }
#'
#' @source \url{https://github.com/MorganLevineLab/CellDRIFT}
"CellDRIFT_CpGs"
