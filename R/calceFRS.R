#' calceFRS
#'
#' @description A function to calculate the 20 CpG eFRS
#'
#' @param DNAm a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param pheno Optional: The sample phenotype data (also with samples as rows) that the clock will be appended to.
#' @param CpGImputation An optional namesd vector for the mean value of each CpG that will be input from another dataset if such values are missing here (from sample cleaning)
#' @param imputation Logical value that will allows you to perform (T)/ skip (F) imputation of mean values for missing CpGs. Warning: when imputation = F if there are missing CpGs, it will automatically ignore these CpGs during calculation, making the clock values less accurate.
#'
#' @return If you added the optional pheno input (preferred) the function appends a column with the clock calculation and returns the dataframe. Otherwise, it will return a vector of calculated clock values in order of the
#' @export
#'
#' @examples calceFRS(exampleBetas, examplePheno, imputation = F)

calceFRS <- function (DNAm, pheno = NULL, CpGImputation = NULL, imputation = F)
{
  CpGCheck <- length(eFRS_CpGs$CpG) == sum(eFRS_CpGs$CpG %in%
                                                 colnames(DNAm))
  if (CpGCheck == F && is.null(CpGImputation) && imputation ==
      T) {
    stop("Need to provide of named vector of CpG Imputations; Necessary CpGs are missing!")
  } else if (CpGCheck == T | imputation == F) {
    present <- eFRS_CpGs$CpG %in% colnames(DNAm)
    betas <- DNAm[, na.omit(match(eFRS_CpGs$CpG, colnames(DNAm)))]
    tt <- sweep(betas, MARGIN = 2, eFRS_CpGs$Beta[present],
                `*`)
    eFRS <- as.numeric(rowSums(tt, na.rm = T) + 0.204)
    if (is.null(pheno)) {
      eFRS
    }
    else {
      pheno$eFRS <- eFRS
      pheno
    }
  } else {
    message("Imputation of mean CpG Values occured for eFRS")
    missingCpGs <- eFRS_CpGs$CpG[!(eFRS_CpGs$CpG %in%
                                         colnames(DNAm))]
    tempDNAm <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs))
    for (j in 1:length(missingCpGs)) {
      meanVals <- CpGImputation[match(missingCpGs[j],
                                      names(CpGImputation))]
      tempDNAm[, j] <- rep(meanVals, dim(DNAm)[1])
    }
    colnames(tempDNAm) <- missingCpGs
    DNAm <- cbind(DNAm, tempDNAm)
    betas <- DNAm[, match(eFRS_CpGs$CpG, colnames(DNAm))]
    tt <- sweep(betas, MARGIN = 2, eFRS_CpGs$Beta, `*`)
    eFRS <- as.numeric(rowSums(tt, na.rm = T) + 0.204)
    if (is.null(pheno)) {
      eFRS
    }
    else {
      pheno$eFRS <- eFRS
      pheno
    }
  }
}

