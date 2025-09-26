#' calcDamAge
#'
#' @description A function to calculate original 513 CpG DamAge
#'
#' @param DNAm a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param pheno Optional: The sample phenotype data (also with samples as rows) that the clock will be appended to.
#' @param CpGImputation An optional namesd vector for the mean value of each CpG that will be input from another dataset if such values are missing here (from sample cleaning)
#' @param imputation Logical value that will allows you to perform (T)/ skip (F) imputation of mean values for missing CpGs. Warning: when imputation = F if there are missing CpGs, it will automatically ignore these CpGs during calculation, making the clock values less accurate.
#'
#' @return If you added the optional pheno input (preferred) the function appends a column with the clock calculation and returns the dataframe. Otherwise, it will return a vector of calculated clock values in order of the
#' @export
#'
#' @examples calcDamAge(exampleBetas, examplePheno, imputation = F)
calcDamAge <- function(DNAm, pheno = NULL, CpGImputation = NULL, imputation = F) {
  CpGCheck <- length(DamAge_CpGs$CpG) == sum(DamAge_CpGs$CpG %in%
    colnames(DNAm))
  if (CpGCheck == F && is.null(CpGImputation) && imputation ==
    T) {
    stop("Need to provide of named vector of CpG Imputations; Necessary CpGs are missing!")
  } else if (CpGCheck == T | imputation == F) {
    present <- DamAge_CpGs$CpG %in% colnames(DNAm)
    betas <- DNAm[, na.omit(match(DamAge_CpGs$CpG, colnames(DNAm)))]
    tt <- sweep(betas,
      MARGIN = 2, DamAge_CpGs$Beta[present],
      `*`
    )
    DamAge <- as.numeric(rowSums(tt, na.rm = T) + 543.4315887)
    if (is.null(pheno)) {
      DamAge
    } else {
      pheno$DamAge <- DamAge
      pheno
    }
  } else {
    message("Imputation of mean CpG Values occured for DamAge")
    missingCpGs <- DamAge_CpGs$CpG[!(DamAge_CpGs$CpG %in%
      colnames(DNAm))]
    tempDNAm <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs))
    for (j in 1:length(missingCpGs)) {
      meanVals <- CpGImputation[match(
        missingCpGs[j],
        names(CpGImputation)
      )]
      tempDNAm[, j] <- rep(meanVals, dim(DNAm)[1])
    }
    colnames(tempDNAm) <- missingCpGs
    DNAm <- cbind(DNAm, tempDNAm)
    betas <- DNAm[, match(DamAge_CpGs$CpG, colnames(DNAm))]
    tt <- sweep(betas, MARGIN = 2, DamAge_CpGs$Beta, `*`)
    DamAge <- as.numeric(rowSums(tt, na.rm = T) + 543.4315887)
    if (is.null(pheno)) {
      DamAge
    } else {
      pheno$DamAge <- DamAge
      pheno
    }
  }
}
