#' SenMortalityAge
#'
#' @description A function to calculate SenMortalityAge
#'
#' @param DNAm a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param pheno Optional: The sample phenotype data (also with samples as rows) that the clock will be appended to.
#' @param CpGImputation An optional namesd vector for the mean value of each CpG that will be input from another dataset if such values are missing here (from sample cleaning)
#' @param imputation Logical value that will allows you to perform (T)/ skip (F) imputation of mean values for missing CpGs. Warning: when imputation = F if there are missing CpGs, it will automatically ignore these CpGs during calculation, making the clock values less accurate.
#'
#' @return If you added the optional pheno input (preferred) the function appends a column with the clock calculation and returns the dataframe. Otherwise, it will return a vector of calculated clock values in order of the
#' @export
#'
#' @examples calcPhenoAge(exampleBetas, examplePheno, imputation = F)
#'
calcSenMortalityAge <- function(datMeth, datPheno) {
  datMeth <- as.data.frame(datMeth)

  if (length(c(SenMortalityAge_CpGs$CpG[!(SenMortalityAge_CpGs$CpG %in% colnames(datMeth))], SenMortalityAge_CpGs$CpG[apply(datMeth[, colnames(datMeth) %in% SenMortalityAge_CpGs$CpG], 2, function(x) all(is.na(x)))])) == 0) {
    message("CellPopAge - No CpGs were NA for all samples")
  } else {
    missingCpGs_SenMortalityAge <- c(SenMortalityAge_CpGs$CpG[!(SenMortalityAge_CpGs$CpG %in% colnames(datMeth))])
    datMeth[, missingCpGs_SenMortalityAge] <- NA
    datMeth <- datMeth[, SenMortalityAge_CpGs$CpG]
    missingCpGs_SenMortalityAge <- SenMortalityAge_CpGs$CpG[apply(datMeth[, SenMortalityAge_CpGs$CpG], 2, function(x) all(is.na(x)))]
    for (i in 1:length(missingCpGs_SenMortalityAge)) {
      if (!is.na(imputeMissingProbes[missingCpGs_SenMortalityAge[i]])) {
        datMeth[, missingCpGs_SenMortalityAge[i]] <- imputeMissingProbes[missingCpGs_SenMortalityAge[i]]
      } else {
        datMeth[, missingCpGs_SenMortalityAge[i]] <- 0
        # warning("Tried to impute a missing CpG column with Hannum dataset, but failed")
        # warning("Possibly hannum is missing that CpG too")
        # warning("Value replaced with 0")
        # warning("If you have questions about this warning - ask Jess and Dan!! :D ")
      }
    }
    message(paste0("SenMortalityAge needed to fill in ", length(missingCpGs_SenMortalityAge), " CpGs..."))
  }

  DNAmAge <- datPheno


  # Calculate the epigenetic age using the dot product of the coefficients and the methylation data
  predictions <- as.matrix(datMeth[, SenMortalityAge_CpGs$CpG, drop = FALSE]) %*% SenMortalityAge_CpGs$coefficient
  # predictions = predictions + intercept

  datPheno$SenMortalityAge <- predictions

  return(datPheno)
}
