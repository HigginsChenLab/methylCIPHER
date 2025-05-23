#' calcZhang2019
#'
#' @description A function to calculate the Zhang 2019, ultra-precise age measure
#'
#' @param DNAm a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param pheno Optional: The sample phenotype data (also with samples as rows) that the clock will be appended to.
#' @param CpGImputation An optional namesd vector for the mean value of each CpG that will be input from another dataset if such values are missing here (from sample cleaning)
#' @param imputation Logical value that will allows you to perform (T)/ skip (F) imputation of mean values for missing CpGs. Warning: when imputation = F if there are missing CpGs, it will automatically ignore these CpGs during calculation, making the clock values less accurate.
#'
#' @return If you added the optional pheno input (preferred) the function appends a column with the clock calculation and returns the dataframe. Otherwise, it will return a vector of calculated clock values in order of the
#' @export
#'
#' @examples calcZhang2019(exampleBetas, examplePheno, imputation = T)
calcZhang2019 <- function(DNAm, pheno = NULL, CpGImputation = NULL, imputation = T){

  #######################
  ### Read in the Data###
  #######################

  #data("Zhang2019_CpGs")

  ###################################################
  ### Check if all necessary CpGs are in the data ###
  ###################################################

  CpGCheck <- length(Zhang2019_CpGs$CpG) == sum(Zhang2019_CpGs$CpG %in% colnames(DNAm))

  ###################################################################################
  ### The calculation will be performed or an error will be thrown as appropriate ###
  ###################################################################################

  if(CpGCheck == F && is.null(CpGImputation) && imputation == T){

    stop("Need to provide of named vector of CpG Imputations; Necessary CpGs are missing!")

  } else if(CpGCheck == T | imputation == F){

    present <- Zhang2019_CpGs$CpG %in% colnames(DNAm)

    betas <- DNAm[,na.omit(match(Zhang2019_CpGs$CpG,colnames(DNAm)))]
    betas2 <- t(apply(betas,1,scale))
    rownames(betas2) <- rownames(betas)
    colnames(betas2) <- colnames(betas)

    tt <- rowSums(sweep(as.matrix(betas2), MARGIN = 2, Zhang2019_CpGs$coef[present], `*`), na.rm = T) + 65.8

    if(is.null(pheno)){
      tt
    } else{
      pheno$Zhang2019 <- tt
      pheno
    }

  } else {
    message("Imputation of mean CpG Values occured for Zhang2019")
    missingCpGs <- Zhang2019_CpGs$CpG[!(Zhang2019_CpGs$CpG %in% colnames(DNAm))]
    tempDNAm <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs))

    for(j in 1:length(missingCpGs)){
      meanVals <- CpGImputation[match(missingCpGs[j],names(CpGImputation))]
      tempDNAm[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    colnames(tempDNAm) <- missingCpGs
    DNAm <- cbind(DNAm,tempDNAm)

    betas <- DNAm[,match(Zhang2019_CpGs$CpG,colnames(DNAm))]
    betas2 <- t(apply(betas,1,scale))
    rownames(betas2) <- rownames(betas)
    colnames(betas2) <- colnames(betas)

    tt <- rowSums(sweep(betas2, MARGIN = 2, Zhang2019_CpGs$Beta, `*`)) + 65.8

    if(is.null(pheno)){
      tt
    } else{
      pheno$Zhang2019 <- tt
      pheno
    }

  }

}
