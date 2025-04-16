#' calcGrimAge2V2
#'
#' @description A function to calculate GrimAgeV2
#'
#' @param DNAm a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param pheno Optional: The sample phenotype data (also with samples as rows) that the clock will be appended to.
#' @return If you added the optional pheno input (preferred) the function appends a column with the clock calculation and returns the dataframe. Otherwise, it will return a vector of calculated clock values in order of the
#' @export
#'
#' @examples calcPhenoAge(exampleBetas, examplePheno)
#'
#'

calcGrimAgeV2 <- function(datMeth, datPheno){
  if(!("Age" %in% variable.names(datPheno))){
    stop("Error: datPheno must have a column named Age")
  }
  if(!("Female" %in% variable.names(datPheno))){
    stop("Error: datPheno must have a column named Female")
  }
  if(sum(startsWith(colnames(datMeth),"cg")) == 0){
    warning("Warning: It looks like you may need to format datMeth using t(datMeth) to get samples as rows!")
  }

  #In datPheno, rows are samples and columns are phenotypic variables.
  #One of the phenotypic variables must be "Age", and another one "Female" (coded as Female = 1, Male = 0; should be a numeric variable as this will be included in PCGrimAge calculation)
  #Also ensure that the order of datMeth sample IDs matches your phenotype data sample IDs, otherwise your data will be scrambled

  #load(file = paste(path_to_GrimAgeV2_directory,"CalcGrimAgeV2.RData", sep = ""))
  #print("Loaded GrimAgeV2 data.")

  #If needed: Fill in missing CpGs needed for calculation of PCs; use mean values from GSE40279 (Hannum 2013; blood)- note that for other tissues you might prefer to use a different one
  # datMeth <- as.data.frame(datMeth)
  # if(length(c(CpGs[!(CpGs %in% colnames(datMeth))],CpGs[apply(datMeth[,colnames(datMeth) %in% CpGs], 2, function(x)all(is.na(x)))])) == 0){
  #   message("No CpGs were NA for all samples")
  # } else{
  #   missingCpGs <- c(CpGs[!(CpGs %in% colnames(datMeth))])
  #   datMeth[,missingCpGs] <- NA
  #   datMeth = datMeth[,CpGs]
  #   missingCpGs <- CpGs[apply(datMeth[,CpGs], 2, function(x)all(is.na(x)))]
  #   for(i in 1:length(missingCpGs)){
  #     datMeth[,missingCpGs[i]] <- imputeMissingCpGs[missingCpGs[i]]
  #   }
  #   message("Any missing CpGs successfully filled in (see function for more details)")
  # }

  datMeth <- as.data.frame(datMeth)
  if(length(c(CpGs_GrimAge2[!(CpGs_GrimAge2 %in% colnames(datMeth))],CpGs_GrimAge2[apply(datMeth[,colnames(datMeth) %in% CpGs_GrimAge2], 2, function(x)all(is.na(x)))])) == 0){
    message("GrimAgeV2 - No CpGs were NA for all samples")
  } else{
    missingCpGs_GrimAge2 <- c(CpGs_GrimAge2[!(CpGs_GrimAge2 %in% colnames(datMeth))])
    datMeth[,missingCpGs_GrimAge2] <- NA
    datMeth = datMeth[,CpGs_GrimAge2]
    missingCpGs_GrimAge2 <- CpGs_GrimAge2[apply(datMeth[,CpGs_GrimAge2], 2, function(x)all(is.na(x)))]
    for(i in 1:length(missingCpGs_GrimAge2)){
      if (!is.na(imputeMissingCpGs_GrimAge2[missingCpGs_GrimAge2[i]])){
        datMeth[,missingCpGs_GrimAge2[i]] <- imputeMissingCpGs_GrimAge2[missingCpGs_GrimAge2[i]]
        #print(imputeMissingCpGs_GrimAge2[missingCpGs_GrimAge2[i]])
      }
      else{
        datMeth[,missingCpGs_GrimAge2[i]] <- 0
      }
    }
    message(paste0("GrimAgeV2 needed to fill in ", length(missingCpGs_GrimAge2), " CpGs..."))
  }

  #Prepare methylation data for calculation of PC Clocks (subset to 78,464 CpGs and perform imputation if needed)
  #datMeth <- datMeth[,CpGs]
  #meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
  #datMeth <- apply(datMeth,2,meanimpute)
  #Note: you may substitute another imputation method of your choice (e.g. KNN), but we have not found the method makes a significant difference.
  #message("Mean imputation successfully completed for any missing CpG values")

  #Initialize a data frame for PC clocks
  DNAmAge <- datPheno

  # Calculate GrimAge Clocks
  temp <- as.matrix(cbind(datMeth, Age = as.numeric(datPheno$Age, Female=datPheno$Female)))
  DNAmAge$DNAmPACKYRS <- as.numeric(temp[,names(CalcGrimAge2$PACKYRS.model)] %*% CalcGrimAge2$PACKYRS.model + CalcGrimAge2$PACKYRS.intercept)
  DNAmAge$DNAmADM <- as.numeric(temp[,names(CalcGrimAge2$ADM.model)] %*% CalcGrimAge2$ADM.model + CalcGrimAge2$ADM.intercept)
  DNAmAge$DNAmB2M <- as.numeric(temp[,names(CalcGrimAge2$B2M.model)] %*% CalcGrimAge2$B2M.model + CalcGrimAge2$B2M.intercept)
  DNAmAge$DNAmCystatinC <- as.numeric(temp[,names(CalcGrimAge2$CystatinC.model)] %*% CalcGrimAge2$CystatinC.model + CalcGrimAge2$CystatinC.intercept)
  DNAmAge$DNAmGDF15 <- as.numeric(temp[,names(CalcGrimAge2$GDF15.model)] %*% CalcGrimAge2$GDF15.model + CalcGrimAge2$GDF15.intercept)
  DNAmAge$DNAmLeptin <- as.numeric(temp[,names(CalcGrimAge2$Leptin.model)] %*% CalcGrimAge2$Leptin.model + CalcGrimAge2$Leptin.intercept)
  DNAmAge$DNAmPAI1 <- as.numeric(temp[,names(CalcGrimAge2$PAI1.model)] %*% CalcGrimAge2$PAI1.model + CalcGrimAge2$PAI1.intercept)
  DNAmAge$DNAmTIMP1 <- as.numeric(temp[,names(CalcGrimAge2$TIMP1.model)] %*% CalcGrimAge2$TIMP1.model + CalcGrimAge2$TIMP1.intercept)
  #New Proteins
  DNAmAge$DNAmlogA1C <- as.numeric(temp[,names(CalcGrimAge2$logA1C.model)] %*% CalcGrimAge2$logA1C.model + CalcGrimAge2$logA1C.intercept)
  DNAmAge$DNAmlogCRP <- as.numeric(temp[,names(CalcGrimAge2$logCRP.model)] %*% CalcGrimAge2$logCRP.model + CalcGrimAge2$logCRP.intercept)

  # DNAmAge$Age<-DNAmAge$cAGE
  # DNAmAge$Female<-DNAmAge$cFEMALE

  DNAmAge$GrimAgeV2 <- as.numeric(as.matrix(DNAmAge[,CalcGrimAge2$components]) %*% CalcGrimAge2$GrimAge.model ) #+ CalcGrimAge2$GrimAge.intercept)
  y<-DNAmAge$GrimAgeV2
  #print(y)
  DNAmAge$GrimAgeV2 <- (((y- CalcGrimAge2$GrimAge.transform[3])/CalcGrimAge2$GrimAge.transform[4])*CalcGrimAge2$GrimAge.transform[2]) + CalcGrimAge2$GrimAge.transform[1]

  return(DNAmAge)

}
