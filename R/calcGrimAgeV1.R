#' CalcGrimAgeV1
#'
#' @description A function to calculate GrimAgeV1
#'
#' @param DNAm a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param pheno Optional: The sample phenotype data (also with samples as rows) that the clock will be appended to.
#' @param CpGImputation An optional namesd vector for the mean value of each CpG that will be input from another dataset if such values are missing here (from sample cleaning)
#' @param imputation Logical value that will allows you to perform (T)/ skip (F) imputation of mean values for missing CpGs. Warning: when imputation = F if there are missing CpGs, it will automatically ignore these CpGs during calculation, making the clock values less accurate.
#'
#' @return If you added the optional pheno input (preferred) the function appends a column with the clock calculation and returns the dataframe. Otherwise, it will return a vector of calculated clock values in order of the
#' @export
#'
#' @examples CalcGrimAgeV1(exampleBetas, examplePheno, imputation = F)
#'
#'
#'
#'
#'


calcGrimAgeV1 <- function(datMeth, datPheno){

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

  #load(file = paste(path_to_GrimAgeV1_directory,"CalcGrimAgeV1.RData", sep = ""))
  #print("Loaded GrimAgeV1 data.")

  #If needed: Fill in missing CpGs needed for calculation of PCs; use mean values from GSE40279 (Hannum 2013; blood)- note that for other tissues you might prefer to use a different one
  datMeth <- as.data.frame(datMeth)
  if(length(c(CpGs_GrimAge1[!(CpGs_GrimAge1 %in% colnames(datMeth))],CpGs_GrimAge1[apply(datMeth[,colnames(datMeth) %in% CpGs_GrimAge1], 2, function(x)all(is.na(x)))])) == 0){
    message("GrimAgeV1 - No CpGs were NA for all samples")
  } else{
    missingCpGs_GrimAge1 <- c(CpGs_GrimAge1[!(CpGs_GrimAge1 %in% colnames(datMeth))])
    datMeth[,missingCpGs_GrimAge1] <- NA
    datMeth = datMeth[,CpGs_GrimAge1]
    missingCpGs_GrimAge1 <- CpGs_GrimAge1[apply(datMeth[,CpGs_GrimAge1], 2, function(x)all(is.na(x)))]
    for(i in 1:length(missingCpGs_GrimAge1)){
      if (!is.na(imputeMissingCpGs_GrimAge1[missingCpGs_GrimAge1[i]])){
        datMeth[,missingCpGs_GrimAge1[i]] <- imputeMissingCpGs_GrimAge1[missingCpGs_GrimAge1[i]]
        #print(imputeMissingCpGs_GrimAge1[missingCpGs_GrimAge1[i]])
      }
      else{
        datMeth[,missingCpGs_GrimAge1[i]] <- 0

      }
    }
    message(paste0("GrimAgeV1 needed to fill in ", length(missingCpGs_GrimAge1), " CpGs..."))
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
  DNAmAge$DNAmPACKYRS <- as.numeric(temp[,names(CalcGrimAgeV1$PACKYRS.model)] %*% CalcGrimAgeV1$PACKYRS.model + CalcGrimAgeV1$PACKYRS.intercept)
  DNAmAge$DNAmADM <- as.numeric(temp[,names(CalcGrimAgeV1$ADM.model)] %*% CalcGrimAgeV1$ADM.model + CalcGrimAgeV1$ADM.intercept)
  DNAmAge$DNAmB2M <- as.numeric(temp[,names(CalcGrimAgeV1$B2M.model)] %*% CalcGrimAgeV1$B2M.model + CalcGrimAgeV1$B2M.intercept)
  DNAmAge$DNAmCystatinC <- as.numeric(temp[,names(CalcGrimAgeV1$CystatinC.model)] %*% CalcGrimAgeV1$CystatinC.model + CalcGrimAgeV1$CystatinC.intercept)
  DNAmAge$DNAmGDF15 <- as.numeric(temp[,names(CalcGrimAgeV1$GDF15.model)] %*% CalcGrimAgeV1$GDF15.model + CalcGrimAgeV1$GDF15.intercept)
  DNAmAge$DNAmLeptin <- as.numeric(temp[,names(CalcGrimAgeV1$Leptin.model)] %*% CalcGrimAgeV1$Leptin.model + CalcGrimAgeV1$Leptin.intercept)
  DNAmAge$DNAmPAI1 <- as.numeric(temp[,names(CalcGrimAgeV1$PAI1.model)] %*% CalcGrimAgeV1$PAI1.model + CalcGrimAgeV1$PAI1.intercept)
  DNAmAge$DNAmTIMP1 <- as.numeric(temp[,names(CalcGrimAgeV1$TIMP1.model)] %*% CalcGrimAgeV1$TIMP1.model + CalcGrimAgeV1$TIMP1.intercept)
  #print("Components Calculated.Calculating GrimAGE...")
  DNAmAge$Age<-DNAmAge$Age
  DNAmAge$Female<-DNAmAge$Female
  DNAmAge$GrimAgeV1 <- as.numeric(as.matrix(DNAmAge[,CalcGrimAgeV1$components]) %*% CalcGrimAgeV1$GrimAge.model ) #+ CalcGrimAgeV1$GrimAge.intercept)
  y<-DNAmAge$GrimAgeV1
  DNAmAge$GrimAgeV1 <- (((y- CalcGrimAgeV1$GrimAge.transform[3])/CalcGrimAgeV1$GrimAge.transform[4])*CalcGrimAgeV1$GrimAge.transform[2]) + CalcGrimAgeV1$GrimAge.transform[1]

  return(DNAmAge)

}
