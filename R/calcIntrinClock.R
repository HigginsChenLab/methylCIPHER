#'Intrin Clock
#' @description A function to calculate IntrinClock
#'
#' @param DNAm a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param pheno The sample phenotype data (also with samples as rows) that the clock will be appended to.
#' @return If you added the optional pheno input (preferred) the function appends a column with the clock calculation and returns the dataframe. Otherwise, it will return a vector of calculated clock values in order of the
#' @export
#'
#' @examples calcEpiDish(exampleBetas, examplePheno)
#'
#'

calcIntrinClock<-function(datMeth,datPheno){

  data(IntrinClockCpGs)

  #Create the function for calculating IntrinClock
  returnAge = function(ages) {
    adult_age = 20
    limit = 0
    for (age in 1:length(ages)) {
      if (ages[age] >= limit) {
        ages[age] = (((adult_age+1) * (ages[age]) + adult_age))
      } else {
        ages[age] = (adult_age+1)*exp(ages[age])-1
      }
    }
    return (ages)
  }

  #Loads the glmnet model
  model = readRDS(paste0(path_to_IntrinClock,"final_model_small.RData"))
  allowed_cpgs = rownames(coef(model))[-1]

  intrin_betas = data.frame(datMeth)

  if(length(c(allowed_cpgs[!(allowed_cpgs %in% colnames(intrin_betas))],allowed_cpgs[apply(intrin_betas[,colnames(intrin_betas) %in% allowed_cpgs], 2, function(x)all(is.na(x)))])) == 0){
    message("IntrinClock - No CpGs were NA for all samples")
  } else{

    missingCpGs <- c(allowed_cpgs[!(if(length(c(coefficients_df$name[!(coefficients_df$name %in% colnames(datMeth))],coefficients_df$name[apply(datMeth[,colnames(datMeth) %in% coefficients_df$name], 2, function(x)all(is.na(x)))])) == 0){
      message("IntrinClock - No CpGs were NA for all samples")
    } else{

      missingCpGs <- c(allowed_cpgs[!(allowed_cpgs %in% colnames(intrin_betas))])
      intrin_betas[,missingCpGs]<-NA
      intrin_betas = intrin_betas[,allowed_cpgs]
      missingCpGs <- allowed_cpgs[apply(intrin_betas[,allowed_cpgs], 2, function(x)all(is.na(x)))]

      for(i in 1:length(missingCpGs)){
        if (!is.na(imputeMissingCpGs_IntrinClock[missingCpGs[i]])){
          intrin_betas[,missingCpGs[i]] <- imputeMissingCpGs_IntrinClock[missingCpGs[i]]
          print("CpG found in Hannum")
        }
        else{
          intrin_betas[,missingCpGs[i]] <- 0
          warning("Tried to impute a missing CpG column with Hannum dataset, but failed: Value replaced with 0")
        }
      }
      message(paste0("IntrinClock had ,", length(missingCpGs), " missing CpGs/columns"))
    } %in% colnames(datMeth))])
    datMeth[,missingCpGs]<-NA
    datMeth = datMeth[,coefficients_df$name]
    missingCpGs <- coefficients_df$name[apply(datMeth[,coefficients_df$name], 2, function(x)all(is.na(x)))]

    for(i in 1:length(missingCpGs)){
      if (!is.na(imputeMissingCpGs_RetroClock[missingCpGs[i]])){
        datMeth[,missingCpGs[i]] <- imputeMissingCpGs_RetroClock[missingCpGs[i]]
        #print(imputeMissingCpGs_GrimAge1[missingCpGs_GrimAge1[i]])
        print("CpG found in Hannum")
      }
      else{
        datMeth[,missingCpGs[i]] <- 0
        warning("Tried to impute a missing CpG column with Hannum dataset, but failed: Value replaced with 0")
        #warning("If you see this message - tell Jess and Dan!! :D ")
      }
    }
    message("Any missing CpGs successfully filled in (see function for more details)")
  }


  #Subset the CpG sites of interest
  intrin_betas = intrin_betas[,allowed_cpgs]
  #rownames(intrin_betas) = sub("X","",rownames(intrin_betas))

  #Calculate IntrinClock
  ages = predict(model,as.matrix(intrin_betas))
  for (i in 1:length(ages)){
    if (isNA(ages[i])){
      ages[i]<-0
    }
  }
  intrin_ages = returnAge(ages)
  datPheno$IntrinClock = as.numeric(intrin_ages)

  return(datPheno)

}
