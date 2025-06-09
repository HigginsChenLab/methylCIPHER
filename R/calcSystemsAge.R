#' calcSystemsAge
#'
#' @description A function to calculate Systems Age
#'
#' @param DNAm a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param pheno Optional: The sample phenotype data (also with samples as rows) that the clock will be appended to.
#' @param imputation Logical value that will allows you to perform (T)/ skip (F) imputation of mean values for missing CpGs. Warning: when imputation = F if there are missing CpGs, it will automatically ignore these CpGs during calculation, making the clock values less accurate.
#' @param RData
#' Either a character string specifying the path to the methylCIPHER folder containing
#' the `SystemsAge_data.RData` file, or a list/environment containing the contents of
#' the `SystemsAge_data.RData` file.
#'
#' @return If you added the optional pheno input (preferred) the function appends a column with the clock calculation and returns the dataframe. Otherwise, it will return a vector of calculated clock values in order of the
#' @export
#'
#' @examples
#' \dontrun{
#' example_SystemsAge <- calcSystemsAge(
#'   exampleBetas,
#'   examplePheno,
#'   imputation = F,
#'   # If `SystemsAge_data.RData` has been downloaded to default folder at $HOME
#'   RData = get_methylCIPHER_path()
#'  )
#' example_SystemsAge
#' }
#'
calcSystemsAge <- function(DNAm, pheno = NULL, imputation = F, RData = NULL) {
  if(!is.matrix(DNAm)) {
    # Temporary method before converting all other input into matrix
    DNAm <- as.matrix(DNAm)
  }
  # Input validation
  checkmate::assert_logical(imputation, len = 1, any.missing = F)
  checkmate::assert(
    checkmate::check_null(RData),
    checkmate::check_character(RData, len = 1, any.missing = FALSE),
    checkmate::check_list(RData, any.missing = FALSE),
    checkmate::check_environment(RData),
    combine = "or"
  )
  check_DNAm(DNAm)

  if (is.null(RData)) {
    if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path)) {
      current_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
    } else {
      current_path <- getwd() # fallback to current working directory
    }

    setwd(current_path)
    # add more print messages to see progress in computation

    # Code for calculating Systems Age
    # Main Authors:
    # Raghav Sehgal, Yale University, raghav.sehgal@yale.edu
    # Albert T. Higgins-Chen, MD, PhD, Yale University, a.higginschen@yale.edu
    # Morgan E. Levine, PhD, Altos Labs, morgan.levine@yale.edu

    # Code last updated May 17, 2023

    # Load packages and data
    # Change the file directory as needed
    pkgTest <- function(x) {
      if (!require(x, character.only = TRUE)) {
        install.packages(x, dep = TRUE)
        if (!require(x, character.only = TRUE)) stop("Package not found")
      }
    }
    pkgTest("stringr")
    library(stringr)
    pkgTest("googledrive")
    library(googledrive)

    googledrive::drive_auth()

    # In pheno, rows are samples and columns are phenotypic variables.
    # One of the phenotypic variables must be "Age", and another one "Female" (coded as Female = 1, Male = 0; should be a numeric variable as this will be included in PCGrimAge calculation)
    # Also ensure that the order of DNAm sample IDs matches your phenotype data sample IDs, otherwise your data will be scrambled

    home_dir <- Sys.getenv("HOME")

    if (file.exists(paste0(home_dir, "/SystemsAge_data.RData"))) {
      cat("SystemsAge Data File Exists in Home Directory. Loading Data...\n")
      load(paste0(home_dir, "/SystemsAge_data.RData"))
    } else {
      cat("SystemsAge Data does not exist in package directory. Downloading data...\n")
      public_file <- drive_get(as_id("14HBIhA-gkp9-RWD0HVAekRGWOHQBqvhc"))
      drive_download(public_file, path = paste0(home_dir, "/SystemsAge_data.RData"), overwrite = TRUE)
      load(paste0(home_dir, "/SystemsAge_data.RData"))
    }
  } else if (is.character(RData)) {
    if(!"SystemsAge_data.RData" %in% list.files(RData)) {
      stop("The file 'SystemsAge_data.RData' is not found in the provided path.")
    }
    cat("Loading Data...\n")
    SystemsAge_env <- new.env()
    load(file = paste0(RData, "/", "SystemsAge_data.RData"), envir = SystemsAge_env)
    attach(SystemsAge_env)
  } else if (is.list(RData) || is.environment(RData)) {
    SystemsAge_env <- RData
    attach(SystemsAge_env)
  }
  stop("stop here")
  ## Access CpG names from PCA rotation matrix
  CpGs <- rownames(DNAmPCA$rotation)

  # Examine missing values
  # Two types of missing values: CpGs are missing entirely, or CpGs are missing for some samples
  # For CpGs missing entirely:
  CpGs_present <- colnames(DNAm) %in% CpGs
  CpGs_missing <- CpGs[!(CpGs %in% colnames(DNAm))]
  if (sum(CpGs_present) == length(CpGs)) {
    message("Systems Age says all CpGs present, proceed to next step")
  } else {
    message(paste0("Systems Age says missing ", length(CpGs_missing), " CpGs, please check variable CpGs_missing for list of missing CpGs"))
  }
  length(CpGs_missing)
  # For CpGs missing for some samples:
  sum(is.na(DNAm))
  # Please note how many CpGs of each type are missing and report it when publishing.
  # If there are too many CpGs missing, e.g. >10%, consider re-processing your data using a different pipeline.

  # Perform mean imputation for CpGs missing for some samples:
  # Note that you may select a different imputation method if you choose
  # meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
  # DNAm <- apply(DNAm,2,meanimpute)

  # For CpGs are missing entirely, re-add them by using mean values from GSE40279 (Hannum 2013; blood)
  # Be wary of applying this imputation to non-blood tissues (note that Systems Age is trained to quantify aging using blood DNAm)
  DNAm <- data.frame(DNAm)
  DNAm[, CpGs_missing] <- NA
  if (length(CpGs_missing) > 0) {
    for (i in 1:length(CpGs_missing)) {
      DNAm[, CpGs_missing[i]] <- imputeMissingCpGs[CpGs_missing[i]]
    }
  }

  ## Subset Methylation matrix to only those CpGs used for calculation of Systems Age
  meth_df <- DNAm[, CpGs]

  ## Calculate methylation PCs
  DNAmPCs <- predict(DNAmPCA, meth_df)

  ## Calculate DNAm system PCs then system scores
  DNAmSystemPCs <- DNAmPCs[, 1:4017] %*% as.matrix(system_vector_coefficients[1:4017, ])
  system_scores <- matrix(nrow = dim(DNAmSystemPCs)[1], ncol = 11)
  i <- 1
  groups <- c("Blood", "Brain", "Cytokine", "Heart", "Hormone", "Immune", "Kidney", "Liver", "Metab", "Lung", "MusculoSkeletal")
  for (group in groups) {
    tf <- str_detect(colnames(DNAmSystemPCs), group)
    sub <- DNAmSystemPCs[, tf]
    sub_system_coefficients <- system_scores_coefficients_scale[tf]
    if (length(sub_system_coefficients) == 1) {
      system_scores[, i] <- sub * -1
    } else {
      system_scores[, i] <- sub %*% sub_system_coefficients
    }
    i <- i + 1
  }
  groups <- c("Blood", "Brain", "Cytokine", "Heart", "Hormone", "Immune", "Kidney", "Liver", "Metab", "Lung", "MusculoSkeletal")
  colnames(system_scores) <- groups

  ## Generate predicted chronological age
  Age_prediction <- as.matrix(DNAmPCs) %*% as.matrix(Predicted_age_coefficients[2:4019]) + Predicted_age_coefficients[1]
  Age_prediction <- Age_prediction * Age_prediction_model[2] + (Age_prediction**2) * Age_prediction_model[3] + Age_prediction_model[1]
  Age_prediction <- Age_prediction / 12
  system_scores <- cbind(system_scores, Age_prediction)
  colnames(system_scores)[12] <- "Age_prediction"

  ## Generating overall system index
  colnames(system_scores) <- c("Blood", "Brain", "Inflammation", "Heart", "Hormone", "Immune", "Kidney", "Liver", "Metabolic", "Lung", "MusculoSkeletal", "Age_prediction")
  system_PCA <- predict(systems_PCA, system_scores)
  pred <- system_PCA %*% Systems_clock_coefficients
  system_scores <- cbind(system_scores, pred)
  colnames(system_scores)[13] <- "SystemsAge"

  ## Scale system ages
  system_ages <- system_scores
  for (i in c(1:13)) {
    y <- system_ages[, i]
    system_ages[, i] <- (((y - transformation_coefs[i, 1]) / transformation_coefs[i, 2]) * transformation_coefs[i, 4]) + transformation_coefs[i, 3]
    system_ages[, i] <- system_ages[, i] / 12
  }

  # ## Append to your data frame (pheno)
  # ## Check first that data is in correct order
  # #all(rownames(pheno$Sample_ID) == rownames(system_ages))
  # DNAmAge <- data.frame(cbind(pheno,system_ages))

  detach(SystemsAge_env)
  if (is.null(pheno)) {
    print(system_ages)
  } else {
    all(rownames(pheno$Sample_ID) == rownames(system_ages))
    DNAmAge <- data.frame(cbind(pheno, system_ages))
    return(DNAmAge)
  }
}
