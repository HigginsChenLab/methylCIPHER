#' calcSystemsAge
#'
#' @description A function to calculate Systems Age
#'
#' @inheritParams param_template
#' @param RData
#' Either a character string specifying the path to the methylCIPHER folder containing
#' the `SystemsAge_data.RData` file, or a list/environment containing the contents of
#' the `SystemsAge_data.RData` file.
#'
#' @return If you added the optional pheno input (preferred) the function appends a column with the clock calculation and returns the data.frame. Otherwise, it will return a data.frame of calculated clock values.
#' @export
#'
#' @examples
#' \dontrun{
#' example_SystemsAge <- calcSystemsAge(
#'   exampleBetas,
#'   examplePheno,
#'   # If `SystemsAge_data.RData` has been downloaded to default folder at $HOME
#'   RData = get_methylCIPHER_path()
#' )
#' example_SystemsAge
#' }
#'
calcSystemsAge <- function(DNAm, pheno = NULL, RData = NULL) {
  if (!is.matrix(DNAm)) {
    # Temporary method before converting all other input into matrix
    DNAm <- as.matrix(DNAm)
  }
  # Input validation
  checkmate::assert(
    checkmate::check_null(RData),
    checkmate::check_character(RData, len = 1, any.missing = FALSE),
    checkmate::check_list(RData, any.missing = FALSE),
    checkmate::check_environment(RData),
    combine = "or"
  )
  checkmate::assert_data_frame(pheno, min.rows = 1, null.ok = TRUE)
  if (!is.null(pheno)) {
    stopifnot("rows in `pheno` must correpsond to the same row as `DNAm`" = nrow(pheno) == nrow(DNAm))
    warning("Currently, checks for rows in `pheno` corresponding to the same row as `DNAm` is not implemented")
  }
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
    if (!"SystemsAge_data.RData" %in% list.files(RData)) {
      stop("The file 'SystemsAge_data.RData' is not found in the provided path.")
    }
    cat("Loading Data...\n")
    SystemsAge_env <- new.env()
    load(file = paste0(RData, "/", "SystemsAge_data.RData"), envir = SystemsAge_env)
  } else if (is.list(RData) || is.environment(RData)) {
    SystemsAge_env <- RData
  }
  attach(SystemsAge_env)
  ## Imputation
  imputed_DNAm <- impute_DNAm(
    DNAm = DNAm,
    method = "mean",
    CpGs = imputeMissingCpGs[intersect(rownames(DNAmPCA$rotation), names(imputeMissingCpGs))],
    subset = TRUE
  )

  ## Calculate methylation PCs
  DNAmPCs <- predict(DNAmPCA, imputed_DNAm)

  ## Calculate DNAm system PCs then system scores
  DNAmSystemPCs <- DNAmPCs[, 1:4017] %*% as.matrix(system_vector_coefficients[1:4017, ])
  system_scores <- matrix(nrow = dim(DNAmSystemPCs)[1], ncol = 11)
  i <- 1
  groups <- c("Blood", "Brain", "Cytokine", "Heart", "Hormone", "Immune", "Kidney", "Liver", "Metab", "Lung", "MusculoSkeletal")
  for (group in groups) {
    tf <- grepl(group, colnames(DNAmSystemPCs))
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
  detach(SystemsAge_env)

  if (is.null(pheno)) {
    return(as.data.frame(system_ages))
  } else {
    return(cbind(pheno, system_ages))
  }
}
