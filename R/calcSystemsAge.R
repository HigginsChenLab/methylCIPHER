#' calcSystemsAge
#'
#' @description A function to calculate Systems Age
#'
#' @inheritParams param_template
#' @param RData
#' Either a character string specifying the path to the methylCIPHER folder containing
#' the `SystemsAge_data.qs2` file or a list containing the contents of
#' the `SystemsAge_data.qs2` file created by [load_SystemsAge_data()].
#'
#' @return If you added the optional pheno input (preferred) the function appends a column with the clock calculation and returns the data.frame. Otherwise, it will return a data.frame of calculated clock values.
#' @export
#'
#' @examples
#' \dontrun{
#' example_SystemsAge <- calcSystemsAge(
#'   exampleBetas,
#'   examplePheno,
#'   RData = get_methylCIPHER_path()
#' )
#' example_SystemsAge
#' }
#'
calcSystemsAge <- function(DNAm, pheno = NULL, RData = NULL) {
  # Code for calculating Systems Age
  # Main Authors:
  # Raghav Sehgal, Yale University, raghav.sehgal@yale.edu
  # Albert T. Higgins-Chen, MD, PhD, Yale University, a.higginschen@yale.edu
  # Morgan E. Levine, PhD, Altos Labs, morgan.levine@yale.edu
  # Code last updated May 17, 2023

  # Input validation
  # Check DNAm
  check_DNAm(DNAm)
  # Check RData
  checkmate::assert(
    checkmate::check_null(RData),
    checkmate::check_character(RData, len = 1, any.missing = FALSE),
    checkmate::check_list(RData, any.missing = FALSE),
    combine = "or"
  )
  # Check Pheno
  if (!is.null(pheno)) {
    check_pheno(pheno)
    # Check Consistent between `pheno` and `DNAm`
    need_align <- if (length(row.names(DNAm)) != length(pheno$Sample_ID)) {
      TRUE
    } else if (any(row.names(DNAm) != pheno$Sample_ID)) {
      TRUE
    } else {
      FALSE
    }
    if (need_align) {
      samples <- intersect(row.names(DNAm), pheno$Sample_ID)
      DNAm <- DNAm[samples, , drop = FALSE]
      pheno <- align_pheno(pheno, samples)
      stopifnot("`DNAm` and `pheno` samples alignment failed. Check `Sample_ID` of pheno and row.names() of `DNAm`" = all.equal(row.names(DNAm), pheno$Sample_ID))
      message("Samples inconsistencies between DNAm and Pheno were detected and corrected.")
    }
  }
  # handle RData
  if (is.null(RData)) {
    RData <- load_SystemsAge_data()
  }
  if (is.character(RData)) {
    info <- file.info(RData)
    # if pass a folder
    if (info$isdir) {
      # then convert path to read into folder/file.qs2
      RData <- file.path(RData, "SystemsAge_data.qs2")
    }
    # then convert RData into the list of objects to be used
    RData <- load_SystemsAge_data(RData)
  }

  ## Imputation
  DNAm <- impute_DNAm(
    DNAm = DNAm,
    method = "mean",
    CpGs = RData$imputeMissingCpGs,
    subset = TRUE
  )

  ## Re-align to make sure things lined up with the object
  DNAm <- DNAm[, names(RData$imputeMissingCpGs), drop = F]
  stopifnot(all(colnames(DNAm) == names(RData$imputeMissingCpGs)))

  ## Calculate methylation PCs
  DNAmPCs <- predict(RData$DNAmPCA, DNAm)

  ## Calculate DNAm system PCs then system scores
  DNAmSystemPCs <- DNAmPCs[, 1:4017] %*% as.matrix(RData$system_vector_coefficients[1:4017, ])
  system_scores <- matrix(nrow = dim(DNAmSystemPCs)[1], ncol = 11)
  i <- 1
  groups <- c("Blood", "Brain", "Cytokine", "Heart", "Hormone", "Immune", "Kidney", "Liver", "Metab", "Lung", "MusculoSkeletal")
  for (group in groups) {
    tf <- grepl(group, colnames(DNAmSystemPCs))
    sub <- DNAmSystemPCs[, tf]
    sub_system_coefficients <- RData$system_scores_coefficients_scale[tf]
    if (length(sub_system_coefficients) == 1) {
      system_scores[, i] <- sub * -1
    } else {
      system_scores[, i] <- sub %*% sub_system_coefficients
    }
    i <- i + 1
  }
  colnames(system_scores) <- groups

  ## Generate predicted chronological age
  Age_prediction <- as.matrix(DNAmPCs) %*% as.matrix(RData$Predicted_age_coefficients[2:4019]) + RData$Predicted_age_coefficients[1]
  Age_prediction <- Age_prediction * RData$Age_prediction_model[2] + (Age_prediction**2) * RData$Age_prediction_model[3] + RData$Age_prediction_model[1]
  Age_prediction <- Age_prediction / 12
  system_scores <- cbind(system_scores, Age_prediction)
  colnames(system_scores)[12] <- "Age_prediction"

  ## Generating overall system index
  colnames(system_scores) <- c("Blood", "Brain", "Inflammation", "Heart", "Hormone", "Immune", "Kidney", "Liver", "Metabolic", "Lung", "MusculoSkeletal", "Age_prediction")
  system_PCA <- predict(RData$systems_PCA, system_scores)
  pred <- system_PCA %*% RData$Systems_clock_coefficients
  system_scores <- cbind(system_scores, pred)
  colnames(system_scores)[13] <- "SystemsAge"

  ## Scale system ages
  system_ages <- system_scores
  for (i in c(1:13)) {
    y <- system_ages[, i]
    system_ages[, i] <- (((y - RData$transformation_coefs[i, 1]) / RData$transformation_coefs[i, 2]) * RData$transformation_coefs[i, 4]) + RData$transformation_coefs[i, 3]
    system_ages[, i] <- system_ages[, i] / 12
  }

  if (is.null(pheno)) {
    results <- as.data.frame(system_ages)
  } else {
    results <- cbind(pheno, system_ages)
  }
  row.names(results) <- NULL
  return(results)
}
