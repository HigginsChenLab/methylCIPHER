#' Calculate PC Clocks
#'
#' @description
#' Calculates PC Clocks
#'
#' @inheritParams param_template_female_age
#' @param RData
#' Default to `NULL`, which is the default path to the downloaded data folder.
#' Either a character string specifying the path to the folder containing
#' the `PCClocks_data.qs2` file, or a list containing the contents of the
#' `PCClocks_data.qs2` file loaded via [load_PCClocks_data()]. See Details.
#'
#' @details
#' Systems Age calculation requires the `PCClocks_data.qs2` object, which can be
#' downloaded using [download_methylCIPHER()]. Provide either the path to the
#' downloaded `PCClocks_data.qs2` file or the object loaded with
#' [load_PCClocks_data()] to the `RData` parameter. The advantage of passing
#' in the loaded data is that if this function were to be run multiple times,
#' then the object won't be loaded in each time.
#'
#' @inherit param_template_female_age return
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Download the external data
#' download_methylCIPHER(clocks = "PCClocks")
#'
#' # Either path to the data
#' RData <- get_methylCIPHER_path()
#' # Or load in the object
#' RData <- load_PCClocks_data(get_methylCIPHER_path())
#'
#' # Calculate Systems Age using example data
#' result <- calcPCClocks(
#'   betas = exampleBetas,
#'   pheno = examplePheno,
#'   RData = RData
#' )
#' result
#' }
calcPCClocks <- function(DNAm, pheno, ID = "Sample_ID", RData = NULL) {
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
  check_pheno(pheno, ID = ID, extra_columns = c("Female", "Age"))
  # Check Consistent between `pheno` and `DNAm`
  need_align <- !isTRUE(all.equal(row.names(DNAm), pheno[[ID]]))
  if (need_align) {
    samples <- intersect(row.names(DNAm), pheno[[ID]])
    if (length(samples) == 0) {
      stop("DNAm and pheno have no ID in common.")
    }
    DNAm <- DNAm[samples, , drop = FALSE]
    pheno <- align_pheno(pheno, samples, ID = ID)
    stopifnot("`DNAm` and `pheno` samples alignment failed. Check ID of pheno and row.names() of `DNAm`" = all.equal(row.names(DNAm), pheno[[ID]]))
    message("Samples inconsistencies between DNAm and Pheno were detected and corrected.")
  }
  # handle RData
  if (is.null(RData)) {
    RData <- load_PCClocks_data()
  } else if (is.character(RData)) {
    RData <- load_PCClocks_data(RData)
  }

  if (rlang::hash(RData) != "46386ec4be2b2a5239cf67b242d7dc24") {
    stop("The downloaded PCClocks data is corrupted or the wrong data (e.g., SystemsAge) was passed. See `?download_methylCIPHER()`.")
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

  message("Calculating PC Clocks now")

  # Calculate PC Clocks
  pheno$PCHorvath1 <- as.numeric(anti.trafo(sweep(DNAm, 2, RData$CalcPCHorvath1$center) %*% RData$CalcPCHorvath1$rotation %*% RData$CalcPCHorvath1$model + RData$CalcPCHorvath1$intercept))
  pheno$PCHorvath2 <- as.numeric(anti.trafo(sweep(DNAm, 2, RData$CalcPCHorvath2$center) %*% RData$CalcPCHorvath2$rotation %*% RData$CalcPCHorvath2$model + RData$CalcPCHorvath2$intercept))
  pheno$PCHannum <- as.numeric(sweep(DNAm, 2, RData$CalcPCHannum$center) %*% RData$CalcPCHannum$rotation %*% RData$CalcPCHannum$model + RData$CalcPCHannum$intercept)
  pheno$PCPhenoAge <- as.numeric(sweep(DNAm, 2, RData$CalcPCPhenoAge$center) %*% RData$CalcPCPhenoAge$rotation %*% RData$CalcPCPhenoAge$model + RData$CalcPCPhenoAge$intercept)
  pheno$PCDNAmTL <- as.numeric(sweep(DNAm, 2, RData$CalcPCDNAmTL$center) %*% RData$CalcPCDNAmTL$rotation %*% RData$CalcPCDNAmTL$model + RData$CalcPCDNAmTL$intercept)
  DNAm <- cbind(sweep(DNAm, 2, RData$CalcPCGrimAge$center) %*% RData$CalcPCGrimAge$rotation, Female = pheno$Female, Age = pheno$Age)
  pheno$PCPACKYRS <- as.numeric(DNAm[, names(RData$CalcPCGrimAge$PCPACKYRS.model)] %*% RData$CalcPCGrimAge$PCPACKYRS.model + RData$CalcPCGrimAge$PCPACKYRS.intercept)
  pheno$PCADM <- as.numeric(DNAm[, names(RData$CalcPCGrimAge$PCADM.model)] %*% RData$CalcPCGrimAge$PCADM.model + RData$CalcPCGrimAge$PCADM.intercept)
  pheno$PCB2M <- as.numeric(DNAm[, names(RData$CalcPCGrimAge$PCB2M.model)] %*% RData$CalcPCGrimAge$PCB2M.model + RData$CalcPCGrimAge$PCB2M.intercept)
  pheno$PCCystatinC <- as.numeric(DNAm[, names(RData$CalcPCGrimAge$PCCystatinC.model)] %*% RData$CalcPCGrimAge$PCCystatinC.model + RData$CalcPCGrimAge$PCCystatinC.intercept)
  pheno$PCGDF15 <- as.numeric(DNAm[, names(RData$CalcPCGrimAge$PCGDF15.model)] %*% RData$CalcPCGrimAge$PCGDF15.model + RData$CalcPCGrimAge$PCGDF15.intercept)
  pheno$PCLeptin <- as.numeric(DNAm[, names(RData$CalcPCGrimAge$PCLeptin.model)] %*% RData$CalcPCGrimAge$PCLeptin.model + RData$CalcPCGrimAge$PCLeptin.intercept)
  pheno$PCPAI1 <- as.numeric(DNAm[, names(RData$CalcPCGrimAge$PCPAI1.model)] %*% RData$CalcPCGrimAge$PCPAI1.model + RData$CalcPCGrimAge$PCPAI1.intercept)
  pheno$PCTIMP1 <- as.numeric(DNAm[, names(RData$CalcPCGrimAge$PCTIMP1.model)] %*% RData$CalcPCGrimAge$PCTIMP1.model + RData$CalcPCGrimAge$PCTIMP1.intercept)
  pheno$PCGrimAge <- as.numeric(as.matrix(subset(pheno, select = RData$CalcPCGrimAge$components)) %*% RData$CalcPCGrimAge$PCGrimAge.model + RData$CalcPCGrimAge$PCGrimAge.intercept)

  message("PC Clocks successfully calculated!")

  return(pheno)
}
