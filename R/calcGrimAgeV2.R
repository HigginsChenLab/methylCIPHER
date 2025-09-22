#' calcGrimAgeV2
#'
#' @description A function to calculate GrimAgeV2
#'
#' @inheritParams param_template_female_age
#'
#' @inherit param_template_female_age return
#'
#' @export
#'
#' @examples
#' \dontrun{
#' GrimAgeV2 <- calcGrimAgeV2(
#'   exampleBetas,
#'   examplePheno
#' )
#' GrimAgeV2
#' }
calcGrimAgeV2 <- function(DNAm, pheno, ID = "Sample_ID") {
  # Input validation
  # Check DNAm
  check_DNAm(DNAm)
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

  ## Imputation
  DNAm <- impute_DNAm(
    DNAm = DNAm,
    method = "mean",
    CpGs = zero_cpgs(CpGs_GrimAge2),
    subset = TRUE
  )

  ## Re-align to make sure things lined up with the object
  DNAm <- DNAm[, CpGs_GrimAge2, drop = F]

  # Calculate GrimAge Clocks
  DNAm_AF <- cbind(DNAm, as.matrix(subset(pheno, select = c("Age"))))
  pheno$DNAmPACKYRS <- as.numeric(DNAm_AF[, names(CalcGrimAge2$PACKYRS.model)] %*% CalcGrimAge2$PACKYRS.model + CalcGrimAge2$PACKYRS.intercept)
  pheno$DNAmADM <- as.numeric(DNAm_AF[, names(CalcGrimAge2$ADM.model)] %*% CalcGrimAge2$ADM.model + CalcGrimAge2$ADM.intercept)
  pheno$DNAmB2M <- as.numeric(DNAm_AF[, names(CalcGrimAge2$B2M.model)] %*% CalcGrimAge2$B2M.model + CalcGrimAge2$B2M.intercept)
  pheno$DNAmCystatinC <- as.numeric(DNAm_AF[, names(CalcGrimAge2$CystatinC.model)] %*% CalcGrimAge2$CystatinC.model + CalcGrimAge2$CystatinC.intercept)
  pheno$DNAmGDF15 <- as.numeric(DNAm_AF[, names(CalcGrimAge2$GDF15.model)] %*% CalcGrimAge2$GDF15.model + CalcGrimAge2$GDF15.intercept)
  pheno$DNAmLeptin <- as.numeric(DNAm_AF[, names(CalcGrimAge2$Leptin.model)] %*% CalcGrimAge2$Leptin.model + CalcGrimAge2$Leptin.intercept)
  pheno$DNAmPAI1 <- as.numeric(DNAm_AF[, names(CalcGrimAge2$PAI1.model)] %*% CalcGrimAge2$PAI1.model + CalcGrimAge2$PAI1.intercept)
  pheno$DNAmTIMP1 <- as.numeric(DNAm_AF[, names(CalcGrimAge2$TIMP1.model)] %*% CalcGrimAge2$TIMP1.model + CalcGrimAge2$TIMP1.intercept)
  # New Proteins
  pheno$DNAmlogA1C <- as.numeric(DNAm_AF[, names(CalcGrimAge2$logA1C.model)] %*% CalcGrimAge2$logA1C.model + CalcGrimAge2$logA1C.intercept)
  pheno$DNAmlogCRP <- as.numeric(DNAm_AF[, names(CalcGrimAge2$logCRP.model)] %*% CalcGrimAge2$logCRP.model + CalcGrimAge2$logCRP.intercept)
  # Calculate GrimAgeV2
  pheno$GrimAgeV2 <- as.numeric(as.matrix(subset(pheno, select = CalcGrimAge2$components)) %*% CalcGrimAge2$GrimAge.model) #+ CalcGrimAge2$GrimAge.intercept)
  y <- pheno$GrimAgeV2
  pheno$GrimAgeV2 <- (((y - CalcGrimAge2$GrimAge.transform[3]) / CalcGrimAge2$GrimAge.transform[4]) * CalcGrimAge2$GrimAge.transform[2]) + CalcGrimAge2$GrimAge.transform[1]

  return(pheno)
}
