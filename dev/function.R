calcGrimAgeV1_1 <- function(DNAm, pheno) {
  # Input validation
  # Check DNAm
  check_DNAm(DNAm)
  # Check Pheno
  check_pheno(pheno, extra_columns = c("Female", "Age"))
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
    if(length(samples) == 0) {
      stop("DNAm and pheno have no Sample_ID in common")
    }
    DNAm <- DNAm[samples, , drop = FALSE]
    pheno <- align_pheno(pheno, samples)
    stopifnot("`DNAm` and `pheno` samples alignment failed. Check `Sample_ID` of pheno and row.names() of `DNAm`" = all.equal(row.names(DNAm), pheno$Sample_ID))
    message("Samples inconsistencies between DNAm and Pheno were detected and corrected.")
  }

  ## Imputation
  DNAm <- impute_DNAm(
    DNAm = DNAm,
    method = "mean",
    CpGs = zero_cpgs(CpGs_GrimAge1),
    subset = TRUE
  )

  ## Re-align to make sure things lined up with the object
  DNAm <- DNAm[, CpGs_GrimAge1, drop = F]

  # Calculate GrimAge Clocks
  temp <- as.matrix(cbind(DNAm, Age = as.numeric(pheno$Age, Female = pheno$Female)))
  pheno$DNAmPACKYRS <- as.numeric(temp[, names(CalcGrimAgeV1$PACKYRS.model)] %*% CalcGrimAgeV1$PACKYRS.model + CalcGrimAgeV1$PACKYRS.intercept)
  pheno$DNAmADM <- as.numeric(temp[, names(CalcGrimAgeV1$ADM.model)] %*% CalcGrimAgeV1$ADM.model + CalcGrimAgeV1$ADM.intercept)
  pheno$DNAmB2M <- as.numeric(temp[, names(CalcGrimAgeV1$B2M.model)] %*% CalcGrimAgeV1$B2M.model + CalcGrimAgeV1$B2M.intercept)
  pheno$DNAmCystatinC <- as.numeric(temp[, names(CalcGrimAgeV1$CystatinC.model)] %*% CalcGrimAgeV1$CystatinC.model + CalcGrimAgeV1$CystatinC.intercept)
  pheno$DNAmGDF15 <- as.numeric(temp[, names(CalcGrimAgeV1$GDF15.model)] %*% CalcGrimAgeV1$GDF15.model + CalcGrimAgeV1$GDF15.intercept)
  pheno$DNAmLeptin <- as.numeric(temp[, names(CalcGrimAgeV1$Leptin.model)] %*% CalcGrimAgeV1$Leptin.model + CalcGrimAgeV1$Leptin.intercept)
  pheno$DNAmPAI1 <- as.numeric(temp[, names(CalcGrimAgeV1$PAI1.model)] %*% CalcGrimAgeV1$PAI1.model + CalcGrimAgeV1$PAI1.intercept)
  pheno$DNAmTIMP1 <- as.numeric(temp[, names(CalcGrimAgeV1$TIMP1.model)] %*% CalcGrimAgeV1$TIMP1.model + CalcGrimAgeV1$TIMP1.intercept)
  pheno$Age <- pheno$Age
  pheno$Female <- pheno$Female
  pheno$GrimAgeV1 <- as.numeric(as.matrix(pheno[, CalcGrimAgeV1$components]) %*% CalcGrimAgeV1$GrimAge.model) #+ CalcGrimAgeV1$GrimAge.intercept)
  y <- pheno$GrimAgeV1
  pheno$GrimAgeV1 <- (((y - CalcGrimAgeV1$GrimAge.transform[3]) / CalcGrimAgeV1$GrimAge.transform[4]) * CalcGrimAgeV1$GrimAge.transform[2]) + CalcGrimAgeV1$GrimAge.transform[1]

  return(pheno)
}
