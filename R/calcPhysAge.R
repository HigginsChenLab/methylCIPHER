PhysAge_surrogates <- function(wf, label, DNAm) {
  intercept <- as.numeric(wf$Beta[1])
  w <- wf[-1, c("CpG", "Beta")]

  # keep only CpGs present in DNAm, and align Beta to those columns
  keep <- intersect(colnames(DNAm), w$CpG)
  if (length(keep) == 0L) {
    # no overlapping CpGs → return NAs for this score
    results <- data.frame(
      x1 = rep(NA_real_, times = nrow(DNAm)),
      x2 = rep(NA_real_, times = nrow(DNAm))
    )
  } else {
    w_aligned <- w$Beta[match(keep, w$CpG)]
    # multiply each CpG column by its weight, then take (non-NA) row means
    raw <- rowMeans(sweep(DNAm[, keep, drop = FALSE], 2, w_aligned, `*`), na.rm = TRUE)
    results <- data.frame(
      x1 = raw + intercept,
      x2 = raw
    )
  }
  names(results) <- paste0(label, c("_final", "_raw"))
  row.names(results) <- NULL
  return(results)
}

#' calcPhysAge
#'
#' @description A function to calculate DNAmPhysAge
#'
#' @inheritParams param_template
#'
#' @inherit param_template return
#'
#' @references
#' Please cite the source when using this clock:
#'
#' Arpawong, T.E., Hernandez, B., Potter, C., Leigh, R.J., Klopack, E.T., Hill, C., Fiorito, G., Smyth, L.J., O’Halloran, A.M., McGuinness, B., Faul, J.D., Kenny, R.A., McKnight, A.J, Crimmins, E.M., and McCrory, C, 2025. Physiological health Age (PhysAge): a novel multi-system molecular timepiece predicts health and mortality in older adults. GeroScience, pp.1-21. \url{https://doi.org/10.1007/s11357-025-01832-1}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' PhysAge <- calcPhysAge(
#'   exampleBetas,
#'   examplePheno
#' )
#' PhysAge
#' }
calcPhysAge <- function(DNAm, pheno = NULL) {
  # Input validation
  check_DNAm(DNAm)
  if (is.null(pheno)) {
    pheno <- data.frame(row.names = row.names(DNAm))
  }
  pheno_cols <- names(pheno)

  ## Imputation
  CpGs <- PhysAge_CpGs$mean
  names(CpGs) <- PhysAge_CpGs$CpG

  DNAm <- impute_DNAm(
    DNAm = DNAm,
    method = "mean",
    CpGs = CpGs,
    subset = TRUE
  )

  ## Re-align to make sure things lined up with the object
  DNAm <- DNAm[, PhysAge_CpGs$CpG, drop = F]

  # Calculate components of PhysAge
  for (i in seq_along(PhysAge_data)) {
    pheno <- cbind(pheno, PhysAge_surrogates(PhysAge_data[[i]], names(PhysAge_data)[i], DNAm))
  }

  # Zscores
  raw_vars <- paste0(names(PhysAge_data), "_raw")
  neg_vars <- c("DNAmPeakflow_raw", "DNAmDHEAS_raw", "DNAmHDL_raw")
  data_to_scale <- as.matrix(pheno[, raw_vars])
  ## Negate
  data_to_scale[, neg_vars] <- data_to_scale[, neg_vars] * -1
  ## Scale and add z-scores to pheno
  z_scores <- scale(data_to_scale)
  colnames(z_scores) <- paste0("z", names(PhysAge_data))
  pheno <- cbind(pheno, as.data.frame(z_scores))

  # DNAmPhysAge
  pheno[["DNAmPhysAge"]] <- rowSums(z_scores, na.rm = TRUE)
  pheno[["zDNAmPhysAge"]] <- scale(pheno[["DNAmPhysAge"]])[, 1]
  pheno[["DNAmPhysAge_years"]] <- (pheno[["zDNAmPhysAge"]] * 9.630436) + 68.14726

  # Interpretable scale
  year_scores <- (z_scores * 9.630436) + 68.14726
  colnames(year_scores) <- paste0(names(PhysAge_data), "_years")
  pheno <- cbind(pheno, as.data.frame(year_scores))

  # Sort
  other_cols <- sort(setdiff(names(pheno), pheno_cols))
  pheno <- pheno[, c(pheno_cols, other_cols)]

  return(pheno)
}
