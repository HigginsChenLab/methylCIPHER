#' Calculate DNAmFitAge
#'
#' @description Calculates the DNAmFitAge biological fitness age indicator,
#'   including six DNAm fitness biomarkers (DNAmGait_noAge, DNAmGrip_noAge,
#'   DNAmVO2max, DNAmGait_wAge, DNAmGrip_wAge, DNAmFEV1_wAge) and the
#'   composite DNAmFitAge and FitAgeAcceleration.
#'
#' @param DNAm A matrix of methylation beta values with samples as rows and
#'   CpGs as columns.
#' @param pheno A data.frame with samples as rows. Must include columns:
#'   \itemize{
#'     \item \code{Female}: Numeric 0/1 (1 = female).
#'     \item \code{Age}: Numeric chronological age.
#'     \item \code{GrimAgeV1}: Pre-calculated GrimAgeV1 (from
#'       \code{calcGrimAgeV1}). If absent, it will be computed automatically.
#'   }
#'
#' @return The pheno data.frame with columns appended: DNAmGait_noAge,
#'   DNAmGrip_noAge, DNAmVO2max, DNAmGait_wAge, DNAmGrip_wAge,
#'   DNAmFEV1_wAge, DNAmFitAge, and FitAgeAcceleration.
#'
#' @references
#' McGreevy, K.M., et al. (2023). DNAmFitAge: biological age indicator
#' incorporating physical fitness. Aging, 15(10), 3904-3938.
#' \url{https://doi.org/10.18632/aging.204538}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' calcDNAmFitAge(exampleBetas, examplePheno)
#' }
calcDNAmFitAge <- function(DNAm, pheno) {
  # Input validation
  check_DNAm(DNAm)
  check_pheno(pheno, extra_columns = c("Female", "Age"))
  if (!"GrimAgeV1" %in% names(pheno)) {
    message("DNAmFitAge: GrimAgeV1 not found in pheno; computing GrimAgeV1 automatically.")
    pheno <- calcGrimAgeV1(DNAm, pheno)
  }

  # Ensure pheno and DNAm are aligned
  if (nrow(pheno) != nrow(DNAm)) {
    stop("pheno and DNAm must have the same number of rows.")
  }

  models <- DNAmFitAge_data

  # --- Imputation: fill missing CpGs with sex-specific medians ---
  all_cpgs <- models$AllCpGs
  present_cpgs <- intersect(all_cpgs, colnames(DNAm))
  missing_cpgs <- setdiff(all_cpgs, colnames(DNAm))

  # Sex indices (based on row position, shared between DNAm and pheno)
  fem_rows <- which(pheno$Female == 1)
  male_rows <- which(pheno$Female == 0)

  if (length(missing_cpgs) > 0) {
    message("DNAmFitAge: ", length(missing_cpgs),
            " missing CpGs imputed with sex-specific training medians.")

    # Build imputation matrix with sex-specific medians
    impute_matrix <- matrix(
      NA_real_,
      nrow = nrow(DNAm),
      ncol = length(missing_cpgs),
      dimnames = list(row.names(DNAm), missing_cpgs)
    )
    for (cpg in missing_cpgs) {
      if (length(fem_rows) > 0 && cpg %in% colnames(models$Female_Medians_All)) {
        impute_matrix[fem_rows, cpg] <- models$Female_Medians_All[[cpg]]
      }
      if (length(male_rows) > 0 && cpg %in% colnames(models$Male_Medians_All)) {
        impute_matrix[male_rows, cpg] <- models$Male_Medians_All[[cpg]]
      }
    }
    DNAm <- cbind(DNAm, impute_matrix)
  }

  # Mean-impute any remaining NAs in the required CpGs
  DNAm_sub <- DNAm[, all_cpgs, drop = FALSE]
  na_idx <- which(is.na(DNAm_sub), arr.ind = TRUE)
  if (nrow(na_idx) > 0) {
    col_means <- colMeans(DNAm_sub, na.rm = TRUE)
    DNAm_sub[na_idx] <- col_means[na_idx[, 2]]
  }

  # --- Calculate 6 fitness biomarkers ---

  # Helper: compute a single model's estimate for a set of sample indices
  # Some models include "Age" as a predictor alongside CpGs
  calc_model <- function(idx, model) {
    if (length(idx) == 0) return(numeric(0))
    terms <- model$term
    coefs <- model$estimate
    intercept <- coefs[1]
    pred_terms <- terms[-1]
    pred_coefs <- coefs[-1]

    # Build predictor matrix: CpGs from DNAm_sub, Age from pheno
    pred_mat <- matrix(NA_real_, nrow = length(idx), ncol = length(pred_terms),
                       dimnames = list(NULL, pred_terms))
    cpg_terms <- pred_terms[pred_terms %in% colnames(DNAm_sub)]
    if (length(cpg_terms) > 0) {
      pred_mat[, cpg_terms] <- DNAm_sub[idx, cpg_terms, drop = FALSE]
    }
    if ("Age" %in% pred_terms) {
      pred_mat[, "Age"] <- pheno$Age[idx]
    }

    as.numeric(pred_mat %*% pred_coefs + intercept)
  }

  # Initialize result columns
  n <- nrow(pheno)
  fitness_cols <- c("DNAmGait_noAge", "DNAmGrip_noAge", "DNAmVO2max",
                    "DNAmGait_wAge", "DNAmGrip_wAge", "DNAmFEV1_wAge")
  for (col in fitness_cols) pheno[[col]] <- NA_real_

  # Female estimates
  pheno$DNAmGait_noAge[fem_rows]  <- calc_model(fem_rows, models$Gait_noAge_Females)
  pheno$DNAmGrip_noAge[fem_rows]  <- calc_model(fem_rows, models$Grip_noAge_Females)
  pheno$DNAmVO2max[fem_rows]      <- calc_model(fem_rows, models$VO2maxModel)
  pheno$DNAmGait_wAge[fem_rows]   <- calc_model(fem_rows, models$Gait_wAge_Females)
  pheno$DNAmGrip_wAge[fem_rows]   <- calc_model(fem_rows, models$Grip_wAge_Females)
  pheno$DNAmFEV1_wAge[fem_rows]   <- calc_model(fem_rows, models$FEV1_wAge_Females)

  # Male estimates
  pheno$DNAmGait_noAge[male_rows] <- calc_model(male_rows, models$Gait_noAge_Males)
  pheno$DNAmGrip_noAge[male_rows] <- calc_model(male_rows, models$Grip_noAge_Males)
  pheno$DNAmVO2max[male_rows]     <- calc_model(male_rows, models$VO2maxModel)
  pheno$DNAmGait_wAge[male_rows]  <- calc_model(male_rows, models$Gait_wAge_Males)
  pheno$DNAmGrip_wAge[male_rows]  <- calc_model(male_rows, models$Grip_wAge_Males)
  pheno$DNAmFEV1_wAge[male_rows]  <- calc_model(male_rows, models$FEV1_wAge_Males)

  # --- Calculate DNAmFitAge (Klemera-Doubal weighted combination) ---
  pheno$DNAmFitAge <- NA_real_

  # Female FitAge
  fi <- fem_rows[complete.cases(pheno[fem_rows, c("Age", "DNAmGait_noAge",
                  "DNAmGrip_noAge", "DNAmVO2max", "GrimAgeV1")])]
  if (length(fi) > 0) {
    pheno$DNAmFitAge[fi] <-
      0.1044232 * ((pheno$DNAmVO2max[fi]      - 46.825091) / (-0.13620215)) +
      0.1742083 * ((pheno$DNAmGrip_noAge[fi]  - 39.857718) / (-0.22074456)) +
      0.2278776 * ((pheno$DNAmGait_noAge[fi]  -  2.508547) / (-0.01245682)) +
      0.4934908 * ((pheno$GrimAgeV1[fi]       -  7.978487) / ( 0.80928530))
  }

  # Male FitAge
  mi <- male_rows[complete.cases(pheno[male_rows, c("Age", "DNAmGait_noAge",
                  "DNAmGrip_noAge", "DNAmVO2max", "GrimAgeV1")])]
  if (length(mi) > 0) {
    pheno$DNAmFitAge[mi] <-
      0.1390346 * ((pheno$DNAmVO2max[mi]      - 49.836389) / (-0.141862925)) +
      0.1787371 * ((pheno$DNAmGrip_noAge[mi]  - 57.514016) / (-0.253179827)) +
      0.1593873 * ((pheno$DNAmGait_noAge[mi]  -  2.349080) / (-0.009380061)) +
      0.5228411 * ((pheno$GrimAgeV1[mi]       -  9.549733) / ( 0.835120557))
  }

  # FitAgeAcceleration: residuals of DNAmFitAge ~ Age
  complete_idx <- which(!is.na(pheno$DNAmFitAge))
  pheno$FitAgeAcceleration <- NA_real_
  if (length(complete_idx) > 1) {
    pheno$FitAgeAcceleration[complete_idx] <- residuals(
      lm(DNAmFitAge ~ Age, data = pheno[complete_idx, ])
    )
  }

  return(pheno)
}
