#' Calculate IntrinClock
#'
#' @description Calculates the IntrinClock epigenetic age estimate using a
#'   glmnet model with developmental age transformation.
#'
#' @inheritParams param_template
#'
#' @inherit param_template return
#'
#' @export
#'
#' @examples
#' \dontrun{
#' calcIntrinClock(exampleBetas, examplePheno)
#' }
calcIntrinClock <- function(DNAm, pheno = NULL) {
  # Input validation
  check_DNAm(DNAm)

  # Load bundled model
  model <- IntrinClock_data$model
  imputation_ref <- IntrinClock_data$imputation
  allowed_cpgs <- rownames(stats::coef(model))[-1]

  # Build imputation reference, falling back to 0
  if (is.null(imputation_ref)) {
    impute_vals <- zero_cpgs(allowed_cpgs)
  } else {
    missing_from_ref <- setdiff(allowed_cpgs, names(imputation_ref))
    impute_vals <- imputation_ref
    if (length(missing_from_ref) > 0) {
      impute_vals <- c(impute_vals, zero_cpgs(missing_from_ref))
    }
  }
  impute_vals <- impute_vals[allowed_cpgs]

  # Imputation
  DNAm <- impute_DNAm(
    DNAm = DNAm,
    method = "mean",
    CpGs = impute_vals,
    subset = TRUE
  )

  # Re-align to make sure things lined up with the object
  DNAm <- DNAm[, allowed_cpgs, drop = FALSE]

  # Extract coefficients directly from model internals
  lambda_idx <- which(model$glmnet.fit$lambda == model$lambda.min)
  intercept <- model$glmnet.fit$a0[lambda_idx]
  betas_coef <- as.numeric(model$glmnet.fit$beta[, lambda_idx])
  ages <- as.matrix(DNAm) %*% betas_coef + intercept
  intrin_ages <- anti.trafo(ages)

  if (is.null(pheno)) {
    result <- data.frame(IntrinClock = as.numeric(intrin_ages))
    row.names(result) <- row.names(DNAm)
    return(result)
  }
  pheno$IntrinClock <- as.numeric(intrin_ages)
  return(pheno)
}
