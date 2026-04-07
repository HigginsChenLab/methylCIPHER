#' Calculate RepliTali
#'
#' @description Calculates the RepliTali mitotic age clock.
#'
#' @inheritParams param_template
#'
#' @inherit param_template return
#'
#' @export
#'
#' @examples
#' \dontrun{
#' calcRepliTali(exampleBetas)
#' }
calcRepliTali <- function(DNAm, pheno = NULL) {
  check_DNAm(DNAm)

  cpg_names <- RepliTali_CpGs$CpG

  DNAm <- impute_DNAm(
    DNAm = DNAm,
    method = "mean",
    CpGs = zero_cpgs(cpg_names),
    subset = TRUE
  )
  DNAm <- DNAm[, cpg_names, drop = FALSE]

  predictions <- as.numeric(
    as.matrix(DNAm) %*% RepliTali_CpGs$coefficient + RepliTali_ref$intercept
  )

  if (is.null(pheno)) {
    result <- data.frame(RepliTali = predictions)
    row.names(result) <- row.names(DNAm)
    return(result)
  }
  pheno$RepliTali <- predictions
  return(pheno)
}

#' Calculate RepliTaliNorm
#'
#' @description Calculates the RepliTaliNorm normalized mitotic age clock.
#'
#' @inheritParams param_template
#'
#' @inherit param_template return
#'
#' @export
#'
#' @examples
#' \dontrun{
#' calcRepliTaliNorm(exampleBetas)
#' }
calcRepliTaliNorm <- function(DNAm, pheno = NULL) {
  check_DNAm(DNAm)

  cpg_names <- RepliTaliNorm_CpGs$CpG

  DNAm <- impute_DNAm(
    DNAm = DNAm,
    method = "mean",
    CpGs = zero_cpgs(cpg_names),
    subset = TRUE
  )
  DNAm <- DNAm[, cpg_names, drop = FALSE]

  predictions <- as.numeric(
    as.matrix(DNAm) %*% RepliTaliNorm_CpGs$coefficient + RepliTaliNorm_ref$intercept
  )

  if (is.null(pheno)) {
    result <- data.frame(RepliTaliNorm = predictions)
    row.names(result) <- row.names(DNAm)
    return(result)
  }
  pheno$RepliTaliNorm <- predictions
  return(pheno)
}
