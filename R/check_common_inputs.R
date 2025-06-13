#' Check DNA Methylation Data Matrix
#'
#' This function performs a series of checks on the DNA methylation data matrix to ensure it is properly formatted and contains valid data.
#'
#' @inheritParams param_template
#' @param missing_allowed A logical value indicating whether missing values are allowed in the matrix. Default is \code{TRUE}.
#'
#' @details
#' The function performs the following checks and stops with an error if any of the following conditions are not met:
#' \itemize{
#'   \item \code{DNAm} is a matrix with double values, with at least one row and one column.
#'   \item If \code{missing_allowed} is \code{FALSE}, \code{DNAm} must not contain any missing values.
#'   \item If \code{missing_allowed} is \code{TRUE}, \code{DNAm} may contain missing values, but no column or row can be entirely \code{NA}.
#'   \item Column names of \code{DNAm} must be unique character strings.
#' }
#' Additionally, the function issues warnings for the following conditions:
#' \itemize{
#'   \item If row names of \code{DNAm} are not set, as they should typically represent sample names.
#'   \item If there are more samples (rows) than cpgs (column), which is unusual for DNA methylation data.
#'   \item No column names start with "cg", suggesting that the matrix might need to be transposed.
#'   \item If missing values are present and \code{missing_allowed} is \code{TRUE}, indicating that imputation will be handled by \code{impute_DNAm}.
#' }
#'
#' @return Returns \code{invisible(TRUE)} if all checks pass.
#' @keywords internal
check_DNAm <- function(DNAm, missing_allowed = TRUE) {
  checkmate::assert_matrix(
    DNAm,
    mode = "double",
    any.missing = missing_allowed,
    min.rows = 1,
    min.cols = 1
  )
  # CpGs names has to be colnames
  checkmate::assert_character(colnames(DNAm), unique = TRUE, null.ok = FALSE)

  if(is.null(row.names(DNAm))) {
    warning("DNAm sample names as row.names are not detected.")
  }
  if(nrow(DNAm) > ncol(DNAm)) {
    warning("DNAm should be formatted in samples * CpG. Currently, DNAm have more rows (samples) than columns (CpGs) which is highly unlikely.")
  }
  if(!any(grepl("^cg", colnames(DNAm)))) {
    warning("Warning: It looks like you may need to format DNAm using t(DNAm) to get samples as rows!")
  }

  if (missing_allowed) {
    missing <- is.na(DNAm)
    if (any(colSums(missing) == nrow(missing))) {
      stop("CpGs with all NA are not allowed.")
    }
    if (any(rowSums(missing) == ncol(missing))) {
      stop("Samples with all NA are not allowed.")
    }
    if (any(missing)) {
      warning("Missing values in DNAm detected. Imputation handled by `impute_DNAm`.")
    }
  }

  invisible(TRUE)
}

#' Check Phenotype Data Frame
#'
#' This function checks the phenotype data frame (`pheno`) to ensure it contains the required columns and that the data types are correct.
#'
#' @inheritParams param_template
#' @param extra_columns An optional character vector specifying additional columns to check.
#'
#' @details
#' The function performs the following checks and stops with an error if any of the following conditions are not met:
#' \itemize{
#'   \item \code{pheno} is a data frame with at least one row.
#'   \item \code{pheno} contains a column named "Sample_ID", which must be a character vector without missing values.
#'   \item If "Female" is specified in \code{extra_columns}, \code{pheno} must contain a column named "Female", which must be an integerish vector with values 0 or 1 and no missing values.
#'   \item If "Age" is specified in \code{extra_columns}, \code{pheno} must contain a column named "Age", which must be a numeric vector with finite values and no missing values.
#' }
#'
#' @return The function performs the checks and stops with an error if any check fails.
#'
#' @keywords internal
check_pheno <- function(pheno, extra_columns = NULL) {
  checkmate::assert_data_frame(pheno, min.rows = 1, null.ok = FALSE)
  # check `Sample_ID`. By default this column must be here if pheno is ever added.
  stopifnot("`pheno` must contain `Sample_ID` column" = "Sample_ID" %in% names(pheno))
  checkmate::assert_character(
    pheno[["Sample_ID"]],
    any.missing = FALSE,
    null.ok = FALSE,
  )
  if("Female" %in% extra_columns) {
    checkmate::assert_integerish(
      pheno[["Female"]],
      lower = 0,
      upper = 1,
      null.ok = FALSE,
      any.missing = FALSE
    )
  }
  if("Age" %in% extra_columns) {
    checkmate::assert_numeric(
      pheno[["Age"]],
      finite = TRUE,
      null.ok = FALSE,
      any.missing = FALSE
    )
  }

  invisible(TRUE)
}

#' Align Pheno
#'
#' Filter and sort pheno given a list of ids. ids are guaranteed to be in pheno.
#'
#' @inheritParams param_template
#'
#' @examples
#' \dontrun{
#' examplePheno$Sample_ID <- examplePheno$name
#' align_pheno(examplePheno, c("7786915023_R02C02", "7786915023_R02C02"))
#' }
#'
#' @keywords internal
align_pheno <- function(pheno, Sample_ID) {
  indices <- match(Sample_ID, pheno$Sample_ID)
  return(pheno[indices, , drop = F])
}
