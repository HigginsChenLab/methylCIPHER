check_DNAm <- function(DNAm, missing_allowed = TRUE) {
  checkmate::assert_matrix(DNAm, mode = "double", any.missing = missing_allowed, min.rows = 1, min.cols = 1)
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
      stop("CpGs with all NA detected.")
    }
    if (any(rowSums(missing) == ncol(missing))) {
      stop("Samples with all NA detected.")
    }
    if (any(missing)) {
      warning("Missing values in DNAm detected.")
    }
  }

  invisible(TRUE)
}
