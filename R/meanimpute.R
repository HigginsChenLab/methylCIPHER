#' Mean imputation for Beta Matrix
#'
#' Performs column mean imputation for missing values in a methylation beta matrix, where rows represent samples and columns represent CpGs.
#'
#' @inheritParams param_template
#' @return A matrix of methylation beta values with missing values imputed using column means.
#'
#' @keywords internal
#'
mean_impute <- function(DNAm) {
  na_indices <- which(is.na(DNAm), arr.ind = TRUE)
  column_means <- colMeans(DNAm, na.rm = TRUE)
  DNAm[na_indices] <- column_means[na_indices[, 2]]
  return(DNAm)
}

meanimpute <- function(x){
  apply(x,2,function(z)ifelse(is.na(z),mean(z,na.rm=T),z))
}
