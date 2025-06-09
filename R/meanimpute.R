#' Mean Imputation for singly missing CpGs
#'
#' @param x A dataframe of CpG Betas with missing values
#'
#' @return Mean Imputed dataframe
#' @export
meanimpute <- function(x){
  apply(x,2,function(z)ifelse(is.na(z),mean(z,na.rm=T),z))
  # na_indices <- which(is.na(x), arr.ind = TRUE)
  # stopifnot("Check input with check_DNAm" = nrow(na_indices) > 0)
  # column_means <- colMeans(x, na.rm = TRUE)
  # x[na_indices] <- column_means[na_indices[, 2]]
  # return(x)
}

meanimpute1 <- function(x){
  na_indices <- which(is.na(x), arr.ind = TRUE)
  stopifnot("Check input with check_DNAm" = nrow(na_indices) > 0)
  column_means <- colMeans(x, na.rm = TRUE)
  x[na_indices] <- column_means[na_indices[, 2]]
  return(x)
}
