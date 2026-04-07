#' Remove All-NA CpG Probes
#'
#' @description Removes CpG columns where all samples are NA.
#'
#' @param DNAm A matrix of methylation beta values with samples as rows and
#'   CpGs as columns.
#'
#' @return The DNAm matrix with all-NA columns removed.
#' @export
#'
#' @examples
#' \dontrun{
#' DNAm <- removeNAprobes(DNAm)
#' }
removeNAprobes <- function(DNAm) {
  all_na_cols <- colSums(is.na(DNAm)) == nrow(DNAm)
  if (any(all_na_cols)) {
    message(paste0("Removing ", sum(all_na_cols), " all-NA CpG columns."))
    DNAm <- DNAm[, !all_na_cols]
  }
  return(DNAm)
}
