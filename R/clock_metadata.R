#' @title Clock Metadata
#'
#' @description Metadata for epigenetic clocks including authorship, training
#'   details, and array compatibility.
#'
#' @format A data frame with 92 rows and 10 variables:
#' \describe{
#'   \item{Clock Name}{Name of the epigenetic clock}
#'   \item{1st Author}{First author of the publication}
#'   \item{Year}{Publication year}
#'   \item{PMID}{PubMed ID}
#'   \item{Trained Phenotype}{The phenotype the clock was trained to predict}
#'   \item{# of CpGs}{Number of CpG probes used}
#'   \item{Cohort Trained}{Training cohort}
#'   \item{Tissues Derived}{Tissue type(s) used}
#'   \item{Age Range Trained}{Age range of training samples}
#'   \item{Array Type Trained}{Methylation array platform}
#' }
"clock_metadata"
