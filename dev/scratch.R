
load_all()

Sys.getenv("HOME")
get_methylCIPHER_path()
list.files(get_methylCIPHER_path())

my_list <- list()
my_env <- new.env()
load(paste0(get_methylCIPHER_path(), "/", "SystemsAge_data.RData"), envir = my_env)

debug(calcSystemsAge)

as.matrix(exampleBetas) |> class()
class()
example_SystemsAge <- calcSystemsAge(
   exampleBetas,
   examplePheno,
   imputation = F,
   # If `SystemsAge_data.RData` has been downloaded to default folder at $HOME
   RData = get_methylCIPHER_path()
)

?calcSystemsAge

impute_DNAm <- function(DNAm, method = c("mean"), CpGs = NULL) {
  suppressWarnings(check_DNAm(DNAm))
  method <- match.arg(method)
  checkmate::assert_character(CpGs, min.len = 1, names = "unique", null.ok = TRUE)
  if(is.null(CpGs)) {
    CpGs <- numeric(length = ncol(DNAm))
    names(CpGs) <- colnames(DNAm)
  }

  # Completely missing CpG
  needed_cpgs <- setdiff(names(CpGs), colnames(DNAm))

  needed_matrix <- matrix(
    CpGs[needed_cpgs],
    ncol = length(needed_cpgs),
    nrow = nrow(DNAm),
    byrow = TRUE,
    dimnames = list(row.names(DNAm), needed_cpgs)
  )

  # Imputed CpG
  matched_cpgs <- intersect(names(CpGs), colnames(DNAm))
  imputed_matrix <- if(sum(is.na(DNAm[, matched_cpgs])) == 0) {
    DNAm[, matched_cpgs]
  } else {
    meanimpute1(DNAm[, matched_cpgs])
  }

  return(cbind(imputed_matrix, matched_cpgs))
}


impute_DNAm(as.matrix(exampleBetas), CpGs = my_env$imputeMissingCpGs[1:2])
