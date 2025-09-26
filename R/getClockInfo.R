#' getClockProbes
#'
#' @param DNAm The methylation Beta values you will calculate epigenetic clocks with, where columns are CpGs, and rows are samples.
#'
#' @return A table which compares the number of probes required by the full clocks in the current package, and the ones in your data. If you are missing many probes for a clock, imputation is necessary, though you might consider using a different clock or dataset if possible.
#' @export
#'
#' @examples getClockProbes(exampleBetas)
getClockInfo <- function() {
  library(googledrive)
  library(googlesheets4)

  gs4_auth()

  home_dir <- Sys.getenv("HOME")

  sheet_id <- "1fqxvyntgDX4AcQLQobCRNeOZjs2TqaigVwn6LuM6tMo"

  df <- read_sheet(sheet_id)

  df <- df[df$MethylCIPHER == "Public", ]

  df$TranslAGE <- NULL
  df$MethylCIPHER <- NULL
  df$`Misc Information about Clock` <- NULL
  df$`Friendly Name` <- NULL
  df$`Alternative Names` <- NULL

  ProbeTable <- as.data.frame(df)

  ProbeTable
}
