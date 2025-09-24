## code to prepare `large_clocks_gdrive` dataset goes here
## These codes define the clocks that need to be downloaded externally because the
## size of the weights are large.

## Add
library(googledrive)
clocks <- c("SystemsAge_data.qs2", "PCClocks_data.qs2", "methylCIPHER_test.csv")

drive_auth()

large_clocks_gdrive <- drive_get(clocks)
names(large_clocks_gdrive$id) <- large_clocks_gdrive$name

large_clocks_gdrive$clock <- stringr::str_remove(large_clocks_gdrive$name, "_data.qs2")
large_clocks_gdrive$name <- NULL
large_clocks_gdrive$drive_resource <- NULL

drive_deauth()

# library(zen4R)
# zenodo <- ZenodoManager$new(logger = "INFO")
# rec <- zenodo$getRecordByDOI("10.5281/zenodo.17162604")
# rec$listFiles(pretty = TRUE)

large_clocks_data <- list(
  "googledrive" = large_clocks_gdrive
  # , "zenodo" = rec
)

usethis::use_data(large_clocks_data, internal = TRUE, overwrite = TRUE)
