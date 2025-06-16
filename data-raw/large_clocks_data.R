## code to prepare `large_clocks_data` dataset goes here
## These codes define the clocks that need to be downloaded externally because the
## size of the weights are large.

## Add
library(googledrive)
clocks <- c("SystemsAge_data.qs2", "PCClocks_data.qs2")

files <- drive_get(clocks)
names(files$id) <- files$name

drive_auth()

large_clocks_data <- as.data.frame(
  rbind(
    c(
      clock = "SystemsAge",
      url = files$id[["SystemsAge_data.qs2"]],
      download_name = "SystemsAge_data.qs2",
      type = "googledrive"
    ),
    c(
      clock = "PCClocks",
      url = files$id[["PCClocks_data.qs2"]],
      download_name = "PCClocks_data.qs2",
      type = "googledrive"
    )
  )
)

drive_deauth()
usethis::use_data(large_clocks_data, internal = TRUE, overwrite = TRUE)
