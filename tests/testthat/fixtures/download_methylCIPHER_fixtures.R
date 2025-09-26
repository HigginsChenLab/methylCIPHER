library(googledrive)

drive_auth()

file_name <- "methylCIPHER_test.csv"
path <- file.path(withr::local_tempdir(), file_name)

set.seed(1234)
test <- data.frame(X1 = rnorm(5))
write.csv(test, path, row.names = FALSE)

folder <- drive_find(q = "name = 'shared_objects' and mimeType = 'application/vnd.google-apps.folder'")
if (nrow(folder) == 0) {
  folder <- drive_mkdir("shared_objects")
}

drive_upload(
  media = path,
  path = folder$id,
  name = file_name,
  overwrite = TRUE
)

drive_get(file_name)

drive_deauth()
unlink(path)
