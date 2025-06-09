# Strategy
# Create a path of the folder that contains the downloaded file
#   check the folder to see if its a folder
#   check the available objects and see if everythings downloaded
#   checksum potentially
# When the calcSystemsAge gets called the first time, objects are read in based
# on the path or list of dependencies
# Then proceed as normal

get_methylCIPHER_default_path <- function() {
  return(paste(Sys.getenv("HOME"), "methylCIPHER", sep = "/"))
}

set_methylCIPHER_path <- function(path) {
  options(methylCIPHER.path = path)
}

get_methylCIPHER_path <- function() {
  path <- getOption("methylCIPHER.path")
  if(is.null(path)) {
    return(get_methylCIPHER_default_path())
  }
  return(path)
}
