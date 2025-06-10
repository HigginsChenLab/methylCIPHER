# Strategy
# Create a path of the folder that contains the downloaded file
#   check the folder to see if its a folder
#   check the available objects and see if everythings downloaded
#   checksum potentially
# When the calcSystemsAge gets called the first time, objects are read in based
# on the path or list of dependencies
# Then proceed as normal

# Copied from cmdstanr github path.R
get_methylCIPHER_default_path <- function() {
  home <- Sys.getenv("HOME")
  if (.Platform$OS.type == "windows") {
    userprofile <- Sys.getenv("USERPROFILE")
    h_drivepath <- file.path(Sys.getenv("HOMEDRIVE"), Sys.getenv("HOMEPATH"))
    win_home <- ifelse(userprofile == "", h_drivepath, userprofile)
    if (win_home != "") {
      home <- win_home
    }
  }
  return(paste0(home, "/", "methylCIPHER"))
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
