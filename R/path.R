# Copied from cmdstanr github path.R


#' Get Default Path for methylCIPHER Data
#'
#' This function determines the default directory path where methylCIPHER data,
#' required for SystemsAge and PCClocks, should be stored.
#'
#' @return The default path to the methylCIPHER data directory.
#' @export
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

#' Set Custom Path for methylCIPHER Data
#'
#' This function allows the user to set a custom directory path for storing
#' methylCIPHER data, which is required for SystemsAge and PCClocks. The specified
#' path is stored in R's options under "methylCIPHER.path".
#'
#' @param path A string specifying the custom path to the methylCIPHER data directory.
#'
#' @return None
#' @export
#'
#' @examples
#' get_methylCIPHER_path() # Default Path
#' set_methylCIPHER_path(".") # Change to "."
#' get_methylCIPHER_path() # New Path
set_methylCIPHER_path <- function(path) {
  options(methylCIPHER.path = path)
}

#' Get Path for methylCIPHER Data
#'
#' This function retrieves the path to the methylCIPHER data directory, which is
#' required for SystemsAge and PCClocks. It first checks if a custom path has been
#' set using [set_methylCIPHER_path()] or [options()]. If no custom path is found, it returns
#' the default path obtained from [get_methylCIPHER_default_path()].
#'
#' @return The path to the methylCIPHER data directory.
#' @export
#'
#' @examples
#' get_methylCIPHER_path() # Default Path
#' set_methylCIPHER_path(".") # Change to "."
#' get_methylCIPHER_path() # New Path
get_methylCIPHER_path <- function() {
  path <- getOption("methylCIPHER.path")
  if(is.null(path)) {
    return(get_methylCIPHER_default_path())
  }
  return(path)
}
