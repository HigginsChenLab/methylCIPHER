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
  return(file.path(home, "methylCIPHER"))
}

#' Get Path for methylCIPHER Data
#'
#' Retrieves the path to the methylCIPHER data directory. The function checks for a custom path set via
#' \code{\link{set_methylCIPHER_path}} or the \code{methylCIPHER.path} option
#' (see \code{\link{options}}). If no custom path is specified, it returns the
#' default path from \code{\link{get_methylCIPHER_default_path}}.
#'
#' Tips: To set a custom path permanently, add
#' \code{options(methylCIPHER.path = "non_default_path")} to your \code{.Rprofile}
#' file.
#'
#' @return A character string specifying the path to the methylCIPHER data directory.
#' @export
#'
#' @examples
#' # Returns the default path
#' get_methylCIPHER_path()
#' # Sets path to current directory
#' set_methylCIPHER_path(".")
#' # Returns the new path
#' get_methylCIPHER_path()
#' @seealso \code{\link{set_methylCIPHER_path}}, \code{\link{get_methylCIPHER_default_path}}
get_methylCIPHER_path <- function() {
  path <- getOption("methylCIPHER.path")
  if (is.null(path)) {
    return(get_methylCIPHER_default_path())
  }
  return(path)
}

#' Set Custom Path for methylCIPHER Data
#'
#' This function allows the user to set a custom directory path for storing
#' methylCIPHER data, which is required for SystemsAge and PCClocks. The specified
#' path is stored in R's options under "methylCIPHER.path".
#'
#' @param path A string specifying the custom path to the methylCIPHER data directory.
#'
#' @export
#'
#' @inherit get_methylCIPHER_path examples seealso
set_methylCIPHER_path <- function(path) {
  options(methylCIPHER.path = path)
}

#' Download methylCIPHER Clock Data
#'
#' This function downloads external data for specific methylCIPHER clocks to a specified directory.
#'
#' @param clocks A character vector specifying which clocks to download. Options
#' are `"all"`, `"SystemsAge"`, and `"PCClocks"`. Default is "all", which includes
#' both `"SystemsAge"` and `"PCClocks"`.
#' @param path The directory path where the downloaded files will be saved.
#' If \code{NULL}, it defaults to the path returned by [get_methylCIPHER_path()].
#' @param force A logical value indicating whether to force the download even if
#' the files already exist. Default is \code{FALSE}.
#'
#' @return Invisibly returns `TRUE` if downloads were attempted or if there was
#' nothing to download.
#'
#' @details Some clocks have large data dependencies that has to be first downloaded
#' before the clocks can be calculated. This function download the requested clock
#' to a specified path or default path returned by [get_methylCIPHER_path()]. SEe
#' [get_methylCIPHER_path()] on instructions on how to set/change the default path.
#'
#' @examples
#' \dontrun{
#' # Download all clocks to the default path
#' download_methylCIPHER()
#'
#' # Download only "SystemsAge" to a custom path
#' download_methylCIPHER(clocks = "SystemsAge", path = "/path/to/directory")
#'
#' # Force re-download of specified clock
#' download_methylCIPHER(clocks = "SystemsAge", force = TRUE)
#' }
#'
#' @export
download_methylCIPHER <- function(
    clocks = c("all", "SystemsAge", "PCClocks"),
    path = NULL,
    force = FALSE) {
  # Input validation
  clocks <- match.arg(clocks, several.ok = TRUE)
  if (is.null(path)) {
    path <- get_methylCIPHER_path()
    if(!dir.exists(path)) {
      dir.create(path, showWarnings = TRUE)
    }
  }
  checkmate::assert_directory_exists(path, access = "rw")
  # "all" = download all clocks
  if ("all" %in% clocks) {
    clocks <- c("SystemsAge", "PCClocks")
  }
  clocks <- unique(clocks)

  # Which files to download?
  download_tbl <- subset(large_clocks_data, clock %in% clocks)
  # Download to path
  download_tbl$download_to <- file.path(path, download_tbl$download_name)
  # Handling of force. If exists then don't download. If exist & force then re-download
  download_tbl$exists <- if (force) {
    FALSE
  } else {
    file.exists(download_tbl$download_to)
  }
  download_tbl <- subset(download_tbl, exists == FALSE)
  if (nrow(download_tbl) == 0) {
    message("Nothing to download. Set `force = TRUE` to re-download existing files.")
    return(invisible(TRUE))
  }
  # Authenticate with Google Drive if any files require it
  if (any(download_tbl$type == "googledrive")) {
    if (!requireNamespace("googledrive", quietly = TRUE)) {
      stop("Please install the 'googledrive' package to download these files (`install.packages('googledrive')`).")
    }
    googledrive::drive_auth()
    on.exit(googledrive::drive_deauth(), add = TRUE)
  }

  # Download
  message(paste("Attempting to download", nrow(download_tbl), "file(s)..."))
  for (i in seq_len(nrow(download_tbl))) {
    tryCatch(
      {
        if (download_tbl$type[i] == "googledrive") {
          googledrive::drive_download(
            file = googledrive::as_id(download_tbl$url[i]),
            path = download_tbl$download_to[i],
            overwrite = TRUE
          )
        } else {
          download.file(url = download_tbl$url[i], destfile = download_tbl$download_to[i])
        }
      },
      error = function(e) {
        warning(
          paste(
            "Failed to download:",
            download_tbl$download_name[i],
            "\nError:",
            e$message
          )
        )
      }
    )
  }
  return(invisible(TRUE))
}

#' Dependency Factory
#'
#' Some functions like [calcPCClocks()] and [calcSystemsAge()] requires large external
#' data to be downloaded and loaded into R. The best way to increase the speed
#' of this process is to use qs2 storage format. This function factory generate
#' functions to load in these data based on a template.
#'
#' @param object_name name of object ie SystemsAge_data
#' @param object_hash hash of object generated in data-raw
#'
#' @keywords internal
load_data_function <- function(object_name, object_hash) {
  function(path = NULL) {
    full_object_name <- paste0(force(object_name), ".qs2")
    # input validation
    checkmate::assert_character(path, len = 1, null.ok = TRUE)
    file_path <- if (!is.null(path)) {
      stopifnot("Provided file/dir doesn't exists" = dir.exists(path) || file.exists(path))
      info <- file.info(path)
      # if pass a folder
      if (info$isdir) {
        # then convert path to read into folder/file.qs2
        file.path(path, full_object_name)
      } else {
        path
      }
    } else {
      file.path(get_methylCIPHER_path(), full_object_name)
    }
    # read file with qs2
    tryCatch(
      {
        obj <- qs2::qs_read(file_path, validate_checksum = TRUE)
      },
      error = function(e) {
        stop(
          sprintf(
            paste(
              "Failed to read '%s'.",
              "Did you download '%s' and passed the path or loaded object to this function?",
              "See function `download_methylCIPHER` to download or re-download the required data.",
              "Function failed: %s",
              sep = "\n"
            ),
            full_object_name,
            full_object_name,
            e
          )
        )
      }
    )

    # check sum
    ## ensure that the object loaded with this function is intended.
    ## generated by data-raw/SystemsAge_data.R or similar scripts
    if ((rlang::hash(obj) != force(object_hash)) && object_hash != "skip") {
      stop(
        sprintf(
          "Data integrity check failed. Try re-downloading the '%s.qs2' file.",
          force(object_name)
        )
      )
    }
    return(obj)
  }
}

#' Load SystemsAge Data
#'
#' A wrapper that loads SystemsAge data from a specified or default location.
#'
#' @param path Character string specifying the file path to the .qs2 file.
#'   If `NULL`, loads from the path at [get_methylCIPHER_path()].
#' @return A list containing objects needed to calculate SystemsAge.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' SystemsAge_data <- load_SystemsAge_data()
#' SystemsAge_data <- load_SystemsAge_data("path/to/SystemsAge_data.qs2")
#' }
load_SystemsAge_data <- load_data_function("SystemsAge_data", "d984914ff6aa17d8a6047fed5f9f6e4d")

#' Load PCClocks Data
#'
#' A wrapper that loads PCClocks data from a specified or default location.
#'
#' @param path Character string specifying the file path to the .qs2 file.
#'   If `NULL`, loads from the path at [get_methylCIPHER_path()].
#' @return A list containing objects needed to calculate PCClocks.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' PCClocks_data <- load_PCClocks_data()
#' PCClocks_data <- load_PCClocks_data("path/to/PCClocks_data.qs2")
#' }
load_PCClocks_data <- load_data_function("PCClocks_data", "46386ec4be2b2a5239cf67b242d7dc24")
