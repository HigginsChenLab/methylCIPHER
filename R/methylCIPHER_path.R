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
#' [set_methylCIPHER_path()] or the `methylCIPHER.path` option
#' (see [options()]). If no custom path is specified, it returns the
#' default path from [get_methylCIPHER_default_path()].
#'
#' Tips: To set a custom path permanently, add
#' `options(methylCIPHER.path = "non_default_path")` to your `.Rprofile`
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
#' @seealso [set_methylCIPHER_path()], [get_methylCIPHER_default_path()]
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
#' This function downloads external data for specific methylCIPHER clocks to a specified
#' directory.
#'
#' @param clocks A character vector specifying which clocks to download. Valid options
#'   include `"all"`, `"SystemsAge"`, and `"PCClocks"`. The default is `"all"`, which
#'   includes both `"SystemsAge"` and `"PCClocks"`. Multiple specific clocks can be
#'   requested, e.g., `c("SystemsAge", "PCClocks")`.
#' @param path The directory path where the downloaded files will be saved.
#'   If `NULL` (the default), it uses the path returned by [get_methylCIPHER_path()].
#' @param source A character string specifying the download source: `"googledrive"` or
#'   `"zenodo"`. The default is `"googledrive"`. Note that `"zenodo"` can sometimes be slow
#'   and time out, while `"googledrive"` requires authentication via web browser.
#' @param force A logical value indicating whether to force the download even if
#'   the files already exist in the target directory. Defaults to `FALSE`.
#' @param ... Additional arguments passed to [googledrive::drive_auth()] or
#'   [zen4R::ZenodoRecord()].
#'
#' @return Invisibly returns `TRUE` if downloads were attempted or if there was
#'   nothing to download.
#'
#' @details Some clocks have large data dependencies that must be downloaded before the
#'   clocks can be calculated. This function downloads the requested clock data to a
#'   specified path or the default path returned by [get_methylCIPHER_path()].
#'   See [get_methylCIPHER_path()] for instructions on how to set or change
#'   the default path.
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
#'
#' # Download from Google Drive (may prompt for authentication)
#' download_methylCIPHER(source = "googledrive")
#' }
#'
#' @export
download_methylCIPHER <- function(
    clocks = c("all", "SystemsAge", "PCClocks", "test"),
    source = c("googledrive", "zenodo"),
    path = NULL,
    force = FALSE,
    ...) {
  # pre-conditioning
  clocks <- match.arg(clocks, several.ok = TRUE)
  source <- match.arg(source)

  if (is.null(path)) {
    path <- get_methylCIPHER_path()
    if (!dir.exists(path)) {
      message("Creating a new folder at ", path)
      dir.create(path, showWarnings = TRUE, recursive = TRUE)
    }
  }

  # Validate path exists and is writable
  checkmate::assert_directory_exists(path, access = "rw", .var.name = "path")

  # "all" = download all clocks
  clocks <- unique(clocks)
  if ("test" %in% clocks) {
    clocks <- "methylCIPHER_test.csv"
    download_name <- clocks
  } else if ("all" %in% clocks) {
    clocks <- c("PCClocks", "SystemsAge")
    download_name <- paste0(clocks, "_data.qs2")
  }
  download_to <- file.path(path, download_name)

  exists <- if (force) {
    rep(FALSE, times = length(download_to))
  } else {
    file.exists(download_to)
  }

  if (any(exists)) {
    message(paste("Skipping", sum(exists), "file(s) that already exist. Use force = TRUE to re-download."))
    clocks <- clocks[!exists]
    download_name <- download_name[!exists]
    download_to <- download_to[!exists]
  }

  if (length(clocks) == 0) {
    message("No files to download.")
    return(invisible(TRUE))
  }

  if (source == "googledrive") {
    if (!requireNamespace("googledrive", quietly = TRUE)) {
      stop("Please install the 'googledrive' package to download these files (`install.packages('googledrive')`).")
    }
    googledrive::drive_auth()
    on.exit(googledrive::drive_deauth(), add = TRUE)
  }

  if (source == "zenodo") {
    if (!requireNamespace("zen4R", quietly = TRUE)) {
      stop("Please install the 'zen4R' package to download these files (`install.packages('zen4R')`).")
    }
    zenodo <- zen4R::ZenodoManager$new(logger = "INFO")
    rec <- zenodo$getRecordByDOI("10.5281/zenodo.17162604")
  }

  message(paste("Attempting to download", length(clocks), "file(s)..."))
  download_success <- logical(length(clocks))
  names(download_success) <- download_name

  for (i in seq_along(clocks)) {
    tryCatch(
      {
        message("Downloading ", clocks[i], "...")

        if (source == "googledrive") {
          file_i <- subset(large_clocks_data$googledrive, path %in% download_name[i])
          file_i <- as.list(file_i)
          googledrive::drive_download(
            file = file_i$id,
            path = download_to[i],
            overwrite = TRUE,
            ...
          )
          download_success[i] <- TRUE
        }

        if (source == "zenodo") {
          rec$downloadFiles(path = path, files = list(download_name[i]), ...)
          download_success[i] <- TRUE
        }

        message("Successfully downloaded ", clocks[i])
      },
      error = function(e) {
        message("Error: ", e$message)
        warning(paste("Failed to download:", download_name[i]), call. = FALSE)
        download_success[i] <- FALSE
      }
    )
  }

  return(invisible(download_success))
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
