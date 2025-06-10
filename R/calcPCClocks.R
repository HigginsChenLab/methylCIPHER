#' calcPCClocks
#'
#' @description A function to calculate PC Clocks
#'
#' @inheritParams param_template
#' @param pheno A data.frame containing the phenotype of the samples where each row must correspond to each row of `DNAm`. Must contains the following columns: `Female` = 1 or 0 and chronological age in years.
#' @param RData
#' Either a character string specifying the path to the methylCIPHER folder containing
#' the `CalcAllPCClocks.RData` file, or a list/environment containing the contents of
#' the `CalcAllPCClocks.RData` file.
#'
#' @return A data.frame of calculated clock values appened to `pheno`.
#' @export
#'
#' @examples
#' \dontrun{
#' example_PCClocks <- calcPCClocks(
#'   exampleBetas,
#'   examplePheno,
#'   # If `PCClocks_data.RData` has been downloaded to default folder at $HOME
#'   RData = get_methylCIPHER_path()
#' )
#' example_PCClocks
#' }
#'
calcPCClocks <- function(DNAm, pheno, RData = NULL) {
  if (!is.matrix(DNAm)) {
    # Temporary method before converting all other input into matrix
    DNAm <- as.matrix(DNAm)
  }
  # Input validation
  checkmate::assert(
    checkmate::check_null(RData),
    checkmate::check_character(RData, len = 1, any.missing = FALSE),
    checkmate::check_list(RData),
    checkmate::check_environment(RData),
    combine = "or"
  )
  checkmate::assert_data_frame(pheno, min.rows = 1, null.ok = FALSE)
  if(any(!c("Female", "Age") %in% names(pheno))) {
    stop("`pheno` must have a `Female` and `Age` column.")
  }
  checkmate::assert_integerish(pheno[["Female"]], lower = 0, upper = 1, null.ok = FALSE, any.missing = FALSE)
  checkmate::assert_numeric(pheno[["Age"]], finite = TRUE, null.ok = FALSE, any.missing = FALSE)
  warning("Currently, checks for rows in `pheno` corresponding to the same row as `DNAm` is not implemented")
  check_DNAm(DNAm)

  if (is.null(RData)) {
    if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path)) {
      current_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
    } else {
      current_path <- getwd() # fallback to current working directory
    }

    # add check that path_to_PCClocks ends with a /

    # path_to_PCClocks should end with a "/"
    # DNAm is a matrix of methylation Beta values, where row names are samples, and
    #   column names are CpGs
    # pheno has rows as samples and columns as phenotype variables. This can also
    #   include the original clocks if you used the Horvath online calculator as well.
    #   It MUST include a column named "Age" and a column named "Female"

    if (!("Age" %in% variable.names(pheno))) {
      stop("Error: pheno must have a column named Age")
    }
    if (!("Female" %in% variable.names(pheno))) {
      stop("Error: pheno must have a column named Female")
    }
    if (sum(startsWith(colnames(DNAm), "cg")) == 0) {
      warning("Warning: It looks like you may need to format DNAm using t(DNAm) to get samples as rows!")
    }

    # Note: this code assumes all your files are in one working directory. Alter the code as needed based on file locations.
    # Load packages
    pkgTest <- function(x) {
      if (!require(x, character.only = TRUE)) {
        install.packages(x, dep = TRUE)
        if (!require(x, character.only = TRUE)) stop("Package not found")
      }
    }
    pkgTest("dplyr")
    library(dplyr)
    pkgTest("tibble")
    library(tibble)
    pkgTest("tidyr")
    library(tidyr)
    pkgTest("googledrive")
    library(googledrive)

    googledrive::drive_auth()

    # In pheno, rows are samples and columns are phenotypic variables.
    # One of the phenotypic variables must be "Age", and another one "Female" (coded as Female = 1, Male = 0; should be a numeric variable as this will be included in PCGrimAge calculation)
    # Also ensure that the order of DNAm sample IDs matches your phenotype data sample IDs, otherwise your data will be scrambled

    home_dir <- Sys.getenv("HOME")

    if (file.exists(paste0(home_dir, "/CalcAllPCClocks.RData"))) {
      cat("PCClocks File Exists in Home Directory. Loading Data...\n")
      load(paste0(home_dir, "/CalcAllPCClocks.RData"))
    } else {
      cat("PCClocks Data does not exist in package directory. Downloading data...\n")
      public_file <- drive_get(as_id("1xhFUMBSrjRta3tgL0OTBLVgNlLnJJNRZ"))
      drive_download(public_file, path = paste0(home_dir, "/CalcAllPCClocks.RData"), overwrite = TRUE)
      load(paste0(home_dir, "/CalcAllPCClocks.RData"))
    }

    # load(file = paste(path_to_PCClocks_directory,"CalcAllPCClocks.RData", sep = ""))
  } else if (is.character(RData)) {
    if (!"CalcAllPCClocks.RData" %in% list.files(RData)) {
      stop("The file 'CalcAllPCClocks.RData' is not found in the provided path.")
    }
    cat("Loading Data...\n")
    PCClocks_env <- new.env()
    load(file = paste0(RData, "/", "CalcAllPCClocks.RData"), envir = PCClocks_env)
  } else if (is.list(RData) || is.environment(RData)) {
    PCClocks_env <- RData
  }
  message("PCClocks Data successfully loaded")

  attach(PCClocks_env)
  # Imputation
  datMeth <- as.data.frame(DNAm)
  if(length(c(CpGs[!(CpGs %in% colnames(datMeth))],CpGs[apply(datMeth[,colnames(datMeth) %in% CpGs], 2, function(x)all(is.na(x)))])) == 0){
    message("No CpGs were NA for all samples")
  } else{
    missingCpGs <- c(CpGs[!(CpGs %in% colnames(datMeth))])
    datMeth[,missingCpGs] <- NA
    datMeth = datMeth[,CpGs]
    missingCpGs <- CpGs[apply(datMeth[,CpGs], 2, function(x)all(is.na(x)))]
    for(i in 1:length(missingCpGs)){
      datMeth[,missingCpGs[i]] <- imputeMissingCpGs[missingCpGs[i]]
    }
    message("Any missing CpGs successfully filled in (see function for more details)")
  }

  #Prepare methylation data for calculation of PC Clocks (subset to 78,464 CpGs and perform imputation if needed)
  datMeth <- datMeth[,CpGs]
  meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
  datMeth <- apply(datMeth,2,meanimpute)
  #Note: you may substitute another imputation method of your choice (e.g. KNN), but we have not found the method makes a significant difference.
  message("Mean imputation successfully completed for any missing CpG values")

  #Initialize a data frame for PC clocks
  DNAmAge <- pheno

  # var = readline(prompt = "Type the column name in datPheno with sample names (or type skip):")
  # if(var != "skip"){
  #   if(sum(DNAmAge[, var] == row.names(datMeth)) != nrow(DNAmAge)){
  #     warning("Warning: It would appear that datPheno and datMeth do not have matching sample order! Check your inputs!")
  #   } else message("datPheno and datMeth sample order verified to match!")
  # }

  message("Calculating PC Clocks now")

  #Calculate PC Clocks
  DNAmAge$PCHorvath1 <- as.numeric(anti.trafo(sweep(as.matrix(datMeth),2,CalcPCHorvath1$center) %*% CalcPCHorvath1$rotation %*% CalcPCHorvath1$model + CalcPCHorvath1$intercept))
  DNAmAge$PCHorvath2 <- as.numeric(anti.trafo(sweep(as.matrix(datMeth),2,CalcPCHorvath2$center) %*% CalcPCHorvath2$rotation %*% CalcPCHorvath2$model + CalcPCHorvath2$intercept))
  DNAmAge$PCHannum <- as.numeric(sweep(as.matrix(datMeth),2,CalcPCHannum$center) %*% CalcPCHannum$rotation %*% CalcPCHannum$model + CalcPCHannum$intercept)
  DNAmAge$PCPhenoAge <- as.numeric(sweep(as.matrix(datMeth),2,CalcPCPhenoAge$center) %*% CalcPCPhenoAge$rotation %*% CalcPCPhenoAge$model + CalcPCPhenoAge$intercept)
  DNAmAge$PCDNAmTL <- as.numeric(sweep(as.matrix(datMeth),2,CalcPCDNAmTL$center) %*% CalcPCDNAmTL$rotation %*% CalcPCDNAmTL$model + CalcPCDNAmTL$intercept)
  temp <- cbind(sweep(as.matrix(datMeth),2,CalcPCGrimAge$center) %*% CalcPCGrimAge$rotation,Female = DNAmAge$Female,Age = DNAmAge$Age)
  DNAmAge$PCPACKYRS <- as.numeric(temp[,names(CalcPCGrimAge$PCPACKYRS.model)] %*% CalcPCGrimAge$PCPACKYRS.model + CalcPCGrimAge$PCPACKYRS.intercept)
  DNAmAge$PCADM <- as.numeric(temp[,names(CalcPCGrimAge$PCADM.model)] %*% CalcPCGrimAge$PCADM.model + CalcPCGrimAge$PCADM.intercept)
  DNAmAge$PCB2M <- as.numeric(temp[,names(CalcPCGrimAge$PCB2M.model)] %*% CalcPCGrimAge$PCB2M.model + CalcPCGrimAge$PCB2M.intercept)
  DNAmAge$PCCystatinC <- as.numeric(temp[,names(CalcPCGrimAge$PCCystatinC.model)] %*% CalcPCGrimAge$PCCystatinC.model + CalcPCGrimAge$PCCystatinC.intercept)
  DNAmAge$PCGDF15 <- as.numeric(temp[,names(CalcPCGrimAge$PCGDF15.model)] %*% CalcPCGrimAge$PCGDF15.model + CalcPCGrimAge$PCGDF15.intercept)
  DNAmAge$PCLeptin <- as.numeric(temp[,names(CalcPCGrimAge$PCLeptin.model)] %*% CalcPCGrimAge$PCLeptin.model + CalcPCGrimAge$PCLeptin.intercept)
  DNAmAge$PCPAI1 <- as.numeric(temp[,names(CalcPCGrimAge$PCPAI1.model)] %*% CalcPCGrimAge$PCPAI1.model + CalcPCGrimAge$PCPAI1.intercept)
  DNAmAge$PCTIMP1 <- as.numeric(temp[,names(CalcPCGrimAge$PCTIMP1.model)] %*% CalcPCGrimAge$PCTIMP1.model + CalcPCGrimAge$PCTIMP1.intercept)
  DNAmAge$PCGrimAge <- as.numeric(as.matrix(DNAmAge[,CalcPCGrimAge$components]) %*% CalcPCGrimAge$PCGrimAge.model + CalcPCGrimAge$PCGrimAge.intercept)

  message("PC Clocks successfully calculated!")
  detach(PCClocks_env)

  return(DNAmAge)
}
