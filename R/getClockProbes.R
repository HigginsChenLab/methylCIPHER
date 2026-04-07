#' getClockProbes
#'
#' @param DNAm The methylation Beta values you will calculate epigenetic clocks with, where columns are CpGs, and rows are samples.
#'
#' @return A table which compares the number of probes required by the full clocks in the current package, and the ones in your data. If you are missing many probes for a clock, imputation is necessary, though you might consider using a different clock or dataset if possible.
#' @export
#'
#' @examples getClockProbes(exampleBetas)
getClockProbes <- function(DNAm) {
  all_na_cols <- colSums(is.na(DNAm)) == nrow(DNAm)
  if (any(all_na_cols)) {
    warning("Warning! Some CpGs are NA in all samples. Suggest using removeNAprobes to remove these.")
  }

  ClockDataList <- as.data.frame(data(package = "methylCIPHER")[3][[1]])$Item
  ClockDataList <- grep("_CpG", ClockDataList, value = TRUE)
  ClockDataList <- ClockDataList[grep("Calc", ClockDataList, invert = TRUE)]

  results <- list()

  for (i in 1:length(ClockDataList)) {
    x <- get(ClockDataList[i])

    if (ClockDataList[i] %in% c(
      "Bocklandt_CpG", "DunedinPACE_CpGs", "DunedinPoAm38_CpGs", "EpiToc_CpGs", "hypoClock_CpGs",
      "Garagnani_CpG", "PCBrainAge_CpGs", "Weidner_CpGs", "PCClocks_CpGs", "SystemsAge_CpGs", "eightfiftykCpGs"
    )) {
      currentCpGList <- x
    } else if (ClockDataList[i] == "EpiToc2_CpGs") {
      currentCpGList <- rownames(x)
    } else if (ClockDataList[i] %in% c(
      "DNAmTL_CpGs", "HRSInCHPhenoAge_CpGs", "Horvath1_CpGs",
      "MiAge_CpGs", "PEDBE_CpGs", "Zhang2019_CpGs"
    )) {
      pull <- function(x1, y1) {
        x1[, if (is.name(substitute(y1))) deparse(substitute(y1)) else y1, drop = FALSE][[1]]
      }

      CpGColumn <- grepl("CpG|Marker|ID|id|name", colnames(x))
      currentCpGList <- pull(x, colnames(x)[CpGColumn])
    } else {
      CpGColumn <- grep("CpG|Marker|ID|id|name", colnames(x), value = TRUE)

      currentCpGList <- x[[CpGColumn[1]]]
    }

    clockName <- sub("_CpGs?$", "", ClockDataList[i])
    total <- length(currentCpGList)
    present <- sum(currentCpGList %in% colnames(DNAm))
    pct <- paste0(round((present / total) * 100, 0), "%")

    results[[length(results) + 1]] <- data.frame(
      Clock = clockName,
      Total.Probes = total,
      Present.Probes = present,
      Percent.Present = pct
    )
  }

  # Add PhysAge subscores from PhysAge_data
  for (nm in names(PhysAge_data)) {
    cpgs <- PhysAge_data[[nm]]$CpG
    cpgs <- cpgs[cpgs != "Intercept"]
    total <- length(cpgs)
    present <- sum(cpgs %in% colnames(DNAm))
    pct <- paste0(round((present / total) * 100, 0), "%")

    results[[length(results) + 1]] <- data.frame(
      Clock = paste0("PhysAge: ", nm),
      Total.Probes = total,
      Present.Probes = present,
      Percent.Present = pct
    )
  }

  # Add GrimAge V2 component subscores (V1 shared components have identical CpGs)
  grim2_components <- c("PACKYRS", "ADM", "B2M", "CystatinC", "GDF15", "Leptin", "PAI1", "TIMP1", "logA1C", "logCRP")
  for (comp in grim2_components) {
    cpgs <- names(CalcGrimAge2[[paste0(comp, ".model")]])
    cpgs <- cpgs[grepl("^cg|^ch", cpgs)]
    total <- length(cpgs)
    present <- sum(cpgs %in% colnames(DNAm))
    pct <- paste0(round((present / total) * 100, 0), "%")
    results[[length(results) + 1]] <- data.frame(
      Clock = paste0("GrimAgeV2: DNAm", comp),
      Total.Probes = total,
      Present.Probes = present,
      Percent.Present = pct
    )
  }

  # Add DNAmFitAge component subscores
  fitage_components <- list(
    "DNAmGait_noAge (F)" = "Gait_noAge_Females",
    "DNAmGait_noAge (M)" = "Gait_noAge_Males",
    "DNAmGrip_noAge (F)" = "Grip_noAge_Females",
    "DNAmGrip_noAge (M)" = "Grip_noAge_Males",
    "DNAmVO2max"         = "VO2maxModel",
    "DNAmGait_wAge (F)"  = "Gait_wAge_Females",
    "DNAmGait_wAge (M)"  = "Gait_wAge_Males",
    "DNAmGrip_wAge (F)"  = "Grip_wAge_Females",
    "DNAmGrip_wAge (M)"  = "Grip_wAge_Males",
    "DNAmFEV1_wAge (F)"  = "FEV1_wAge_Females",
    "DNAmFEV1_wAge (M)"  = "FEV1_wAge_Males"
  )
  for (label in names(fitage_components)) {
    model_name <- fitage_components[[label]]
    terms <- DNAmFitAge_data[[model_name]]$term
    cpgs <- terms[grepl("^cg|^ch", terms)]
    total <- length(cpgs)
    present <- sum(cpgs %in% colnames(DNAm))
    pct <- paste0(round((present / total) * 100, 0), "%")
    results[[length(results) + 1]] <- data.frame(
      Clock = paste0("DNAmFitAge: ", label),
      Total.Probes = total,
      Present.Probes = present,
      Percent.Present = pct
    )
  }

  # Add MethylCIPHERplus clocks if installed
  if (requireNamespace("MethylCIPHERplus", quietly = TRUE)) {
    plus_clocks <- list()

    # DNAmEMRAge: named weight vector, drop "Age" feature
    plus_clocks[["DNAmEMRAge"]] <- tryCatch({
      d <- getExportedValue("MethylCIPHERplus", "DNAmEMRAge_data")
      nms <- names(d$model)
      nms[grepl("^cg|^ch\\.", nms)]
    }, error = function(e) NULL)

    # OMICmAge: named weight vector, keep only CpG probes
    plus_clocks[["OMICmAge"]] <- tryCatch({
      d <- getExportedValue("MethylCIPHERplus", "OMICmAge_data")
      nms <- names(d$model)
      nms[grepl("^cg|^ch\\.", nms)]
    }, error = function(e) NULL)

    # FIAge: sparse matrix, row names minus intercept
    plus_clocks[["FIAge"]] <- tryCatch({
      d <- getExportedValue("MethylCIPHERplus", "FI_450K_Matrix")
      cpgs <- d@Dimnames[[1]]
      cpgs[cpgs != "(Intercept)"]
    }, error = function(e) NULL)

    # DNAmSAGE clocks: data frames with CpG column
    for (sage in c("DNAmSAGE_cGAS", "DNAmSAGE_p14Arf", "DNAmSAGE_p16Ink4a", "DNAmSAGE_p21")) {
      plus_clocks[[sage]] <- tryCatch({
        d <- getExportedValue("MethylCIPHERplus", paste0(sage, "_CpGs"))
        cpgs <- d$CpG
        cpgs[cpgs != "intercept"]
      }, error = function(e) NULL)
    }

    for (nm in names(plus_clocks)) {
      cpgs <- plus_clocks[[nm]]
      if (is.null(cpgs)) next
      total <- length(cpgs)
      present <- sum(cpgs %in% colnames(DNAm))
      pct <- paste0(round((present / total) * 100, 0), "%")

      results[[length(results) + 1]] <- data.frame(
        Clock = paste0("MCP+ ", nm),
        Total.Probes = total,
        Present.Probes = present,
        Percent.Present = pct
      )
    }
  }

  ProbeTable <- do.call(rbind, results)
  rownames(ProbeTable) <- NULL
  ProbeTable
}
