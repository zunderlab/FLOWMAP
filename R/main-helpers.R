
CheckModeOne <- function(files) {
  # "files" variable could be one of the following:
  # a single fcs file path
  # a single folder path containing 1 fcs files
  # a vector of 1 fcs file path
  fail.flag <- TRUE
  guide <- NULL
  if (length(files) == 1) {
    if (grepl(pattern = "\\.fcs", files)) {
      fail.flag <- FALSE
      guide <- "FCS"
    } else if (length(list.files(files)) == 1) {
      if (grepl(pattern = "\\.fcs", list.files(files))) {
        fail.flag <- FALSE
        guide <- "folder"
      }
    }
  }
  return(c(fail.flag, guide))
}

CheckModeSingle <- function(files) {
  fail.flag <- TRUE
  guide <- NULL
  # "files" variable could be one of the following:
  # a single folder path containing 2+ fcs files
  # a vector of fcs file paths
  if (length(files) == 1) {
    if (length(list.files(files)) > 1) {
      if (sum(grepl(pattern = "\\.fcs", list.files(files))) > 1) {
        # if (all(grepl(pattern = "\\.fcs", list.files(files)))) {
        fail.flag <- FALSE
        guide <- "folder"
      }
    }
  } else if (length(files) > 1) {
    if (all(grepl(pattern = "\\.fcs", files))) {
      fail.flag <- FALSE
      guide <- "FCS"
    }
  }
  return(c(fail.flag, guide))
}

CheckModeMulti <- function(files) {
  fail.flag <- TRUE
  guide <- NULL
  # "files" variable could be one of the following:
  # a single folder path, containing subfolders,
  # which each contain 2+ fcs files
  # a list name by time, each element is a vector
  # of 1+ fcs file paths from each treatment/condition
  if (length(files) == 1) {
    if (length(list.files(files)) > 1) {
      subfolders <- list.files(files, full.names = TRUE)
      temp.fail.flag <- FALSE
      for (n in subfolders) {
        if (!(sum(grepl(pattern = "\\.fcs", list.files(n))) >= 1)) {
          temp.fail.flag <- TRUE
        }
      }
      fail.flag <- temp.fail.flag
      guide <- "subfolder"
    }
  } else if (length(files) > 1) {
    multiple.condition.flag <- FALSE
    temp.fail.flag <- FALSE
    if (is.list(files)) {
      for (i in 1:length(files)) {
        if (!all(grepl(pattern = "\\.fcs", files[[i]]))) {
          temp.fail.flag <- TRUE
        }
        if (length(files[[i]]) >= 1) {
          multiple.condition.flag <- TRUE
        }
      }
    }
    if (multiple.condition.flag & !temp.fail.flag) {
      fail.flag <- FALSE
      guide <- "list"
    }
  }
  return(c(fail.flag, guide))
}

ParseTimes <- function(fcs.file.names, name.sort) {
  times <- c()
  for (i in 1:length(fcs.file.names)) {
    this.name <- fcs.file.names[i]
    this.name <- basename(this.name)
    this.name <- gsub("\\.fcs", "", this.name)
    this.name <- unlist(strsplit(this.name, ""))
    this.name <- this.name[suppressWarnings(!is.na(as.numeric(this.name)))]
    if (length(this.name) > 1) {
      this.name <- paste(this.name, collapse = "")
    }
    times <- c(times, this.name)
    rm(this.name)
  }
  if (name.sort) {
    times <- times[order(as.numeric(times))]
  }
  print(paste0("Timepoints: ", times))
  return(times)
}

MultiFolderParseTimes <- function(files, fcs.file.names, name.sort) {
  times <- c()
  for (i in 1:length(fcs.file.names)) {
    if (length(fcs.file.names[[i]]) > 1) {
      warning("Multiple FCS files per timepoint, only using the first for time label!")
      times <- c(times, ParseTimes(fcs.file.names[[i]][1], name.sort))
    } else {
      times <- c(times, ParseTimes(fcs.file.names[[i]], name.sort))
    }
  }
  alt.times <- ParseTimes(list.files(files), name.sort)
  if (!identical(alt.times, times)) {
    warning("Times from subfolder names do not match times from FCS file names! Using times from FCS file names.")
  }
  return(times)
}

MultiListParseTimes <- function(fcs.file.names, name.sort) {
  times <- c()
  for (i in 1:length(fcs.file.names)) {
    times <- c(times, ParseTimes(fcs.file.names[[i]], name.sort))
  }
  times.temp <- unique(times)
  alt.times <- names(fcs.file.names)
  if (name.sort) {
    alt.times <- alt.times[order(as.numeric(alt.times))]
  }
  if (!identical(alt.times, times)) {
    warning("Times from list names do not match times from FCS file names! Using times from subfolder names.")
  }
  return(alt.times)
}

ProcessConditions <- function(list.of.clean.FCS.files, fcs.file.names) {
  label.key.special <- c()
  for (i in 1:length(list.of.clean.FCS.files)) {
    this.name <- basename(fcs.file.names[i])
    this.name <- gsub(".fcs", "", this.name)
    if (grepl("-", this.name)) {
      this.name <- unlist(strsplit(this.name, split = "-"))[1]
    }
    if (grepl("\\.", this.name)) {
      this.name <- unlist(strsplit(this.name, split = "\\."))[1]
    }
    label.key.special <- c(label.key.special, this.name)
    list.of.clean.FCS.files[[i]] <- cbind(list.of.clean.FCS.files[[i]],
                                          Condition = rep(i, times = nrow(list.of.clean.FCS.files[[i]])))
  }
  results <- list(fixed.files = list.of.clean.FCS.files,
                  label.key.special = label.key.special)
  return(results)
}

#' ConstructVarAnnotate
#'
#' \code{ConstructVarAnnotate} returns a list that can be used during FLOWMAPR
#' analysis to change channel names from the default "name" field used in a
#' flowFrame to the names provided in the "desc" fields. Generally in CyTOF/mass
#' cytometry or flow cytometry data, the original "name" field corresponds to the
#' measurement channel (e.g. FITC or Sm154Di) while the "desc" field will correspond
#' to the name of the marker entered by the user at the machine (e.g. cCaspase3 or
#' CyclinB1).
#'
#' @param FCS.file.name The full file path to one FCS file intended for analysis
#' @return List mapping channel names (taken from the "name" field of a
#' flowFrame) to alternative names from the flowFrame "desc" fields
#' @examples
#' \dontrun{ConstructVarAnnotate(FCS.file.name = system.file("extdata/SingleFLOWMAP/d1.fcs",package = "FLOWMAPR"))}
#' @export
ConstructVarAnnotate <- function(FCS.file.name) {
  fcs.file <- read.FCS(FCS.file.name)
  fcs.file.matrix <- exprs(fcs.file)
  channel.names <- c()
  marker.names <- c()
  var.annotate <- list()
  for (i in 1:length(names(colnames(fcs.file.matrix)))) {
    this.name <- names(colnames(fcs.file.matrix))[i]
    channel.names <- c(channel.names, description(fcs.file)[[this.name]])
    this.name <- gsub(this.name, pattern = "N", replacement = "S")
    marker.names <- c(marker.names, description(fcs.file)[[this.name]])
    rm(this.name)
  }
  for (i in 1:length(marker.names)) {
    if (marker.names[i] == " ") {
      marker.names[i] <- channel.names[i]
    }
    temp <- gsub(marker.names[i], pattern = "-", replacement = "_")
    temp <- gsub(temp, pattern = " ", replacement = "_")
    temp <- gsub(temp, pattern = "/", replacement = "_")
    var.annotate[[channel.names[i]]] <- temp
  }
  return(var.annotate)
}

#' SuggestClusteringVar
#'
#' \code{SuggestClusteringVar} returns a vector of user-specified length that contains
#' suggested names of channels to be used as \code{clustering.var} in FLOWMAPR analysis.
#' Chosen channels are determined based on analysis of the provided FCS files, as the
#' variance is calculated within and between timepoints. The top varying markers are noted
#' and tabulated, and the markers that consistently vary the most are selected.
#'
#' @param fcs.file.names A vector of full file paths to the FCS files to be used
#' @param mode FLOWMAPR mode to use in analysis based on starting input,
#' available options include \code{c("single", "multi", "one")}
#' @param var.annotate List mapping channel names to user-specified names to properly
#' annotate all FCS file data
#' @param var.remove Vector naming channels to be removed from all loaded FCS data
#' @param top.num Numeric specifying the number of variables in the vector to be returned
#' @return Vector naming channels that vary the most within and between FCS files
#' @examples
#' fcs.file.names <- c("Desktop/1.fcs", "Desktop/2.fcs")
#' var.annotate <- ConstructVarAnnotate[1]
#' var.remove <- c("Channel3", "Channel4")
#'
#' SuggestClusteringVar(fcs.file.names, mode = "single", var.annotate,
#' var.remove, top.num = 20)
#' @export
SuggestClusteringVar <- function(fcs.file.names, mode, var.annotate,
                                 var.remove, top.num) {
  suggested.clustering.var <- c()
  combined.fcs.files <- c()
  if (mode == "one") {
    fcs.file <- LoadCleanFCS(fcs.file.names = fcs.file.names, channel.remove = var.remove,
                             channel.annotate = var.annotate, subsamples = FALSE)
    combined.fcs.files <- fcs.file[[1]]
  } else if (mode == "single") {
    fcs.files <- LoadCleanFCS(fcs.file.names = fcs.file.names, channel.remove = var.remove,
                              channel.annotate = var.annotate, subsamples = FALSE)
    for (i in 1:length(fcs.files)) {
      combined.fcs.files <- rbind(combined.fcs.files, fcs.files[[i]])
    }
  } else if (mode == "multi") {
    fcs.files <- LoadMultiCleanFCS(list.of.file.names = fcs.file.names, channel.remove = var.remove,
                                   channel.annotate = var.annotate, subsamples = FALSE)
    for (i in 1:length(fcs.files)) {
      for (j in 1:length(fcs.files[[i]])) {
        combined.fcs.files <- rbind(combined.fcs.files, fcs.files[[i]][[j]])
      }
    }
  } else {
    stop("User-specified mode not recognized!")
  }

  if (ncol(combined.fcs.files) < top.num) {
    stop("Requesting more suggested clustering var than available in data!")
  }

  all.var <- apply(combined.fcs.files, 2, var)
  if (mode == "one") {
    top.selected.var <- sort(all.var, decreasing = TRUE)[1:top.num]
  } else if (mode == "single") {
    # get variance within each timepoint
    var.over.time <- c()
    for (i in 1:length(fcs.files)) {
      temp.short.df <- fcs.files[[i]][, 2:ncol(fcs.files[[i]])]
      this.var <- apply(temp.short.df, 2, var)
      var.over.time <- rbind(var.over.time, this.var)
    }
    # get variance within adjacent timepoints
    var.cross.time <- c()
    for (i in 1:(length(fcs.files) - 1)) {
      temp.short.df1 <- fcs.files[[i]][, 2:ncol(fcs.files[[i]])]
      temp.short.df2 <- fcs.files[[(i + 1)]][, 2:ncol(fcs.files[[i]])]
      temp.short.df <- rbind(temp.short.df1, temp.short.df2)
      cross.var <- apply(temp.short.df, 2, var)
      var.cross.time <- rbind(var.cross.time, cross.var)
    }
    all.vars.time <- rbind(var.over.time[1, ], var.cross.time[1, ], var.over.time[2, ],
                           var.cross.time[2, ], var.over.time[3, ], var.cross.time[3, ],
                           var.over.time[4, ])
    median.var.time <- apply(all.vars.time, 2, median)
    # take top most varying markers
    top.selected.var <- sort(median.var.time, decreasing = TRUE)[1:top.num]
  } else if (mode == "multi") {
    var.over.time <- c()
    # get variance within each timepoint
    for (i in 1:length(fcs.files)) {
      temp.combined <- c()
      for (j in 1:length(fcs.files[[i]])) {
        temp.combined <- rbind(temp.combined, fcs.files[[i]][[j]])
      }
      this.var <- apply(temp.combined, 2, var)
      var.over.time <- rbind(var.over.time, this.var)
    }
    # get variance within adjacent timepoints
    var.cross.time <- c()
    for (i in 1:(length(fcs.files) - 1)) {
      temp.combined.1 <- c()
      temp.combined.2 <- c()
      for (j in 1:length(fcs.files[[i]])) {
        temp.combined.1 <- rbind(temp.combined.1, fcs.files[[i]][[j]])
      }
      for (j in 1:length(fcs.files[[i + 1]])) {
        temp.combined.2 <- rbind(temp.combined.2, fcs.files[[i + 1]][[j]])
      }
      temp.combined <- rbind(temp.combined.1, temp.combined.2)
      cross.var <- apply(temp.combined, 2, var)
      var.cross.time <- rbind(var.cross.time, cross.var)
    }

    all.vars.time <- c()
    for (i in 1:(nrow(var.over.time) - 1)) {
      all.vars.time <- rbind(all.vars.time, var.over.time[i, ])
      all.vars.time <- rbind(all.vars.time, var.cross.time[i, ])
    }
    all.vars.time <- rbind(all.vars.time, var.over.time[(nrow(var.over.time)), ])
    median.var.time <- apply(all.vars.time, 2, median)
    # take top most varying markers
    top.selected.var <- sort(median.var.time, decreasing = TRUE)[1:top.num]
  }
  suggested.clustering.var <- names(top.selected.var)
  return(suggested.clustering.var)
}

#' SuggestVarRemove
#'
#' \code{SuggestVarRemove} returns a vector of channel names that suggest markers
#' to remove during FLOWMAPR analysis based on the annotated channel list \code{var.annotate}
#' and \code{var.to.remove} a user-supplied list of substrings that designate
#' the names of blank channels.
#'
#' @param var.annotate List mapping channel names to user-specified names to properly
#' annotate all FCS file data
#' @param var.to.remove Vector of substrings to, default is set to \code{NULL}
#' @return Vector naming channels suggested to be removed from downstream FLOWMAPR analysis
#' @examples
#' var.annotate <- ConstructVarAnnotate("Desktop/A.fcs")
#' SuggestVarRemove(var.annotate, var.to.remove = "Blank_Channel")
#' @export
SuggestVarRemove <- function(var.annotate, var.to.remove = NULL) {
  suggested.var.remove <- c()
  final.var.names <- unname(unlist(var.annotate))
  if (is.null(var.to.remove)) {
    usual.var.remove <- c("bead", "DNA", "BC", "Event", "length", "Time")
    channel.blank <- c("Dd", "Di")
    var.to.remove <- c(usual.var.remove, channel.blank)
  }
  for (i in var.to.remove) {
    suggested.var.remove <- c(suggested.var.remove,
                              final.var.names[grepl(pattern = i,
                                                    x = final.var.names)])
  }
  return(suggested.var.remove)
}
