
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
  alt.times <- names(fcs.file.names)
  if (name.sort) {
    alt.times <- alt.times[order(as.numeric(alt.times))]
  }
  if (!identical(alt.times, times)) {
    warning("Times from list names do not match times from FCS file names! Using times from FCS file names.")
  }
  return(times)
}

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
  print("channel.names")
  print(channel.names)
  print("marker.names")
  print(marker.names)
  for (i in 1:length(marker.names)) {
    if (marker.names[i] == " ") {
      marker.names[i] <- channel.names[i]
    }
    var.annotate[[channel.names[i]]] <- marker.names[i]
  }
  print("var.annotate")
  print(var.annotate)
  return(var.annotate)
}

SuggestClusteringVar <- function(fcs.file.names, var.annotate, var.remove) {
  
  suggested.clustering.var <- c()
  
  return(suggested.clustering.var)
}

SuggestVarRemove <- function(var.annotate) {
  usual.var.remove <- c("bead", "DNA", "BC", "Event", "length", "Time")
  channel.blank <- c("Dd", "Di")
  
  suggested.var.remove <- c()
  
  return(suggested.var.remove)
}
