
ShuffleCells <- function(fcs.files) {
  x <- c()
  for (i in 1:length(fcs.files)) {
    subsamp <- nrow(fcs.files[[i]])
    df1 <- fcs.files[[i]]
    x <- c(x, nrow(df1))
    df2 <- df1[sample(nrow(df1)), ]
    fcs.files[[i]] <- df2
    rownames(fcs.files[[i]]) <- seq(1:subsamp)
  }
  return(fcs.files)
}

MultiShuffleCells <- function(fcs.files) {
  x <- c()
  for (n in 1:length(fcs.files)) {
    for (i in 1:length(fcs.files[[n]])) {
      subsamp <- nrow(fcs.files[[n]][[i]])
      df1 <- fcs.files[[n]][[i]]
      x <- c(x, nrow(df1))
      df2 <- df1[sample(nrow(df1)), ]
      fcs.files[[n]][[i]] <- df2
      rownames(fcs.files[[n]][[i]]) <- seq(1:subsamp)
    }
  }
  return(fcs.files)
}

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
  alt.times <- names(files)
  if (name.sort) {
    alt.times <- alt.times[order(as.numeric(alt.times))]
  }
  if (!identical(alt.times, times)) {
    warning("Times from list names do not match times from FCS file names! Using times from FCS file names.")
  }
  return(times)
}
