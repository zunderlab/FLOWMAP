AsinhNik <- function(value) {
  # Define a function for arcsinh transformation (based on definition from Nikolay).
  value <- value - 1
  for(i in 1:length(value)) {
    if ((value[i] < 0) | is.na(value[i])) {
      value[i] <- rnorm(1, mean = 0, sd = 0.01)
    }
    value <- value / 5
    value <- asinh(value) # value <- log(value + sqrt(value^2 + 1))  
    return(value)
  }
}


GetFCSNames <- function(folder, file.format, sort = TRUE) {
  # get FCS files
  fcs.files = list.files(path = folder, pattern = file.format,
                         recursive = FALSE, full.names = TRUE)
  if (sort) {
    # sort to organize by hour
    fcs.files <- sort(fcs.files)
  }
  return(fcs.files)
}


GetMultiFCSNames <- function(folder, file.format, sort = TRUE) {
  # get FCS files
  subfolders <- list.files(folder)
  list.of.treat.file.names <- list()
  for (treat in subfolders) {
    folder.name <- treat
    setwd(treat)
    x <- getwd()
    list.of.treat.file.names[[treat]] <- GetFCSNames(x, file.format, sort = TRUE)
    setwd(folder)
  }
  return(list.of.treat.file.names)
}


LoadCleanFCS <- function(fcs.file.names, channel.remove, channel.annotate,
                         subsample = 10000, subsample.rand = FALSE,
                         transform = TRUE, scale = FALSE) {
  clean.fcs.files <- list()
  for (i in 1:length(fcs.file.names)) {
    # print FCS file that is being read
    current.file <- tail(strsplit(fcs.file.names[i],"/")[[1]],n=1)
    cat("Reading FCS file data from:", current.file, "\n")
    # store currently read FCS file
    if (!subsample) {
      tmp.FCS1 <- read.FCS(fcs.file.names[i], transformation = "linearize", which.lines = NULL,
                           alter.names = FALSE, column.pattern = NULL, invert.pattern = FALSE,
                           decades = 0, ncdf = FALSE, min.limit = NULL, dataset = NULL,
                           emptyValue = TRUE)  
      tmp.FCS1 <- head(tmp.FCS1, n = unname(dim(tmp.FCS1)[1]))
    }
    else { 
      cat("Subsampling", current.file, "to", subsample, "cells\n")
      if (subsample.rand) {
        subsamp.FCS1 <- read.FCS(fcs.file.names[i], transformation = "linearize", which.lines = subsample,
                                 alter.names = FALSE, column.pattern = NULL,
                                 invert.pattern = FALSE, decades = 0, ncdf = FALSE,
                                 min.limit = NULL, dataset = NULL, emptyValue = TRUE)
        tmp.FCS1 <- head(subsamp.FCS1, n = subsample)
      }
      else {
        full.FCS1 <- read.FCS(fcs.file.names[i], transformation = "linearize", which.lines = NULL,
                              alter.names = FALSE, column.pattern = NULL,
                              invert.pattern = FALSE, decades = 0, ncdf = FALSE,
                              min.limit = NULL, dataset = NULL, emptyValue = TRUE)
        tmp.FCS1 <- head(full.FCS1, n = subsample)
      }
    }
    # rename variables with protein marker measured instead of metal channel
    cat("Fixing channel names from:", current.file, "\n")
    for (x in 1:length(colnames(tmp.FCS1))) {
      if (exists(colnames(tmp.FCS1)[x], where = channel.annotate)) {
        colnames(tmp.FCS1)[x] <- channel.annotate[[colnames(tmp.FCS1)[x]]]
      }
    }
    # remove unneeded variables
    cat("Removing unnecessary channel names from:", current.file, "\n")
    tmp.FCS2 <- subset(tmp.FCS1, select = colnames(tmp.FCS1)[!colnames(tmp.FCS1) %in% channel.remove])
    if (transform) {
      cat("Transforming data from:",current.file,"\n")
      tmp.FCS3 <- apply(tmp.FCS2, 2, AsinhNik) 
    }
    if (scale) {
      print("no scale implemented yet")
    }
    tmp.FCS4 <- as.data.frame(tmp.FCS3)
    clean.fcs.files[[i]] <- tmp.FCS4
    rm(tmp.FCS1, tmp.FCS2, tmp.FCS3, tmp.FCS4)
  }  
  return(clean.fcs.files)
}


LoadMultiCleanFCS <- function(list.of.treat.file.names, channel.remove, channel.annotate,
                              subsample = 10000, subsample.rand = FALSE, transform = TRUE, scale = FALSE) {
  list.of.treat.clean.FCS.files <- list()
  for (treat in names(list.of.treat.file.names)) {
    cat("treat is", treat, "\n")
    fcs.file.names <- list.of.treat.file.names[[treat]]
    cat("fcs.file.names are", fcs.file.names, "\n")
    list.of.treat.clean.FCS.files[[treat]] <- LoadCleanFCS(fcs.file.names,
                                                           channel.remove,
                                                           channel.annotate,
                                                           subsample,
                                                           subsample.rand = FALSE,
                                                           transform, scale)
    for (i in 1:length(list.of.treat.clean.FCS.files[[treat]])) {
      Treat <- rep(treat, times = dim(list.of.treat.clean.FCS.files[[treat]][[i]])[1])
      list.of.treat.clean.FCS.files[[treat]][[i]] <- cbind(list.of.treat.clean.FCS.files[[treat]][[i]],
                                                           Treat)
      rm(Treat)
    }
  }
  return(list.of.treat.clean.FCS.files)
}

ConvertNumericLabel <- function(list.of.clean.FCS.files.with.labels) {
  label.key <- list()
  fixed.list.FCS.files <- list()
  for (i in 1:length(list.of.clean.FCS.files.with.labels)) {
    label.key[[i]] <- names(list.of.clean.FCS.files.with.labels)[i]
    tmp.label.fcs <- list.of.clean.FCS.files.with.labels[[i]]
    for (n in 1:length(tmp.label.fcs)) {
      inds <- which(tmp.label.fcs[[n]] == label.key[[i]], arr.ind = TRUE)
      col.ind <- unname(inds[1, 2])
      tmp.label.fcs[[n]][, col.ind] <- i
    }
    fixed.list.FCS.files[[i]] <- tmp.label.fcs
  }
  names(fixed.list.FCS.files) <- names(list.of.clean.FCS.files.with.labels)
  return(list(fixed.list.FCS.files = fixed.list.FCS.files,
              label.key = label.key))
}
