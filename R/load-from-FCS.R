
Asinh <- function(value) {
  # Define a function for arcsinh transformation (based on definition from Nikolay).
  value <- value - 1
  for(i in 1:length(value)) {
    if ((value[i] < 0) | is.na(value[i])) {
      value[i] <- rnorm(1, mean = 0, sd = 0.01)
    }
  }
  value <- value / 5
  value <- asinh(value) # value <- log(value + sqrt(value^2 + 1))  
  return(value)
}

GetFCSNames <- function(folder, sort = TRUE) {
  # get FCS files
  fcs.files <- list.files(path = folder, pattern = "\\.fcs",
                          recursive = FALSE, full.names = TRUE)
  if (sort) {
    # sort to organize by hour
    fcs.files <- sort(fcs.files)
  }
  return(fcs.files)
}

GetMultiFCSNames <- function(folder, sort = TRUE) {
  # get FCS files
  subfolders <- list.files(folder, full.names = TRUE)
  list.of.time.file.names <- list()
  for (folder in subfolders) {
    fcs.files <- list.files(path = folder, pattern = "\\.fcs",
                            recursive = FALSE, full.names = TRUE)
    time <- basename(folder)
    if (sort) {
      # sort to organize by hour
      fcs.files <- sort(fcs.files)
    }
    list.of.time.file.names[[time]] <- fcs.files
  }
  return(list.of.time.file.names)
}

LoadCleanFCS <- function(fcs.file.names, channel.remove, channel.annotate,
                         subsamples = 1000, subsample.rand = TRUE,
                         transform = TRUE, scale = FALSE) {
  clean.fcs.files <- list()
  print("subsamples")
  print(subsamples)
  print("T")
  print(T)
  print("class(T)")
  print(class(T))
  
  if (length(subsamples) == 1 & subsamples != FALSE) {
    cat("Subsampling all files to:", subsamples, "\n")
    subsample.new <- rep(subsamples, times = length(fcs.file.names))
    subsamples <- subsample.new
  }
  for (i in 1:length(fcs.file.names)) {
    current.file <- tail(strsplit(fcs.file.names[i], "/")[[1]], n = 1)
    cat("Reading FCS file data from:", current.file, "\n")
    # store currently read FCS file
    print("a")
    print("!subsamples")
    print(!subsamples)
    if (!subsamples) {
      print("b")
      tmp.FCS1 <- read.FCS(fcs.file.names[i], transformation = "linearize", which.lines = NULL,
                           alter.names = FALSE, column.pattern = NULL, invert.pattern = FALSE,
                           decades = 0, ncdf = FALSE, min.limit = NULL, dataset = NULL,
                           emptyValue = TRUE)  
      tmp.FCS1 <- head(tmp.FCS1, n = unname(dim(tmp.FCS1)[1]))
    } else { 
      cat("Subsampling", current.file, "to", subsamples[i], "cells\n")
      if (subsample.rand) {
        subsamp.FCS1 <- read.FCS(fcs.file.names[i], transformation = "linearize",
                                 which.lines = subsamples[i],
                                 alter.names = FALSE, column.pattern = NULL,
                                 invert.pattern = FALSE, decades = 0, ncdf = FALSE,
                                 min.limit = NULL, dataset = NULL, emptyValue = TRUE)
        tmp.FCS1 <- head(subsamp.FCS1, n = subsamples[i])
      }
      else {
        full.FCS1 <- read.FCS(fcs.file.names[i], transformation = "linearize", which.lines = NULL,
                              alter.names = FALSE, column.pattern = NULL,
                              invert.pattern = FALSE, decades = 0, ncdf = FALSE,
                              min.limit = NULL, dataset = NULL, emptyValue = TRUE)
        tmp.FCS1 <- head(full.FCS1, n = subsamples[i])
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
      tmp.FCS3 <- apply(tmp.FCS2, 2, Asinh) 
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

LoadMultiCleanFCS <- function(list.of.time.file.names, channel.remove, channel.annotate,
                              subsamples = 1000, subsample.rand = TRUE, transform = TRUE, scale = FALSE) {
  list.of.time.clean.FCS.files <- list()
  subsamp.orig <- subsamples
  for (time in names(list.of.time.file.names)) {
    # cat("time is", time, "\n")
    fcs.file.names <- list.of.time.file.names[[time]]
    # cat("fcs.file.names are", basename(fcs.file.names), "\n")
    if (length(subsamp.orig) == 1 & subsamp.orig != FALSE) {
      cat("Subsampling all files to:", subsamp.orig, "\n")
      subsample.new <- rep(subsamp.orig, times = length(fcs.file.names))
      subsamples <- subsample.new
    } else {
      subsamples <- subsamp.orig[[time]]
    }
    list.of.time.clean.FCS.files[[time]] <- LoadCleanFCS(fcs.file.names,
                                                         channel.remove,
                                                         channel.annotate,
                                                         subsamples,
                                                         subsample.rand,
                                                         transform, scale)
    f.names <- c() 
    for (i in 1:length(list.of.time.clean.FCS.files[[time]])) {
      Time <- rep(as.numeric(time), times = dim(list.of.time.clean.FCS.files[[time]][[i]])[1])
      list.of.time.clean.FCS.files[[time]][[i]] <- cbind.data.frame(list.of.time.clean.FCS.files[[time]][[i]],
                                                                    Time, stringsAsFactors = FALSE)
      this.name <- basename(fcs.file.names[i])
      this.name <- gsub(".fcs", "", this.name)
      if (grepl("-", this.name)) {
        this.name <- unlist(strsplit(this.name, split = "-"))[1]
      }
      if (grepl("\\.", this.name)) {
        this.name <- unlist(strsplit(this.name, split = "\\."))[1]
      }
      Condition <- rep(this.name, times = dim(list.of.time.clean.FCS.files[[time]][[i]])[1])
      list.of.time.clean.FCS.files[[time]][[i]] <- cbind.data.frame(list.of.time.clean.FCS.files[[time]][[i]],
                                                                    Condition, stringsAsFactors = FALSE)
      f.names <- c(f.names, this.name)
      rm(Time, Condition)
    }
    names(list.of.time.clean.FCS.files[[time]]) <- f.names
  }
  return(list.of.time.clean.FCS.files)
}

ConvertVariables <- function(clustering.var, var.annotate) {
  ind <- match(clustering.var, var.annotate)
  fixed.clustering.var <- names(var.annotate)[ind]
  return(fixed.clustering.var)
}

# SPADE
# density-dependent downsampling

DownsampleFCS <- function(fcs.files, clustering.var, distance.metric,
                          exclude.pctile = 0.01, target.pctile = 0.99,
                          target.number = NULL, target.percent = 0.1,
                          ...) {
  
  # Eli: In this study, the SPADE package was used
  # to perform density-dependent downsampling and hierarchical
  # clustering as previously described (Linderman et al., 2012;
  # Qiu et al., 2011). Parameters used for clustering were
  # Downsampling Exclude Percentile = 0.01, Downsampling Target
  # Percentile = 0.99, and Target Clusters = 200. 
  
  # optional variables
  # transforms
  
  downsample.files <- c()
  for (file.name in fcs.files) {
    base.name <- unlist(strsplit(basename(file.name), "\\."))[1]
    infilename <- paste(base.name, "density.fcs", sep = "_")
    spade::SPADE.addDensityToFCS(file.name, infilename,
                                 cols = ConvertVariables(clustering.var, var.annotate), comp = TRUE,
                                 transforms = flowCore::arcsinhTransform(a = 0, b = 0.2))
    outfilename <- paste(base.name, "downsample.fcs", sep = "_")
    spade::SPADE.downsampleFCS(infilename = infilename, outfilename,
                               exclude_pctile = exclude.pctile,
                               target_pctile = target.pctile,
                               target_number = target.number,
                               target_percent = target.percent)
    downsample.files <- c(downsample.files, outfilename)
  }
  return(downsample.files)
}

ConvertNumericLabel <- function(list.of.clean.FCS.files.with.labels) {
  label.key <- list()
  fixed.list.FCS.files <- list()
  for (i in 1:length(list.of.clean.FCS.files.with.labels)) {
    label.key[[i]] <- names(list.of.clean.FCS.files.with.labels[[i]])
    tmp.label.fcs <- list.of.clean.FCS.files.with.labels[[i]]
    for (n in 1:length(tmp.label.fcs)) {
      inds <- which(tmp.label.fcs[[n]] == label.key[[i]], arr.ind = TRUE)
      col.ind <- unname(inds[1, 2])
      tmp.label.fcs[[n]][, col.ind] <- as.numeric(n)
    }
    fixed.list.FCS.files[[i]] <- tmp.label.fcs
  }
  names(fixed.list.FCS.files) <- names(list.of.clean.FCS.files.with.labels)
  return(list(fixed.list.FCS.files = fixed.list.FCS.files,
              label.key = label.key))
}

ConvertCharacterLabel <- function(data.frame.with.numeric.labels, label.key) {
  data.frame.with.character.labels <- data.frame.with.numeric.labels
  print("label.key")
  print(label.key)
  times <- unique(data.frame.with.numeric.labels[, "Time"])
  for (t in times) {
    this.label <- label.key[[t]]
    this.ind <- which(data.frame.with.numeric.labels[, "Time"] == t)
    for (i in 1:length(this.label)) {
      fix.ind <- which(data.frame.with.numeric.labels[, "Condition"] == i)
      use.ind <- intersect(fix.ind, this.ind)
      data.frame.with.character.labels[use.ind, "Condition"] <- this.label[i]
    }
  }
  return(data.frame.with.character.labels)
}
