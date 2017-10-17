
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

RemoveExistingTimeVar <- function(fcs.file) {
  fcs.file <- as.data.frame(fcs.file)
  ind <- grep(colnames(fcs.file), pattern = "Time")
  fcs.file[, ind] <- NULL
  return(fcs.file)
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
                         subsamples = 1000, transform = TRUE) {
  clean.fcs.files <- list()
  if (length(subsamples) == 1 & subsamples != FALSE) {
    cat("Subsampling all files to:", subsamples, "\n")
    subsample.new <- rep(subsamples, times = length(fcs.file.names))
    subsamples <- subsample.new
  }
  for (i in 1:length(fcs.file.names)) {
    current.file <- tail(strsplit(fcs.file.names[i], "/")[[1]], n = 1)
    cat("Reading FCS file data from:", current.file, "\n")
    # store currently read FCS file
    if (subsamples == FALSE) {
      fcs.file <- read.FCS(fcs.file.names[i])  
      fcs.file <- as.data.frame(exprs(fcs.file))
    } else { 
      cat("Subsampling", current.file, "to", subsamples[i], "cells\n")
      fcs.file <- read.FCS(fcs.file.names[i], which.lines = subsamples[i])
      fcs.file <- as.data.frame(exprs(fcs.file))
    }
    # rename variables with protein marker measured instead of metal channel
    cat("Fixing channel names from:", current.file, "\n")
    for (x in 1:length(colnames(fcs.file))) {
      if (exists(colnames(fcs.file)[x], where = channel.annotate)) {
        colnames(fcs.file)[x] <- channel.annotate[[colnames(fcs.file)[x]]]
      }
    }
    # remove unneeded variables
    cat("Removing unnecessary channel names from:", current.file, "\n")
    fcs.file <- subset(fcs.file, select = colnames(fcs.file)[!colnames(fcs.file) %in% channel.remove])
    if (transform) {
      cat("Transforming data from:", current.file, "\n")
      fcs.file <- apply(fcs.file, 2, Asinh) 
    }
    fcs.file <- as.data.frame(fcs.file)
    fcs.file <- RemoveExistingTimeVar(fcs.file) 
    clean.fcs.files[[i]] <- fcs.file
    rm(fcs.file)
  }
  return(clean.fcs.files)
}

LoadMultiCleanFCS <- function(list.of.file.names, channel.remove, channel.annotate,
                              subsamples = 1000, transform = TRUE) {
  list.of.FCS.files <- list()
  subsamp.orig <- subsamples
  print("subsamp.orig")
  print(subsamp.orig)
  for (t in 1:length(list.of.file.names)) {
    fcs.file.names <- list.of.file.names[[t]]
    if (length(subsamp.orig) == 1 & subsamp.orig != FALSE) {
      cat("Subsampling all files to:", subsamp.orig, "\n")
      subsample.new <- rep(subsamp.orig, times = length(fcs.file.names))
      subsamples <- subsample.new
    } else if (length(subsamp.orig) > 1 & subsamp.orig != FALSE) {
      subsamples <- subsamp.orig[[t]]
    }
    list.of.FCS.files[[t]] <- LoadCleanFCS(fcs.file.names,
                                           channel.remove,
                                           channel.annotate,
                                           subsamples,
                                           transform)
    f.names <- c() 
    for (i in 1:length(list.of.FCS.files[[t]])) {
      Time <- rep(as.numeric(t), times = dim(list.of.FCS.files[[t]][[i]])[1])
      list.of.FCS.files[[t]][[i]] <- cbind.data.frame(list.of.FCS.files[[t]][[i]],
                                                      Time, stringsAsFactors = FALSE)
      this.name <- basename(fcs.file.names[i])
      this.name <- gsub(".fcs", "", this.name)
      if (grepl("-", this.name)) {
        this.name <- unlist(strsplit(this.name, split = "-"))[1]
      }
      if (grepl("\\.", this.name)) {
        this.name <- unlist(strsplit(this.name, split = "\\."))[1]
      }
      Condition <- rep(this.name, times = dim(list.of.FCS.files[[t]][[i]])[1])
      list.of.FCS.files[[t]][[i]] <- cbind.data.frame(list.of.FCS.files[[t]][[i]],
                                                      Condition, stringsAsFactors = FALSE)
      f.names <- c(f.names, this.name)
      rm(Time, Condition)
    }
    names(list.of.FCS.files[[t]]) <- f.names
  }
  return(list.of.FCS.files)
}

ConvertVariables <- function(clustering.var, var.annotate) {
  ind <- match(clustering.var, var.annotate)
  fixed.clustering.var <- names(var.annotate)[ind]
  return(fixed.clustering.var)
}

#' SPADE: density-dependent downsampling
#' 
#' This code was adapted from the spade R package,
#' available here: https://github.com/nolanlab/spade.
#' 
#' Specifically, the code in this file is adapted from downsample.R.
#' 
#' This code, authored by M. Linderman, P. Qiu, E. Simonds, Z. Bjornson,
#' and maintained by Michael Linderman <michael.d.linderman@gmail.com>,
#' was used under the GNU General Public License v. 2.0
#' (https://github.com/nolanlab/spade/blob/master/LICENSE) available
#' here: https://opensource.org/licenses/GPL-2.0. In accordance with these
#' license rules, our code is available under GPL-3.0.

DownsampleFCS <- function(fcs.file.names, clustering.var, channel.annotate,
                          channel.remove, exclude.pctile = 0.01, target.pctile = 0.99,
                          target.number = NULL, target.percent = 0.1,
                          transform = TRUE) {
  cat("exclude.pctile is", exclude.pctile, "\n")
  cat("target.pctile is", target.pctile, "\n")
  cat("target.number is", target.number, "\n")
  cat("target.percent is", target.percent, "\n")
  downsample.data <- list()
  for (file.name in fcs.file.names) {
    transforms <- flowCore::arcsinhTransform(a = 0, b = 0.2)
    SPADE.removeExistingDensityAndClusterColumns(file.name)
    current.file <- tail(strsplit(file.name, "/")[[1]], n = 1)
    cat("Reading FCS file data from:", current.file, "\n")
    fcs.file <- read.FCS(file.name)  
    fcs.file <- as.data.frame(exprs(fcs.file))
    cat("Fixing channel names from:", current.file, "\n")
    for (x in 1:length(colnames(fcs.file))) {
      if (exists(colnames(fcs.file)[x], where = channel.annotate)) {
        colnames(fcs.file)[x] <- channel.annotate[[colnames(fcs.file)[x]]]
      }
    }
    cat("Removing unnecessary channel names from:", current.file, "\n")
    fcs.file <- subset(fcs.file, select = colnames(fcs.file)[!colnames(fcs.file) %in% channel.remove])
    if (transform) {
      cat("Transforming data from:", current.file, "\n")
      fcs.file <- apply(fcs.file, 2, Asinh) 
    }
    fcs.file <- RemoveExistingTimeVar(fcs.file) 
    print(dim(fcs.file))
    cat("Calculating density for:", current.file, "\n")
    density <- SPADE.density(fcs.file[, clustering.var], kernel_mult = 5.0, apprx_mult = 1.5, med_samples = 2000)
    if (max(density) == 0.0) {
      warning(paste(current.file, "has degenerate densities, possibly due to many identical observations", sep = " "))
    }
    fcs.file <- cbind(fcs.file, "density" = density)
    boundary <- quantile(fcs.file[, "density"], c(exclude.pctile, target.pctile), names = FALSE)
    print(dim(fcs.file))
    cat("Removing outliers for:", current.file, "\n")
    fcs.file2 <- subset(fcs.file, fcs.file[, "density"] > boundary[1])    
    print(dim(fcs.file2))
    if (!is.null(target.percent)) {
      target.number = round(target.percent * nrow(fcs.file2))
      cat("Targeting", target.number, "events for", current.file, "\n")
    }
    density <- fcs.file2[, "density"]
    if (is.null(target.number)) {
      cat("Downsampling for:", current.file, "\n")
      boundary <- boundary[2]
      fcs.file3 <- subset(fcs.file2, boundary/density > runif(nrow(fcs.file2)))
      print(dim(fcs.file3))
    } else if (target.number < nrow(fcs.file2)) {
      sorted.density <- sort(density)
      cdf <- rev(cumsum(1.0 / rev(sorted.density)))
      boundary <- target.number/cdf[1] 
      if (boundary > sorted.density[1]) {  # Boundary actually falls amongst densities present
        targets <- (target.number - 1:length(sorted.density)) / cdf 
        boundary <- targets[which.min(targets - sorted.density > 0)]
      }
      cat("Downsampling for:", current.file, "\n")
      fcs.file3 <- subset(fcs.file2, boundary/density > runif(length(density)))
      print(dim(fcs.file3))
    } else if (target.number > nrow(fcs.file2)) {
      stop("More events requested than present in file")
    }
    downsample.data[[current.file]] <- fcs.file3
    rm(fcs.file, fcs.file2, fcs.file3)
  }
  return(downsample.data)
}


MultiDownsampleFCS <- function(fcs.file.names, clustering.var, channel.annotate,
                               channel.remove, exclude.pctile = 0.01, target.pctile = 0.99,
                               target.number = NULL, target.percent = 0.1,
                               transform = TRUE) {
  list.downsample.data <- list()
  for (t in 1:length(list.of.file.names)) {
    fcs.file.names.t <- fcs.file.names[[t]]
    downsample.data <- DownsampleFCS(fcs.file.names.t, clustering.var, channel.annotate,
                                     channel.remove, exclude.pctile, target.pctile,
                                     target.number, target.percent, transform)
    list.downsample.data[[t]] <- downsample.data
  }
  return(list.downsample.data)
}

#' SPADE: density-dependent downsampling
#' 
#' This code was adapted from the spade R package,
#' available here: https://github.com/nolanlab/spade.
#' 
#' Specifically, the code in this file is adapted from downsample.R.
#' 
#' This code, authored by M. Linderman, P. Qiu, E. Simonds, Z. Bjornson,
#' and maintained by Michael Linderman <michael.d.linderman@gmail.com>,
#' was used under the GNU General Public License v. 2.0
#' (https://github.com/nolanlab/spade/blob/master/LICENSE) available
#' here: https://opensource.org/licenses/GPL-2.0. In accordance with these
#' license rules, our code is available under GPL-3.0.

DownsampleFCS_OLD2 <- function(fcs.file.names, clustering.var, var.annotate,
                               distance.metric, exclude.pctile = 0.01,
                               target.pctile = 0.99,
                               target.number = NULL,
                               target.percent = 0.1) {
  cat("exclude_pctile is", exclude.pctile, "\n")
  cat("target_pctile is", target.pctile, "\n")
  cat("target_number is", target.number, "\n")
  cat("target_percent is", target.percent, "\n")
  downsample.file.names <- c()
  for (file.name in fcs.file.names) {
    base.name <- unlist(strsplit(basename(file.name), "\\."))[1]
    infilename <- paste(base.name, "density.fcs", sep = "_")
    transforms <- flowCore::arcsinhTransform(a = 0, b = 0.2)
    new.cols <- c()
    for (i in clustering.var) {
      new.name <- names(which(var.annotate == i))
      new.cols <- c(new.cols, new.name)
    }
    SPADE.addDensityToFCS(file.name, infilename, cols = new.cols, comp = TRUE, transforms = transforms)
    print("infilename")
    print(infilename)
    outfilename <- paste(base.name, "downsample.fcs", sep = "_")
    SPADE.downsampleFCS(infilename = infilename, outfilename,
                        exclude_pctile = exclude.pctile,
                        target_pctile = target.pctile,
                        target_number = target.number,
                        target_percent = target.percent)
    downsample.file.names <- c(downsample.file.names, outfilename)
  }
  return(downsample.file.names)
}

MultiDownsampleFCS_OLD2 <- function(list.of.file.names, clustering.var,
                                    distance.metric, exclude.pctile = exclude.pctile,
                                    target.pctile = target.pctile,
                                    target.number = target.number,
                                    target.percent = target.percent) {
  list.downsample.file.names <- list()
  for (t in 1:length(list.of.file.names)) {
    fcs.file.names <- list.of.file.names[[t]]
    downsample.file.names <- DownsampleFCS(fcs.file.names, clustering.var,
                                           distance.metric, exclude.pctile = exclude.pctile,
                                           target.pctile = target.pctile,
                                           target.number = target.number,
                                           target.percent = target.percent)
    list.downsample.file.names[[t]] <- downsample.file.names
  }
  return(list.downsample.file.names)
}

#' SPADE: density-dependent downsampling
#' 
#' This code was adapted from the spade R package,
#' available here: https://github.com/nolanlab/spade.
#' 
#' Specifically, the code in this file is adapted from downsample.R, but
#' has been modified to not produce FCS files as a result.
#' 
#' This code, authored by M. Linderman, P. Qiu, E. Simonds, Z. Bjornson,
#' and maintained by Michael Linderman <michael.d.linderman@gmail.com>,
#' was used under the GNU General Public License v. 2.0
#' (https://github.com/nolanlab/spade/blob/master/LICENSE) available
#' here: https://opensource.org/licenses/GPL-2.0. In accordance with these
#' license rules, our code is available under GPL-3.0.

DownsampleFCS_OLD <- function(fcs.file.names, clustering.var, channel.annotate,
                              channel.remove, exclude.pctile = 0.01, target.pctile = 0.99,
                              target.number = NULL, target.percent = 0.1,
                              transform = TRUE) {
  downsample.files <- list()
  cat("Downsampling FCS files.", "\n")
  for (i in 1:length(fcs.file.names)) {
    out.file <- read.FCS(fcs.file.names[i], transformation = "linearize", which.lines = NULL,
                         alter.names = FALSE, column.pattern = NULL,
                         invert.pattern = FALSE, decades = 0, ncdf = FALSE,
                         min.limit = NULL, dataset = NULL, emptyValue = TRUE)
    out.file <- head(out.file, n = dim(out.file)[1])
    # rename variables with protein marker measured instead of metal channel
    for (x in 1:length(colnames(out.file))) {
      if (exists(colnames(out.file)[x], where = channel.annotate)) {
        colnames(out.file)[x] <- channel.annotate[[colnames(out.file)[x]]]
      }
    }
    # remove unneeded variables
    out.file <- subset(out.file, select = colnames(out.file)[!colnames(out.file) %in% channel.remove])
    # Compute the density
    idxs <- match(clustering.var, colnames(out.file))
    # density <- spade::SPADE.density(out.file[, idxs], ...)
    density <- spade::SPADE.density(out.file[, idxs])
    if (max(density) == 0.0) { stop("File has degenerate densities, possibly due to many identical observations!") }
    # Add column named "density" to the FCS file
    if (transform) {
      out.file <- apply(out.file, 2, Asinh) 
    }
    out.file <- cbind(out.file, density = density)
    boundary <- quantile(density, c(exclude.pctile, target.pctile), names = FALSE)
    out.file <- subset(out.file, out.file[, "density"] > boundary[1]) # Exclusion
    if (!is.null(target.percent)) {
      target.number = round(target.percent * nrow(out.file))
      cat("Targeting", target.number, "events for file", i, "\n")
    }
    if (is.null(target.number)) {
      boundary <- boundary[2]
      out.file <- subset(out.file, ((boundary / out.file[, "density"]) > runif(nrow(out.file))))
    } else if (target.number < nrow(out.file)) {
      density.sort <- sort(out.file[, "density"])
      cdf <- rev(cumsum(1.0 / rev(density.sort)))
      # Default solution if target density smaller than any present
      boundary <- target.number/cdf[1]
      if (boundary > density.sort[1]) {  # Boundary actually falls amongst densities present
        targets <- (target.number-1:length(density.sort)) / cdf
        boundary <- targets[which.min(targets - density.sort > 0)]
      }
      out.file  <- subset(out.file, ((boundary / out.file[, "density"]) > runif(length(out.file[, "density"]))))
    } else if (target.number > nrow(out.file)) {
      stop("More events requested than present in file!")
    }
    downsample.files[[i]] <- out.file
  }
  return(downsample.files)
}

MultiDownsampleFCS_OLD <- function(list.of.file.names, clustering.var, channel.annotate,
                                   channel.remove, exclude.pctile = 0.01, target.pctile = 0.99,
                                   target.number = NULL, target.percent = 0.1,
                                   transform = TRUE) {
  list.of.downsample.files <- list()
  for (t in 1:length(list.of.file.names)) {
    fcs.file.names <- list.of.file.names[[t]]
    fcs.files <- DownsampleFCS(fcs.file.names, clustering.var, channel.annotate,
                               channel.remove, exclude.pctile, target.pctile,
                               target.number, target.percent, transform) 
    list.of.downsample.files[[t]] <- fcs.files
    f.names <- c() 
    for (i in 1:length(list.of.downsample.files[[t]])) {
      Time <- rep(as.numeric(t), times = dim(list.of.downsample.files[[t]][[i]])[1])
      list.of.downsample.files[[t]][[i]] <- cbind.data.frame(list.of.downsample.files[[t]][[i]],
                                                             Time, stringsAsFactors = FALSE)
      this.name <- basename(fcs.file.names[i])
      this.name <- gsub(".fcs", "", this.name)
      if (grepl("-", this.name)) {
        this.name <- unlist(strsplit(this.name, split = "-"))[1]
      }
      if (grepl("\\.", this.name)) {
        this.name <- unlist(strsplit(this.name, split = "\\."))[1]
      }
      Condition <- rep(this.name, times = dim(list.of.downsample.files[[t]][[i]])[1])
      list.of.downsample.files[[t]][[i]] <- cbind.data.frame(list.of.downsample.files[[t]][[i]],
                                                             Condition, stringsAsFactors = FALSE)
      f.names <- c(f.names, this.name)
      rm(Time, Condition)
    }
    names(list.of.downsample.files[[t]]) <- f.names
  }
  return(list.of.downsample.files)
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
