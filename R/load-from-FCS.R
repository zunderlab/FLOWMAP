
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

#' GetFCSNames
#'
#' \code{GetFCSNames} returns the FCS file names from a given folder.
#'
#' @param folder the full path for a folder that contains subfolders, each
#' containing the FCS files from a specific timepoint to be used
#' @param sort Logical specifying whether FCS files should be sorted in alphanumeric order
#' @return Vector with the file paths of FCS files found in the given folder
#' @examples
#' folder <- "Desktop/FCS_Files"
#' GetFCSNames(folder, sort = TRUE)
#' @export
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

#' GetMultiFCSNames
#'
#' \code{GetFCSNames} returns the FCS file names as a list from a given folder,
#' that contains subfolders corresponding to each timepoint.
#'
#' @param folder the full path for a folder that contains the FCS files to be used
#' @param sort Logical specifying whether FCS files should be sorted in alphanumeric order
#' @return List with each member containing a vector with the file paths
#' of FCS files found in the given subfolder
#' @examples
#' folder <- "Desktop/FCS_Files"
#' GetMultiFCSNames(folder, sort = TRUE)
#' @export
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

#' LoadCleanFCS
#'
#' \code{LoadCleanFCS} returns a list of FCS file data, where each member is
#' a dataframe from a different timepoint. In the mode \code{"single"}, the list has
#' only a single element, a dataframe from the one FCS file path provided.
#'
#' @param fcs.file.names A vector of full file paths to the FCS files to be used
#' @param channel.remove Vector naming channels to be removed from all loaded FCS data
#' @param channel.annotate List mapping channel names to user-specified names to properly
#' annotate all FCS file data
#' @param subsamples Numeric or vector of numerics specifying how many cells to sample
#' from each FCS file, default is set to \code{1000}
#' @param transform Logical specifying whether to transform the data using an Asinh
#' transform typical of CyTOF/mass cytometry datasets, default is set to \code{TRUE}
#' @return a list where each member is a dataframe containing the single-cell data
#' @examples
#' fcs.file.names <- c("Desktop/A.fcs", "Desktop/B.fcs")
#' var.remove <- c("Channel3", "Channel4")
#' var.annotate <- list("c1" = "Channel1", "c2" = "Channel2",
#' "c3" = "Channel3", "c4" = "Channel4")
#'
#' LoadCleanFCS(fcs.file.names, var.remove, var.annotate, subsamples = 100, transform = TRUE)
#' @export
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

#' LoadMultiCleanFCS
#'
#' \code{LoadMultiCleanFCS} returns a list of FCS file data, where each member is
#' a list that contains dataframe from different conditions from different timepoints.
#'
#' @param list.of.file.names A list where each member contains a vector of full file paths
#' to the FCS files to be used
#' @param channel.remove Vector naming channels to be removed from all loaded FCS data
#' @param channel.annotate List mapping channel names to user-specified names to properly
#' annotate all FCS file data
#' @param subsamples Numeric or vector of numerics specifying how many cells to sample
#' from each FCS file, default is set to \code{1000}
#' @param transform Logical specifying whether to transform the data using an Asinh
#' transform typical of CyTOF/mass cytometry datasets, default is set to \code{TRUE}
#' @return a list where each member is a list containing multiple dataframes containing
#' single-cell data from the same timepoint, but multiple conditions
#' @examples
#' list.of.file.names <- list("1" = c("Desktop/A-1.fcs", "Desktop/B-1.fcs"),
#' "2" = c("Desktop/A-2.fcs", "Desktop/B-2.fcs"))
#' var.remove <- c("Channel3", "Channel4")
#' var.annotate <- list("c1" = "Channel1", "c2" = "Channel2",
#' "c3" = "Channel3", "c4" = "Channel4")
#'
#' LoadMultiCleanFCS(list.of.file.names, var.remove, var.annotate, subsamples = 100, transform = TRUE)
#' @export
LoadMultiCleanFCS <- function(list.of.file.names, channel.remove, channel.annotate,
                              subsamples = 1000, transform = TRUE) {
  list.of.FCS.files <- list()
  subsamp.orig <- subsamples
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
#' available here: \url{https://github.com/nolanlab/spade}.
#'
#' Specifically, the code in this file is adapted from downsample.R.
#'
#' This code, authored by M. Linderman, P. Qiu, E. Simonds, Z. Bjornson,
#' and maintained by Michael Linderman <michael.d.linderman@gmail.com>,
#' was used under the GNU General Public License v. 2.0
#' (\url{https://github.com/nolanlab/spade/blob/master/LICENSE}) available
#' here: \url{https://opensource.org/licenses/GPL-2.0}. In accordance with these
#' license rules, our code is available under GPL-3.0.
#'
#' UPDATE 12/17/19: Spade is no longer maintained
#' We have included relevant function here to maintain functionality of FLOWMAP
#'
#' @export

DownsampleFCS <- function(fcs.file.names, clustering.var, channel.annotate,
                          channel.remove, exclude.pctile = 0.01, target.pctile = 0.99,
                          target.number = NULL, target.percent = 0.1,
                          transform = TRUE, k=15) {
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
    cat("Calculating density for:", current.file, "\n")
    #############density <- SPADE.density(fcs.file[, clustering.var], kernel_mult = 5.0, apprx_mult = 1.5, med_samples = 2000)
    nns <- RANN::nn2(data=fcs.file[, clustering.var], k=k+1, searchtype="priority", eps=0.1)
    temp_nnids.df <- as.data.frame(nns$nn.idx)
    temp_nndists.df <- as.data.frame(nns$nn.dists)
    nn.ids.df <- temp_nnids.df[,2:length(temp_nnids.df)]
    nn.dists.df <- temp_nndists.df[,2:length(temp_nndists.df)]
    numcluster <- nrow(clusters)
    #TODO What to set for k?
    density <- KnnDensity(k=k, min, max, n=0,nn.ids.df = nn.ids.df,
                          nn.dists.df = nn.dists.df,numcluster = numcluster,
                          table.breaks = NULL,offset = 0)
    #############density <- SPADE.density(fcs.file[, clustering.var], kernel_mult = 5.0, apprx_mult = 1.5, med_samples = 2000)
    if (max(density) == 0.0) {
      warning(paste(current.file, "has degenerate densities, possibly due to many identical observations", sep = " "))
    }
    fcs.file <- cbind(fcs.file, "density" = density)
    boundary <- quantile(fcs.file[, "density"], c(exclude.pctile, target.pctile), names = FALSE)
    cat("Removing outliers for:", current.file, "\n")
    fcs.file2 <- subset(fcs.file, fcs.file[, "density"] > boundary[1])
    if (!is.null(target.percent)) {
      target.number = round(target.percent * nrow(fcs.file2))
      cat("Targeting", target.number, "events for", current.file, "\n")
    }
    density <- fcs.file2[, "density"]
    if (is.null(target.number)) {
      cat("Downsampling for:", current.file, "\n")
      boundary <- boundary[2]
      fcs.file3 <- subset(fcs.file2, boundary/density > runif(nrow(fcs.file2)))
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
  for (t in 1:length(fcs.file.names)) {
    fcs.file.names.t <- fcs.file.names[[t]]
    downsample.data <- DownsampleFCS(fcs.file.names.t, clustering.var, channel.annotate,
                                     channel.remove, exclude.pctile, target.pctile,
                                     target.number, target.percent, transform)
    list.downsample.data[[t]] <- downsample.data

    f.names <- c()
    for (i in 1:length(list.downsample.data[[t]])) {
      Time <- rep(as.numeric(t), times = dim(list.downsample.data[[t]][[i]])[1])
      list.downsample.data[[t]][[i]] <- cbind.data.frame(list.downsample.data[[t]][[i]],
                                                         Time, stringsAsFactors = FALSE)
      this.name <- basename(fcs.file.names.t[i])
      this.name <- gsub(".fcs", "", this.name)
      if (grepl("-", this.name)) {
        this.name <- unlist(strsplit(this.name, split = "-"))[1]
      }
      if (grepl("\\.", this.name)) {
        this.name <- unlist(strsplit(this.name, split = "\\."))[1]
      }
      Condition <- rep(this.name, times = dim(list.downsample.data[[t]][[i]])[1])
      list.downsample.data[[t]][[i]] <- cbind.data.frame(list.downsample.data[[t]][[i]],
                                                         Condition, stringsAsFactors = FALSE)
      f.names <- c(f.names, this.name)
      rm(Time, Condition)
    }
    names(list.downsample.data[[t]]) <- f.names
  }
  return(list.downsample.data)
}

SPADE.removeExistingDensityAndClusterColumns <- function(file) {
  # Do not comp or transform ... make this step invisible.
  input_file <- suppressWarnings(read.FCS(file))

  input_file_names <- names(input_file)

  if ("<cluster> cluster" %in% input_file_names ||
      "<density> density" %in% input_file_names) {
    # Drop those columns
    cleaned <- input_file[,!(input_file_names %in% c("<cluster> cluster", "<density> density"))]

    # Rename the original file. Increment the suffix if it already exists.
    suffix <- as.integer(regmatches(file, gregexpr("(\\d+)$", file, perl=TRUE)))

    if (is.na(suffix)) {
      suffix <- ".orig1"
      new_file_name <- paste(file, suffix, sep = "")
      # NB: file.rename is a namespaced function, not a method of the argument.
      file.rename(file, new_file_name)
    } else {
      suffix <- paste(".orig", suffix + 1, sep = "")
      new_file_name <- sub("(.orig\\d+)$", suffix, file);
      file.rename(file, new_file_name)
    }

    # Save with the original file name
    write.FCS(cleaned, file)
  }
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
  for (t in 1:length(times)) {
    # for (t in times) {
    # this.label <- label.key[[t]]
    this.label <- label.key[[t]]
    this.ind <- which(data.frame.with.numeric.labels[, "Time"] == times[t])
    # this.ind <- which(data.frame.with.numeric.labels[, "Time"] == t)
    for (i in 1:length(this.label)) {
      fix.ind <- which(data.frame.with.numeric.labels[, "Condition"] == i)
      use.ind <- intersect(fix.ind, this.ind)
      data.frame.with.character.labels[use.ind, "Condition"] <- this.label[i]
    }
  }
  return(data.frame.with.character.labels)
}

ConvertCharacterLabelSpecial <- function(data.frame.with.numeric.labels, label.key.special) {
  data.frame.with.character.labels <- data.frame.with.numeric.labels
  for (i in 1:length(label.key.special)) {
    fix.ind <- which(data.frame.with.numeric.labels[, "Condition"] == i)
    data.frame.with.character.labels[fix.ind, "Condition"] <- label.key.special[i]
  }
  return(data.frame.with.character.labels)
}
