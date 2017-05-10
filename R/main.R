
ShuffleCells <- function(fcs.files, subsamples) {
  x <- c()
  for (i in 1:length(fcs.files)) {
    if (length(subsamples) > 1) {
      subsamp <- subsamples[i]
    } else if (subsamples == FALSE) {
      subsamp <- nrow(fcs.files[[i]])
    } else {
      subsamp <- subsamples
    }
    df1 <- fcs.files[[i]]
    x <- c(x, nrow(df1))
    df2 <- df1[sample(nrow(df1)), ]
    fcs.files[[i]] <- df2
    rownames(fcs.files[[i]]) <- seq(1:subsamp)
  }
  return(fcs.files)
}

MultiShuffleCells <- function(fcs.files, subsamples) {
  x <- c()
  for (n in 1:length(fcs.files)) {
    for (i in 1:length(fcs.files[[n]])) {
      if (length(subsamples) > 1) {
        subsamp <- subsamples[[n]][[i]]
      } else if (subsamples == FALSE) {
        subsamp <- nrow(fcs.files[[n]][[i]])
      } else {
        subsamp <- subsamples
      }
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
      if (all(grepl(pattern = "\\.fcs", list.files(files)))) {
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
      subfolders <- list.files(files)
      temp.fail.flag <- FALSE
      for (n in subfolders) {
        if (!all(grepl(pattern = "\\.fcs", list.files(n)))) {
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
        if (!all(grepl(pattern = "\\.fcs", files[[n]]))) {
          temp.fail.flag <- TRUE
        }
        if (length(files[[n]]) > 1) {
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


OutputRFile <- function() {
  
}


#' @export
FLOWMAP <- function(seed.X, files, var.remove, var.annotate, clustering.var,
                    cluster.numbers, subsamples, distance.metric,
                    minimum, maximum, per, save.folder, mode = c("single", "multi", "one"),
                    name.sort = TRUE, downsample = TRUE,
                    savePDFs = TRUE, which.palette = "bluered", ...) {
  # optional variables
  # transform 
  # scale
  # starting.files = c("FCS", "cluster_matrix")
  set.seed(seed.X)
  setwd(save.folder)
  if (mode == "single") {
    check <- CheckModeSingle(files)
    cat("check", check, "\n")
    if (check[1]) {
      stop("Unknown 'files' format provided for specified mode!")
    }
    runtype <- "SingleFLOWMAP"
    output.folder <- MakeOutFolder(runtype = runtype)
    setwd(output.folder)
    if (check[2] == "FCS") {
      fcs.file.names <- files
    }
    if (check[2] == "folder") {
      fcs.file.names <- GetFCSNames(folder = files, sort = name.sort)
    }
    orig.times <- ParseTimes(fcs.file.names, name.sort = name.sort)
    file.name <- fcs.file.names[1]
    if (downsample) {
      cat("Downsampling all files using SPADE downsampling", "\n")
      fcs.files <- DownsampleFCS(fcs.file.names, clustering.var, channel.annotate = var.annotate,
                                 channel.remove = var.remove, exclude.pctile = exclude.pctile,
                                 target.pctile = target.pctile, target.number = target.number,
                                 target.percent = target.percent)
    } else {
      fcs.files <- LoadCleanFCS(fcs.file.names = fcs.file.names, channel.remove = var.remove,
                                channel.annotate = var.annotate, subsamples = subsamples)
    }
    fcs.files <- ShuffleCells(fcs.files, subsamples)
    file.clusters <- ClusterFCS(fcs.files = fcs.files, clustering.var = clustering.var,
                                numcluster = cluster.numbers, distance.metric = distance.metric)
    if (downsample) {
      cat("Upsampling all clusters to reflect Counts prior to SPADE downsampling", "\n")
      file.clusters <- Upsample(file.clusters)
    }
    results <- BuildFLOWMAP(FLOWMAP.clusters = file.clusters, per = per, min = minimum,
                            max = maximum, distance.metric = distance.metric, cellnum = subsamples,
                            clustering.var = clustering.var)
    graph <- results$output.graph
  } else if (mode == "multi") {
    check <- CheckModeMulti(files)
    cat("check", check, "\n")
    if (check[1]) {
      stop("Unknown 'files' format provided for specified mode!")
    }
    runtype <- "MultiFLOWMAP"
    output.folder <- MakeOutFolder(runtype = runtype)
    setwd(output.folder)
    if (check[2] == "list") {
      fcs.file.names <- files
      orig.times <- MultiListParseTimes(fcs.file.names, name.sort)
    }
    if (check[2] == "subfolder") {
      fcs.file.names <- GetMultiFCSNames(folder = files, sort = name.sort)
      orig.times <- MultiFolderParseTimes(files, fcs.file.names, name.sort = name.sort)
    }
    file.name <- fcs.file.names[[1]][1]
    if (downsample) {
      cat("Downsampling all files using SPADE downsampling", "\n")
      fcs.files <- MultiDownsampleFCS(fcs.file.names, clustering.var, channel.annotate = var.annotate,
                                      channel.remove = var.remove, exclude.pctile = exclude.pctile,
                                      target.pctile = target.pctile, target.number = target.number,
                                      target.percent = target.percent)
    } else {
      fcs.files <- LoadMultiCleanFCS(fcs.file.names, var.remove, var.annotate,
                                     subsamples = subsamples)
    }
    fcs.files <- MultiShuffleCells(fcs.files, subsamples)
    fcs.files.conversion <- ConvertNumericLabel(fcs.files)
    fixed.fcs.files <- fcs.files.conversion$fixed.list.FCS.files
    label.key <- fcs.files.conversion$label.key
    file.clusters <- MultiClusterFCS(fixed.fcs.files, clustering.var = clustering.var, numcluster = cluster.numbers,
                                     distance.metric = distance.metric)
    if (downsample) {
      cat("Upsampling all clusters to reflect Counts prior to SPADE downsampling", "\n")
      file.clusters <- MultiUpsample(file.clusters)
    }
    graph <- BuildMultiFLOWMAP(file.clusters, per = per, min = minimum,
                               max = maximum, distance.metric = distance.metric, cellnum = subsamples,
                               label.key = label.key)
  } else if (mode == "one") {
    check <- CheckModeOne(files)
    cat("check", check, "\n")
    if (check[1]) {
      stop("Unknown 'files' format provided for specified mode!")
    }
    if (check[2] == "FCS") {
      fcs.file.names <- files
    }
    if (check[2] == "folder") {
      fcs.file.names <- GetFCSNames(folder = files, sort = name.sort)
    }
    file.name <- fcs.file.names
    runtype <- "OneTimepoint"
    output.folder <- MakeOutFolder(runtype = runtype)
    setwd(output.folder)
    fcs.file <- LoadCleanFCS(fcs.file.names = fcs.file.names, channel.remove = var.remove,
                             channel.annotate = var.annotate, subsamples = subsamples)
    file.clusters <- ClusterFCS(fcs.files = fcs.file, clustering.var = clustering.var,
                                numcluster = cluster.numbers, distance.metric = distance.metric)
    first.results <- BuildFirstFLOWMAP(FLOWMAP.clusters = file.clusters,
                                       per = per, min = minimum, max = maximum,
                                       distance.metric = distance.metric,
                                       clustering.var = clustering.var)
    output.graph <- first.results$output.graph
    output.graph <- AnnotateGraph(output.graph = output.graph,
                                  FLOWMAP.clusters = file.clusters,
                                  cellnum = subsamples)
    graph <- output.graph
  } else {
    stop("Unknown mode!")
  }
  file.name <- paste(unlist(strsplit(basename(file.name), "\\."))[1], "FLOW-MAP", sep = "_")
  ConvertToGraphML(output.graph = graph, file.name = file.name)
  graph.xy <- ForceDirectedXY(graph = graph)
  file.name.xy <- paste(file.name, "xy", sep = "_")
  final.file.name <- ConvertToGraphML(output.graph = graph.xy, file.name = file.name.xy)
  fixed.file.name <- paste(file.name.xy, "orig_time", sep = "_")
  graph.with.fixed.times <- ConvertOrigTime(graph.xy, orig.times)
  fixed.file <- ConvertToGraphML(output.graph = graph.with.fixed.times, file.name = fixed.file.name)
  PrintSummary(mode, files, var.annotate, var.remove,
               clustering.var, distance.metric, per, minimum,
               maximum, subsamples, cluster.numbers, seed.X)
  MakeFLOWMAPRFile(...)
  if (savePDFs) {
    cat("Printing pdfs.", "\n")
    ConvertToPDF(graphml.file = final.file.name, which.palette = which.palette, orig.times = orig.times)
  }
  return(graph.xy)
}
