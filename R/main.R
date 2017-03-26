#' FLOWMAPR
#' @name FLOWMAPR
#' @docType package
#' @import flowCore
#' @import Rclusterpp
#' @import robustbase
#' @import igraph
#' @import scaffold
#' @import SDMTools

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
        subsamp <- subsamples[[n]][i]
      } else if (subsamples == FALSE) {
        subsamp <- nrow(fcs.files[[n]][i])
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

#' @export
FLOWMAP <- function(files, var.remove, var.annotate, clustering.var,
                    cluster.numbers, subsamples, distance.metric,
                    minimum, maximum, per, save.folder, mode = c("single", "multi", "one"),
                    starting.files = c("FCS", "cluster_matrix"),
                    shuffle = TRUE, name.sort = TRUE, downsample = TRUE, ...) {
  # "files" variable could be one of the following:
  # a single fcs file path
  # a single folder path containing 1+ fcs files
  # a vector of fcs file paths
  # a single folder path, containing subfolders,
  # which each contain 2+ fcs files
  # a list name by time, each element is a vector
  # of 1+ fcs file paths from each treatment/condition
  
  # optional variables
  # transform 
  # scale
  # subsample.rand
  # exclude.pctile
  # target.pctile 
  # target.number 
  # target.percent
  # starting.files = c("FCS", "cluster_matrix")
  
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
    file.name <- fcs.file.names[1]
    # file.name <- fcs.file.names[1]
    if (downsample) {
      cat("Downsampling all files using SPADE functions", "\n")
      fcs.file.names <- DownsampleFCS(fcs.file.names, clustering.var,
                                      distance.metric, ...)
    }
    fcs.files <- LoadCleanFCS(fcs.file.names = fcs.file.names, channel.remove = var.remove,
                              channel.annotate = var.annotate, subsamples = subsamples)
    if (shuffle) {
      fcs.files <- ShuffleCells(fcs.files, subsamples)
    }
    file.clusters <- ClusterFCS(fcs.files = fcs.files, clustering.var = clustering.var,
                                numcluster = cluster.numbers, distance.metric = distance.metric)
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
    }
    if (check[2] == "subfolder") {
      fcs.file.names <- GetMultiFCSNames(folder = files, sort = name.sort)
    }
    file.name <- fcs.file.names[[1]][1]
    if (downsample) {
      cat("Downsampling all files using SPADE functions", "\n")
      fcs.file.names <- DownsampleFCS(fcs.file.names, clustering.var,
                                      distance.metric, ...)
    }
    fcs.files <- LoadMultiCleanFCS(fcs.file.names, var.remove, var.annotate,
                                   subsamples = subsamples)
    if (shuffle) {
      fcs.files <- MultiShuffleCells(fcs.files, subsamples)
    }
    fcs.files.conversion <- ConvertNumericLabel(fcs.files)
    fixed.fcs.files <- fcs.files.conversion$fixed.list.FCS.files
    label.key <- fcs.files.conversion$label.key
    file.clusters <- MultiClusterFCS(fixed.fcs.files, clustering.var = clustering.var, numcluster = cluster.numbers,
                                     distance.metric = distance.metric)
    graph <- BuildMultiFLOWMAP(file.clusters, per = per, min = minimum,
                               max = maximum, distance.metric = distance.metric, cellnum = subsamples,
                               label.key = label.key)
  } else if (mode == "one") {
    check <- CheckModeOne(files)
    cat("check", check, "\n")
    if (check[1]) {
      stop("Unknown 'files' format provided for specified mode!")
    }
    fcs.file.names <- files # only one FCS file given
    file.name <- fcs.file.names
    runtype <- "SingleTimepoint"
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
  PrintSummary()
  ConvertToPDF(graphml.file = final.file.name, edge.color = "#FF000000")
  return(graph.xy)
}
