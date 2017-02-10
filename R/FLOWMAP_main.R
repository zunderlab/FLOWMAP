#' FLOWMAPR
#' @name FLOWMAPR
#' @docType package
#' @import flowCore
#' @import Rclusterpp
#' @import robustbase
#' @import igraph
#' @import scaffold
#' @import SDMTools

DefineRuntype <- function(files, file.format, name.sort) {
  num.files <- 0
  single.flag <- length(files) <= 1
  fcs.file.flag <- length(grep(pattern = file.format, x = files)) > 0
  vector.flag <- is.vector(files)
  list.flag <- is.list(files)
  if (!fcs.file.flag & single.flag) {
    temp.files <- list.files(files, full.names = TRUE)
    temp.file <- temp.files[1]
    deeper.fcs.file.flag <- length(grep(pattern = file.format, x = temp.file)) > 0
    if (!deeper.fcs.file.flag) {
      temp.folder <- temp.file
      temp.files <- list.files(temp.folder, full.names = TRUE)
      temp.file <- temp.files[1]
      deepest.fcs.file.flag <- length(grep(pattern = file.format, x = temp.file)) > 0
      if (deepest.fcs.file.flag) {
        for (i in (1:length(temp.files))) {
          temp.folder <- temp.files[i]
          num.files <- num.files + length(list.files(temp.folder, pattern = file.format))
        }
      }
    } else {
      num.files <- length(list.files(files, pattern = file.format))
    }
  } else if (!fcs.file.flag & list.flag) {
    temp.files <- files[[1]]
    temp.file <- temp.files[1]
    deeper.fcs.file.flag <- length(grep(pattern = file.format, x = temp.file)) > 0
    if (deeper.fcs.file.flag) {
      for (i in length(files)) {
        num.files <- num.files + length(files[[i]])
      }
    }
  }
  if (single.flag) {
    if (fcs.file.flag) {
      fcs.file.names <- files # only one FCS file given
      runtype <- "SingleTimepoint"
      num.files <- 1
    } else if (deeper.fcs.file.flag) {
      fcs.file.names <- GetFCSNames(folder = files, file.format = file.format, sort = name.sort)
      runtype <- "SingleFLOWMAP"
      num.files <- length(fcs.file.names)
    } else if (deepest.fcs.file.flag) {
      # folder containing subfolders which contain 
      # FCS files provided, load FCS file paths
      fcs.file.names <- GetMultiFCSNames(folder = files, file.format = file.format, sort = name.sort)
      runtype <- "MultiFLOWMAP"
      num.files <- 0
      for (i in 1:length(fcs.file.names)) {
        num.files <- num.files + length(fcs.file.names[[i]])
      }
    } else {
      stop("Unknown 'files' variable type provided!")
    }
  } else if (vector.flag & fcs.file.flag) {
    fcs.file.names <- files # FCS file paths provided by user
    num.files <- length(files)
    runtype <- "SingleFLOWMAP"
  } else if (list.flag & deeper.fcs.file.flag) {
    # a list of FCS file paths provided by user for MultiFLOWMAP 
    fcs.file.names <- files # FCS file paths provided by user
    for (i in 1:length(fcs.file.names)) {
      num.files <- num.files + length(fcs.file.names[[i]])
    }
    runtype <- "MultiFLOWMAP"
  } else {
    stop("Unknown 'files' variable type provided!")
  }
  cat("runtype is", runtype, "\n")
  return(list(runtype = runtype,
              num.files = num.files,
              fcs.file.names = fcs.file.names))
}


#' @export
FLOWMAP <- function(files, file.format, var.remove, var.annotate,
                    clustering.var, cluster.numbers, subsamples, distance.metric,
                    minimum, maximum, per, save.folder, shuffle = FALSE,
                    name.sort = TRUE, downsample = TRUE) {
  # "files" variable could be one of the following:
  # a single fcs file path
  # a single folder path containing 2+ fcs files
  # a vector of fcs file paths
  # a single folder path, containing subfolders,
  # which each contain 2+ fcs files
  # a list name by treatment/conditions, each element is a 
  # vector of 2+ fcs file paths
  runtype.results <- DefineRuntype(files, file.format, name.sort)
  runtype <- runtype.results$runtype
  num.files <- runtype.results$num.files
  fcs.file.names <- runtype.results$fcs.file.names
  setwd(save.folder)
  output.folder <- MakeOutFolder(runtype = runtype)
  setwd(output.folder)
  
  if (downsample) {
    cat("Downsampling all files using SPADE functions")
    fcs.files <- DownsampleFCS(fcs.file.names, clustering.var,
                               distance.metric, exclude.pctile = 0.01,
                               target.pctile = NULL,
                               target.number = NULL,
                               target.percent = 0.1)
  }
  
  if (length(subsamples) > 1 & length(subsamples) != num.files & subsamples != FALSE) {
    stop("Number to subsample not specified for all files!")
  } else if (length(subsamples) == 1) {
    cat("Subsampling all files to:", subsamples, "\n")
  } 
  # if not the above, then length(subsample) = length(num.files)
  # NOTE(Jordan): Separate function for each different runtype.
  # NOTE(Jordan): Can SingleFLOWMAP be a special case of MultiFLOWMAP? Some of the code is the same.
  if (runtype == "SingleFLOWMAP") {
    fcs.files <- LoadCleanFCS(fcs.file.names = fcs.file.names, channel.remove = var.remove,
                              channel.annotate = var.annotate, subsamples = subsamples, subsample.rand = TRUE)
    if (shuffle) {
      # NOTE(Jordan): Helper function.
      for (i in 1:length(fcs.files)) {
        if (length(subsamples) > 1) {
          subsamp <- subsamples[i]
        } else {
          subsamp <- subsamples
        }
        df1 <- fcs.files[[i]]
        df2 <- df1[sample(nrow(df1)), ]
        fcs.files[[i]] <- df2
        rownames(fcs.files[[i]]) <- seq(1:subsamp)
      }
    }
    if (cluster.numbers <= 0 || cluster.numbers == FALSE) {
      all.cells <- data.frame()
      for (i in 1:length(fcs.files)) {
        all.cells <- rbind(all.cells, fcs.files[[i]])
      }
      inf.flag <- sum(duplicated(all.cells))
      cat("inf.flag is", inf.flag, "\n")
      if (inf.flag != 0) {
        # This method doesn't work perfectly yet.
        # Doesn't work for cells in different time points
        # that are identical to each other because
        # they can't be clustered together.
        cat("Cannot use all cells as individual nodes due to duplicate values, clustering to eliminate these", "\n")
        cluster.numbers <- nrow(fcs.files[[1]]) - inf.flag
        file.clusters <- ClusterFCS(fcs.files = fcs.files, clustering.var = clustering.var,
                                    numcluster = cluster.numbers, distance.metric = distance.metric)
      } else {
        cat("Use all cells as individual nodes, no clustering", "\n")
        full.clusters <- data.frame()
        table.breaks <- c()
        table.lengths <- c()
        cluster.medians <- list()
        cluster.counts <- list()
        cell.assgn <- list()
        for (i in 1:length(fcs.files)) {
          full.clusters <- rbind(full.clusters, fcs.files[[i]])
          cluster.medians[[i]] <- fcs.files[[i]]
          cluster.counts[[i]] <- as.data.frame(as.integer(rep(1, times = nrow(fcs.files[[i]]))))
          colnames(cluster.counts[[i]]) <- c("Counts")
          table.lengths <- append(table.lengths, nrow(cluster.medians[[i]]))
          table.breaks <- append(table.breaks, nrow(full.clusters))
          cell.assgn[[i]] <- as.data.frame(seq(1, nrow(cluster.medians[[i]])))
          colnames(cell.assgn[[i]]) <- c("Cluster")
        }
        for (i in 1:length(cluster.medians)) {
          rownames(cluster.medians[[i]]) <- seq(1, table.lengths[i])
        }
        rownames(full.clusters) <- seq(1, dim(full.clusters)[1])
        file.clusters <- FLOWMAPcluster(full.clusters, table.breaks, table.lengths,
                                        cluster.medians, cluster.counts, cell.assgn)
      }
    } else {
      file.clusters <- ClusterFCS(fcs.files = fcs.files, clustering.var = clustering.var,
                                  numcluster = cluster.numbers, distance.metric = distance.metric)
    }
    results <- BuildFLOWMAP(FLOWMAP.clusters = file.clusters, per = per, min = minimum,
                            max = maximum, distance.metric = distance.metric, cellnum = subsamples,
                            clustering.var = clustering.var)
    graph <- results$output.graph
  } else if (runtype == "MultiFLOWMAP") {
    fcs.files <- LoadMultiCleanFCS(fcs.file.names, var.remove, var.annotate,
                                   subsamples = subsamples, subsample.rand)
    if (shuffle) {
      for (n in 1:length(fcs.files)) {
        for (i in 1:length(fcs.files[[n]])) {
          if (length(subsamples) > 1) {
            subsamp <- subsamples[[n]][i]
          } else {
            subsamp <- subsamples
          }
          df1 <- fcs.files[[n]][[i]]
          df2 <- df1[sample(nrow(df1)), ]
          fcs.files[[n]][[i]] <- df2
          rownames(fcs.files[[n]][[i]]) <- seq(1:subsamp)
        }
      }
    }
    fcs.files.conversion <- ConvertNumericLabel(fcs.files)
    fixed.fcs.files <- fcs.files.conversion$fixed.list.FCS.files
    label.key <- fcs.files.conversion$label.key
    
    if (cluster.numbers <= 0 || cluster.numbers == FALSE) {
      stop("Not implemented yet!")
      # file.clusters <- FLOWMAPcluster(full.clusters, table.breaks, table.lengths,
      #                                 cluster.medians, cluster.counts, cell.assgn)
    } else {
      file.clusters <- MultiClusterFCS(fixed.fcs.files, clustering.var = clustering.var, numcluster = cluster.numbers,
                                       distance.metric = distance.metric)
    }
    graph <- BuildMultiFLOWMAP(file.clusters, per = per, min = minimum,
                               max = maximum, distance.metric = distance.metric, cellnum = subsamples,
                               label.key = label.key)
  } else if (runtype == "SingleTimepoint") {
    fcs.file <- LoadCleanFCS(fcs.file.names = fcs.file.names, channel.remove = var.remove,
                             channel.annotate = var.annotate, subsamples = subsamples, subsample.rand = TRUE)
    if (cluster.numbers <= 0 || cluster.numbers == FALSE) {
      stop("Not implemented yet!")
      # file.clusters <- FLOWMAPcluster(full.clusters, table.breaks, table.lengths,
      #                                 cluster.medians, cluster.counts, cell.assgn)
    } else {
      file.clusters <- ClusterFCS(fcs.files = fcs.file, clustering.var = clustering.var,
                                  numcluster = cluster.numbers, distance.metric = distance.metric)
    }
    first.results <- BuildFirstFLOWMAP(FLOWMAP.clusters = file.clusters,
                                       per = per, min = minimum, max = maximum,
                                       distance.metric = distance.metric,
                                       clustering.var = clustering.var)
    output.graph <- first.results$output.graph
    output.graph <- AnnotateGraph(output.graph = output.graph,
                                  FLOWMAP.clusters = file.clusters,
                                  cellnum = subsamples)
    graph <- output.graph
  }
  if (length(files) > 1) {
    file.name <- files[1]
  }
  file.name <- paste(unlist(strsplit(basename(files), "\\."))[1], "FLOW-MAP", sep = "_")
  ConvertToGraphML(output.graph = graph, file.name = file.name)
  graph.xy <- ForceDirectedXY(graph = graph)
  file.name.xy <- paste(file.name, "xy", sep = "_")
  final.file.name <- ConvertToGraphML(output.graph = graph.xy, file.name = file.name.xy)
  PrintSummary()
  ConvertToPDF(graphml.file = final.file.name, edge.color = "#FF000000")
  return(graph.xy)
}
