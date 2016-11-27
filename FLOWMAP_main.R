library(igraph)
library(Rclusterpp)
library(flowCore)
library(robustbase)
library(cluster)
library(fastcluster)
library(Biobase)
library(huge)
library(proxy)
library(scaffold)
library(spade)
library(Rclusterpp)

source(paste(prefolder, "FLOWMAP_loadFromFCS.R", sep = ""))
source(paste(prefolder, "FLOWMAP_clusterCells.R", sep = ""))
source(paste(prefolder, "FLOWMAP_buildGraph.R", sep = ""))
source(paste(prefolder, "FLOWMAP_forceDirected.R", sep = ""))
source(paste(prefolder, "FLOWMAP_saveGraphs.R", sep = ""))
source(paste(prefolder, "FLOWMAP_buildMultiGraph.R", sep = ""))


FLOWMAP <- function(files, file.format, var.remove, var.annotate,
                    clustering.var, cluster.number, subsample, distance.metric,
                    minimum, maximum, per, save.folder, shuffle = FALSE,
                    name.sort = TRUE) {
  # "files" variable could be one of the following:
  # a single fcs file path
  # a single folder path containing 2+ fcs files
  # a vector of fcs file paths
  # a single folder path, containing subfolders,
  # which each contain 2+ fcs files
  # a list name by treatment/conditions, each element is a 
  # vector of 2+ fcs file paths
  num.files <- 0
  single.flag <- length(files) <= 1
  fcs.file.flag <- length(grep(pattern = file.format, x = files)) > 0
  vector.flag <- is.vector(files)
  list.flag <- is.list(files)
  cat("single.flag is", single.flag, "\n")
  cat("fcs.file.flag is", fcs.file.flag, "\n")
  cat("vector.flag is", vector.flag, "\n")
  cat("list.flag is", list.flag, "\n")
  
  if (!fcs.file.flag & single.flag) {
    temp.files <- list.files(files, full.names = TRUE)
    temp.file <- temp.files[1]
    deeper.fcs.file.flag <- length(grep(pattern = file.format, x = temp.file)) > 0
    cat("deeper.fcs.file.flag is", deeper.fcs.file.flag, "\n")
    if (!deeper.fcs.file.flag) {
      temp.folder <- temp.file
      print("temp.folder")
      print(temp.folder)
      temp.files <- list.files(temp.folder, full.names = TRUE)
      temp.file <- temp.files[1]
      print("temp.file")
      print(temp.file)
      deepest.fcs.file.flag <- length(grep(pattern = file.format, x = temp.file)) > 0
      cat("deepest.fcs.file.flag is", deepest.fcs.file.flag, "\n")
      if (deepest.fcs.file.flag) {
        for (i in (1:length(temp.files))) {
          temp.folder <- temp.files[i]
          num.files <- num.files + length(list.files(temp.folder, pattern = file.format))
        }
        cat("number of FCS files is", num.files, "\n")
      }
    } else {
      num.files <- length(list.files(files, pattern = file.format))
      cat("number of FCS files is", num.files, "\n")
    }
  } else if (!fcs.file.flag & list.flag) {
    temp.files <- files[[1]]
    temp.file <- temp.files[1]
    deeper.fcs.file.flag <- length(grep(pattern = file.format, x = temp.file)) > 0
    cat("deeper.fcs.file.flag is", deeper.fcs.file.flag, "\n")
    if (deeper.fcs.file.flag) {
      for (i in length(files)) {
        num.files <- num.files + length(files[[i]])
      }
      cat("number of FCS files is", num.files, "\n")
    }
  }
  
  if (single.flag) {
    if (fcs.file.flag) {
      fcs.file.names <- files # only one FCS file given
      runtype <- "SingleTimepoint"
      num.files <- 1
      cat("number of FCS files is", num.files, "\n")
    } else if (deeper.fcs.file.flag) {
      # folder containing FCS files provided, load FCS file paths
      fcs.file.names <- GetFCSNames(folder = files, file.format = file.format, sort = name.sort)
      runtype <- "SingleFLOWMAP"
      num.files <- length(fcs.file.names)
      cat("number of FCS files is", num.files, "\n")
    } else if (deepest.fcs.file.flag) {
      # folder containing subfolders which contain 
      # FCS files provided, load FCS file paths
      fcs.file.names <- GetMultiFCSNames(folder = files, file.format = file.format, sort = name.sort)
      runtype <- "MultiFLOWMAP"
      num.files <- 0
      for (i in 1:length(fcs.file.names)) {
        num.files <- num.files + length(fcs.file.names[[i]])
      }
      cat("number of FCS files is", num.files, "\n")
    } else {
      stop("Unknown 'files' variable type provided!")
    }
  } else if (vector.flag & fcs.file.flag) {
    fcs.file.names <- files # FCS file paths provided by user
    num.files <- length(files)
    runtype <- "SingleFLOWMAP"
    cat("number of FCS files is", num.files, "\n")
  } else if (list.flag & deeper.fcs.file.flag) {
    # a list of FCS file paths provided by user for MultiFLOWMAP 
    fcs.file.names <- files # FCS file paths provided by user
    for (i in 1:length(fcs.file.names)) {
      num.files <- num.files + length(fcs.file.names[[i]])
    }
    runtype <- "MultiFLOWMAP"
    cat("number of FCS files is", num.files, "\n")
  } else {
    stop("Unknown 'files' variable type provided!")
  }
  
  cat("runtype is", runtype, "\n")
  
  setwd(save.folder)
  output.folder <- MakeOutFolder(runtype = runtype)
  setwd(output.folder)
  
  # num.files <- 
  cat("number of FCS files is", num.files, "\n")
  if (length(subsample) > 1 & length(subsample) != num.files & subsample != FALSE) {
    stop("Number to subsample not specified for all files!")
  } else if (length(subsample) == 1) {
    cat("Subsampling all files to:", subsample, "\n")
    # subsample.new <- rep(subsample, times = num.files)
    # subsample <- subsample.new
  } 
  # if not the above, then length(subsample) = length(num.files)
  
  if (runtype == "SingleFLOWMAP") {
    fcs.files <- LoadCleanFCS(fcs.file.names = fcs.file.names, channel.remove = var.remove,
                              channel.annotate = var.annotate, subsample = subsample, subsample.rand = TRUE)
    if (shuffle) {
      for (i in 1:length(fcs.files)) {
        df1 <- fcs.files[[i]]
        df2 <- df1[sample(nrow(df1)), ]
        fcs.files[[i]] <- df2
        rownames(fcs.files[[i]]) <- seq(1:subsample)
      }
    }
    if (cluster.number <= 0 || cluster.number == FALSE) {
      all.cells <- data.frame()
      for (i in 1:length(fcs.files)) {
        all.cells <- rbind(all.cells, fcs.files[[i]])
      }
      inf.flag <- sum(duplicated(all.cells))
      cat("inf.flag is", inf.flag, "\n")
      if (inf.flag != 0) {
        #
        # This method doesn't work perfectly yet.
        # Doesn't work for cells in different time points
        # that are identical to each other because
        # they can't be clustered together.
        #
        cat("Cannot use all cells as individual nodes due to duplicate values, clustering to eliminate these", "\n")
        cluster.number <- nrow(fcs.files[[1]]) - inf.flag
        print("cluster.number")
        print(cluster.number)
        file.clusters <- ClusterFCS(fcs.files = fcs.files, clustering.var = clustering.var,
                                    numcluster = cluster.number, distance.metric = distance.metric)
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
                                  numcluster = cluster.number, distance.metric = distance.metric)
    }
    results <- BuildFLOWMAP(FLOWMAP.clusters = file.clusters, per = per, min = minimum,
                            max = maximum, distance.metric = distance.metric, cellnum = subsample,
                            clustering.var = clustering.var)
    graph <- results$output.graph
    # print("graph")
    # print(graph)
    # print("name attribute for all vertices in graph")
    # print(get.vertex.attribute(graph, "name", index = V(graph)))
    # print("class of name attribute for all vertices in graph")
    # print(class(get.vertex.attribute(graph, "name", index = V(graph))))
  } else if (runtype == "MultiFLOWMAP") {
    fcs.files <- LoadMultiCleanFCS(fcs.file.names, var.remove, var.annotate,
                                   subsample = subsample, subsample.rand)
    if (shuffle) {
      for (n in 1:length(fcs.files)) {
        for (i in 1:length(fcs.files[[n]])) {
          df1 <- fcs.files[[n]][[i]]
          df2 <- df1[sample(nrow(df1)), ]
          fcs.files[[n]][[i]] <- df2
          rownames(fcs.files[[n]][[i]]) <- seq(1:subsample)
        }
      }
    }
    fcs.files.conversion <- ConvertNumericLabel(fcs.files)
    fixed.fcs.files <- fcs.files.conversion$fixed.list.FCS.files
    label.key <- fcs.files.conversion$label.key
    if (cluster.number <= 0 || cluster.number == FALSE) {
      stop("Not implemented yet!")
      # full.clusters
      # table.breaks
      # table.lengths
      # cluster.medians
      # cluster.counts
      # cell.assgn
      # file.clusters <- FLOWMAPcluster(full.clusters, table.breaks, table.lengths,
      #                                 cluster.medians, cluster.counts, cell.assgn)
    } else {
      file.clusters <- MultiClusterFCS(fixed.fcs.files, clustering.var = clustering.var, numcluster = cluster.number,
                                       distance.metric = distance.metric)
    }
    graph <- BuildMultiFLOWMAP(file.clusters, per = per, min = minimum,
                               max = maximum, distance.metric = distance.metric, cellnum = subsample)
  } else if (runtype == "SingleTimepoint") {
    fcs.file <- LoadCleanFCS(fcs.file.names = fcs.file.names, channel.remove = var.remove,
                             channel.annotate = var.annotate, subsample = subsample, subsample.rand = TRUE)
    if (cluster.number <= 0 || cluster.number == FALSE) {
      stop("Not implemented yet!")
      # full.clusters
      # table.breaks
      # table.lengths
      # cluster.medians
      # cluster.counts
      # cell.assgn
      file.clusters <- FLOWMAPcluster(full.clusters, table.breaks, table.lengths,
                                      cluster.medians, cluster.counts, cell.assgn)
    } else {
      file.clusters <- ClusterFCS(fcs.files = fcs.file, clustering.var = clustering.var,
                                  numcluster = cluster.number, distance.metric = distance.metric)
    }
    first.results <- BuildFirstFLOWMAP(FLOWMAP.clusters = file.clusters,
                                       per = per, min = minimum, max = maximum,
                                       distance.metric = distance.metric,
                                       clustering.var = clustering.var)
    output.graph <- first.results$output.graph
    output.graph <- AnnotateGraph(output.graph = output.graph,
                                  FLOWMAP.clusters = file.clusters,
                                  cellnum = subsample)
    graph <- output.graph
  }
  
  file.name <- paste(basename(folder), "FLOW-MAP", sep = "_")
  # print("graph")
  # print(graph)
  # print("edge weight attribute")
  # print(get.edge.attribute(graph, "weight", index = E(graph)))
  # problem <- which(get.edge.attribute(graph, "weight", index = E(graph)) == Inf)
  # print(problem)
  # print(E(graph)[problem])
  # v1 <- get.edgelist(graph)[problem, ][1]
  # v2 <- get.edgelist(graph)[problem, ][2]
  # print(get.vertex.attribute(graph, "marker1", index = v1))
  # print(get.vertex.attribute(graph, "marker2", index = v1))
  # print(get.vertex.attribute(graph, "marker1", index = v2))
  # print(get.vertex.attribute(graph, "marker2", index = v2))
  ConvertToGraphML(output.graph = graph, file.name = file.name)
  graph.xy <- ForceDirectedXY(graph = graph)
  print("graph.xy")
  print(graph.xy)
  file.name.xy <- paste(basename(folder), "FLOW-MAP", "xy", sep = "_")
  final.file.name <- ConvertToGraphML(output.graph = graph.xy, file.name = file.name.xy)
  PrintSummary()
  ConvertToPDF(graphml.file = final.file.name, edge.color = "#FF000000")
  return(graph.xy)
}






# SingleFLOWMAP <- function(files, file.format, var.remove, var.annotate,
#                           clustering.var, cluster.number, subsample, distance.metric,
#                           minimum, maximum, per, save.folder, shuffle = FALSE,
#                           name.sort = TRUE) {
#   # SingleFLOWMAP <- function(folder, file.format, var.remove, var.annotate,
#   if (length(files) <= 1) {
#     print("boop")
#     if (length(grep(pattern = file.format, x = files)) > 0) {
#       print("beep")
#       fcs.file.names <- files # only one FCS file given
#     } else {
#       print("bop")
#       # folder containing FCS files provided, load FCS file paths
#       fcs.file.names <- GetFCSNames(folder = files, file.format = file.format, sort = name.sort)
#     }
#   } else {
#     print("blah")
#     fcs.file.names <- files # FCS file paths provided by user
#   }
#   setwd(save.folder)
#   output.folder <- MakeOutFolder(runtype = "singleFLOWMAP")
#   setwd(output.folder)
#   fcs.files <- LoadCleanFCS(fcs.file.names = fcs.file.names, channel.remove = var.remove,
#                             channel.annotate = var.annotate, subsample = subsample, subsample.rand = TRUE)
#   if (shuffle) {
#     for (i in 1:length(fcs.files)) {
#       temp.fcs.file1 <- fcs.files[[i]]
#       temp.fcs.file2 <- df1[sample(nrow(temp.fcs.file1)), ]
#       fcs.files[[i]] <- temp.fcs.file2
#       rownames(fcs.files[[i]]) <- seq(1:subsample)
#       rm(temp.fcs.file1, temp.fcs.file2)
#     }
#   }
#   file.clusters <- ClusterFCS(fcs.files = fcs.files, clustering.var = clustering.var,
#                               numcluster = cluster.number, distance.metric = distance.metric)
#   results <- BuildFLOWMAP(FLOWMAP.clusters = file.clusters, per = per, min = minimum,
#                           max = maximum, distance.metric = distance.metric, cellnum = subsample,
#                           clustering.var = clustering.var)
#   graph <- results$output.graph
#   file.name <- paste(basename(folder), "FLOW-MAP", sep = "_")
#   ConvertToGraphML(output.graph = graph, file.name = file.name)
#   graph.xy <- ForceDirectedXY(graph = graph)
#   file.name.xy <- paste(basename(folder), "FLOW-MAP", "xy", sep = "_")
#   final.file.name <- ConvertToGraphML(output.graph = graph.xy, file.name = file.name.xy)
#   PrintSummary()
#   ConvertToPDF(graphml.file = final.file.name, edge.color = "#FF000000")
#   return(graph.xy)
# }
# 
# 
# MultiFLOWMAP <- function(folder, file.format, var.remove,
#                          var.annotate, clustering.var, cluster.number,
#                          subsample, distance.metric, minimum,
#                          maximum, per, save.folder, subsampleRand = TRUE,
#                          shuffle = TRUE) {
#   fcs.file.names <- GetMultiFCSNames(folder, file.format)
#   setwd(save.folder)
#   output.folder <- MakeOutFolder(runtype = "multiFLOWMAP")
#   setwd(output.folder)
#   fcs.files <- LoadMultiCleanFCS(fcs.file.names, var.remove, var.annotate,
#                                  subsample = subsample, subsample.rand)
#   # which.add.noise <- 1
#   # for (i in 1:length(fcs.files)) {
#   #   tmp <- fcs.files[[i]][[which.add.noise]]
#   #   tmp <- AddNoise(tmp, factor = 1, amount = 0)
#   #   fcs.files[[i]][[which.add.noise]] <- tmp
#   #   rm(tmp)
#   # }
#   if (shuffle) {
#     for (n in 1:length(fcs.files)) {
#       for (i in 1:length(fcs.files[[n]])) {
#         df1 <- fcs.files[[n]][[i]]
#         df2 <- df1[sample(nrow(df1)), ]
#         fcs.files[[n]][[i]] <- df2
#         rownames(fcs.files[[n]][[i]]) <- seq(1:subsample)
#       }
#     }
#   }
#   fcs.files.conversion <- ConvertNumericLabel(fcs.files)
#   fixed.fcs.files <- fcs.files.conversion$fixed.list.FCS.files
#   label.key <- fcs.files.conversion$label.key
#   file.clusters <- MultiClusterFCS(fixed.fcs.files, clustering.var = clustering.var, numcluster = cluster.number,
#                                    distance.metric = distance.metric)
#   graph <- BuildMultiFLOWMAP(file.clusters, per = per, min = minimum,
#                              max = maximum, distance.metric = distance.metric, cellnum = subsample)
#   file.name <- paste(basename(folder), "FLOW-MAP", sep = "_")
#   ConvertToGraphML(graph, file.name)
#   graph.xy <- ForceDirectedXY(graph)
#   file.name.xy <- paste(basename(folder), "FLOW-MAP", "xy", sep = "_")
#   final.file.name <- ConvertToGraphML(graph.xy, file.name.xy)
#   PrintSummary()
#   ConvertToPDF(final.file.name, edge.color = "#FF000000")
# }
# 
# 
# AddNoise <- function(fcs.file, factor, amount) {
#   if (!is.na(match("Treat", colnames(fcs.file)))) {
#     tmp <- fcs.file[, -match("Treat", colnames(fcs.file))]
#     Treat <- fcs.file[, "Treat"]
#   }
#   tmp <- as.matrix(tmp)
#   tmp <- jitter(tmp, factor = factor, amount = amount)
#   tmp <- as.data.frame(tmp)
#   if (!is.na(match("Treat", colnames(fcs.file)))) {
#     tmp <- cbind(tmp, Treat)
#   }
#   return(tmp)
# }

