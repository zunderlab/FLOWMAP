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
                    name.sort = TRUE, no.cluster = FALSE) {
  # "files" variable could be one of the following:
  # a single fcs file path
  # a single folder path containing 2+ fcs files
  # a vector of fcs file paths
  # a single folder path, containing subfolders,
  # which each contain 2+ fcs files
  # a list name by treatment/conditions, each element is a 
  # vector of 2+ fcs file paths
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
      temp.files <- list.files(temp.folder)
      temp.file <- temp.files[1]
      deepest.fcs.file.flag <- length(grep(pattern = file.format, x = temp.file)) > 0
      cat("deepest.fcs.file.flag is", deepest.fcs.file.flag, "\n")
    }
  } else if (!fcs.file.flag & list.flag) {
    temp.files <- files[[1]]
    temp.file <- temp.files[1]
    deeper.fcs.file.flag <- length(grep(pattern = file.format, x = temp.file)) > 0
    cat("deeper.fcs.file.flag is", deeper.fcs.file.flag, "\n")
  }
  
  if (single.flag) {
    if (fcs.file.flag) {
      fcs.file.names <- files # only one FCS file given
      runtype <- "SingleTimepoint"
    } else if (deeper.fcs.file.flag) {
      # folder containing FCS files provided, load FCS file paths
      fcs.file.names <- GetFCSNames(folder = files, file.format = file.format, sort = name.sort)
      runtype <- "SingleFLOWMAP"
    } else if (deepest.fcs.file.flag) {
      # folder containing subfolders which contain 
      # FCS files provided, load FCS file paths
      fcs.file.names <- GetMultiFCSNames(folder = files, file.format = file.format, sort = name.sort)
      runtype <- "MultiFLOWMAP"
    } else {
      stop("Unknown 'files' variable type provided!")
    }
  } else if (list.flag & deeper.fcs.file.flag) {
    # a list of FCS file paths provided by user for MultiFLOWMAP 
    fcs.file.names <- files # FCS file paths provided by user
    runtype <- "MutliFLOWMAP"
  } else if (vector.flag & fcs.file.flag) {
    fcs.file.names <- files # FCS file paths provided by user
    runtype <- "SingleFLOWMAP"
  } else {
    stop("Unknown 'files' variable type provided!")
  }

  cat("runtype is", runtype, "\n")
  
  setwd(save.folder)
  output.folder <- MakeOutFolder(runtype = runtype)
  setwd(output.folder)
  
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
    file.clusters <- ClusterFCS(fcs.files = fcs.files, clustering.var = clustering.var,
                                numcluster = cluster.number, distance.metric = distance.metric)
    results <- BuildFLOWMAP(FLOWMAP.clusters = file.clusters, per = per, min = minimum,
                            max = maximum, distance.metric = distance.metric, cellnum = subsample,
                            clustering.var = clustering.var)
    graph <- results$output.graph
    
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
    file.clusters <- MultiClusterFCS(fixed.fcs.files, clustering.var = clustering.var, numcluster = cluster.number,
                                     distance.metric = distance.metric)
    graph <- BuildMultiFLOWMAP(file.clusters, per = per, min = minimum,
                               max = maximum, distance.metric = distance.metric, cellnum = subsample)
  } else if (runtype == "SingleTimepoint") {
    fcs.file <- LoadCleanFCS(fcs.file.names = fcs.file.names, channel.remove = var.remove,
                              channel.annotate = var.annotate, subsample = subsample, subsample.rand = TRUE)
    file.clusters <- ClusterFCS(fcs.files = fcs.file, clustering.var = clustering.var,
                                numcluster = cluster.number, distance.metric = distance.metric)
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
  ConvertToGraphML(output.graph = graph, file.name = file.name)
  graph.xy <- ForceDirectedXY(graph = graph)
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

