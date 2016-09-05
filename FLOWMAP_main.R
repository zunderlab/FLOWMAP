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

SingleFLOWMAP <- function(folder, file.format, var.remove, var.annotate,
                          clustering.var, cluster.number, subsample, distance.metric,
                          minimum, maximum, per, save.folder, shuffle = FALSE) {
  fcs.file.names <- GetFCSNames(folder = folder, file.format = file.format)
  setwd(save.folder)
  output.folder <- MakeOutFolder(runtype = "singleFLOWMAP")
  setwd(output.folder)
  keep.folder <- getwd()
  print(keep.folder)
  fcs.files <- LoadCleanFCS(fcs.file.names = fcs.file.names, channel.remove = var.remove,
                            channel.annotate = var.annotate, subsample = subsample, subsample.rand = TRUE)
  # print("colnames(fcs.files[[1]])")
  # print(colnames(fcs.files[[1]]))
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
  # edgelist.save <- results$edgelist.save
  # total.edges <- 0
  # for (i in 1:length(edgelist.save)) {
  #   num.edges <- dim(edgelist.save[[i]])[1]
  #   total.edges <- total.edges + num.edges
  # }
  file.name <- paste(basename(folder), "original_edge_choice", sep = "_")
  ConvertToGraphML(output.graph = graph, file.name = file.name)
  graph.xy <- ForceDirectedXY(graph = graph)
  file.name.xy <- paste(basename(folder), "original_edge_choice", "xy", sep = "_")
  final.file.name <- ConvertToGraphML(output.graph = graph.xy, file.name = file.name.xy)
  ConvertToPDF(graphml.file = final.file.name, edge.color = "#FF000000")
  print(getwd())
  setwd(keep.folder)
  PrintSummary()
  return(graph.xy)
}


MultiFLOWMAP <- function(folder, file.format, var.remove,
                         var.annotate, clustering.var, cluster.number,
                         subsample, distance.metric, minimum,
                         maximum, per, save.folder, subsampleRand = TRUE,
                         shuffle = TRUE) {
  fcs.file.names <- GetMultiFCSNames(folder, file.format)
  setwd(save.folder)
  output.folder <- MakeOutFolder(runtype = "multiFLOWMAP")
  setwd(output.folder)
  keep.folder <- getwd()
  print(keep.folder)
  fcs.files <- LoadMultiCleanFCS(fcs.file.names, var.remove, var.annotate,
                                 subsample = subsample, subsample.rand)
  which.add.noise <- 1
  for (i in 1:length(fcs.files[[which.add.noise]])) {
    tmp <- fcs.files[[which.add.noise]][[i]]
    tmp <- AddNoise(tmp, factor = 1, amount = 0)
    fcs.files[[which.add.noise]][[i]] <- tmp
    rm(tmp)
  }
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
  file.name <- paste(basename(folder), "original_edge_choice", sep = "_")
  ConvertToGraphML(graph, file.name)
  graph.xy <- ForceDirectedXY(graph)
  file.name.xy <- paste(basename(folder), "original_edge_choice", "xy", sep = "_")
  final.file.name <- ConvertToGraphML(graph.xy, file.name.xy)
  ConvertToPDF(final.file.name, edge.color = "#FF000000")
  # visibleChannels <- c("Timepoint", "Treatment")
  # ConvertToPDF(in_folder, file_pattern, listOfTreatments = listOfTreatments) 
  # convertToPDF(in_folder, file_pattern, listOfTreatments = listOfTreatments,
  #             treatInvisible = TRUE, timeInvisible = TRUE, visibleChannels = visibleChannels) 
  print(getwd())
  setwd(keep.folder)
  PrintSummary()
}


AddNoise <- function(fcs.file, factor, amount) {
  if (!is.na(match("Treat", colnames(fcs.file)))) {
    tmp <- fcs.file[, -match("Treat", colnames(fcs.file))]
    Treat <- fcs.file[, "Treat"]
  }
  tmp <- as.matrix(tmp)
  tmp <- jitter(tmp, factor = factor, amount = amount)
  tmp <- as.data.frame(tmp)
  if (!is.na(match("Treat", colnames(fcs.file)))) {
    tmp <- cbind(tmp, Treat)
  }
  return(tmp)
}

