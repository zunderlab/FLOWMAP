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

source(paste(prefolder, "FLOWMAP_loadFromFCS.R", sep = ""))
source(paste(prefolder, "FLOWMAP_clusterCells.R", sep = ""))
source(paste(prefolder, "FLOWMAP_buildGraph.R", sep = ""))
source(paste(prefolder, "FLOWMAP_forceDirected.R", sep = ""))
source(paste(prefolder, "FLOWMAP_saveGraphs.R", sep = ""))

SingleFLOWMAP <- function(folder, file.format, var.remove, var.annotate,
                          clustering.var, cluster.number, subsample, distance.metric,
                          minimum, maximum, per, shuffle = FALSE) {
  fcs.file.names <- GetFCSNames(folder, file.format)
  output.folder <- MakeOutFolder(runtype = "singleFLOWMAP")
  setwd(output.folder)
  save.folder <- getwd()
  print(save.folder)
  fcs.files <- LoadCleanFCS(fcs.file.names, var.remove, var.annotate, subsample = subsample, subsample.rand = TRUE)
  if (shuffle) {
    for (i in 1:length(fcs.files)) {
      df1 <- fcs.files[[i]]
      df2 <- df1[sample(nrow(df1)), ]
      fcs.files[[i]] <- df2
    }
  }
  file.clusters <- ClusterFCS(fcs.files, clustering.var = clustering.var, numcluster = cluster.number,
                              distance.metric = distance.metric)
  graph <- BuildFLOWMAP(file.clusters, per = per, min = minimum,
                        max = maximum, distance.metric, cellnum = subsample,
                        clustering.var = clustering.var)
  file.name <- paste(basename(folder), "original_edge_choice", sep = "_")
  ConvertToGraphML(graph, file.name)
  graph.xy <- ForceDirectedXY(graph)
  file.name.xy <- paste(basename(folder), "original_edge_choice", "xy", sep = "_")
  final.file.name <- ConvertToGraphML(graph.xy, file.name.xy)
  ConvertToPDF(final.file.name, edge.color = "#FF000000")
  print(getwd())
  setwd(save.folder)
  printSummary()
}


MultiFLOWMAP <- function(folder, file.format, var.remove,
                         var.annotate, clustering.var, cluster.number, saveGRAPHML = TRUE,
                         savePDFS = TRUE, subsampleRand = TRUE) {
  fcs.file.names <- GetMultiFCSNames(folder, file.format)
  output.folder <- MakeOutFolder(runtype = "multiFLOWMAP")
  setwd(output.folder)
  save.folder <- getwd()
  print(save.folder)
  fcs.files <- LoadMultiCleanFCS(fcs.file.names, var.remove, var.annotate,
                                 subsample = subsample, subsample.rand)
  file.clusters <- MultiClusterFCS(fcs.files, channel.cluster = clustering.var, numcluster = cluster.number)
  graph <- BuildMultiFLOWMAP(file.clusters, per = per, min = minimum,
                             max = maximum, distance.metric = distance.metric, cellnum = subsample)
  file.name <- paste(basename(folder), "original_edge_choice", sep = "_")
  ConvertToGraphML(graph, file.name)
  graph.xy <- forceDirectedXY(graph)
  file.name.xy <- paste(basename(folder), "original_edge_choice", "xy", sep = "_")
  final.file.name <- ConvertToGraphML(graph.xy, file.name.xy)
  ConvertToPDF(final.file.name, edge.color = "#FF000000")
  # visibleChannels <- c("Timepoint", "Treatment")
  ConvertToPDF(in_folder, file_pattern, listOfTreatments = listOfTreatments) 
  # convertToPDF(in_folder, file_pattern, listOfTreatments = listOfTreatments,
  #             treatInvisible = TRUE, timeInvisible = TRUE, visibleChannels = visibleChannels) 
  print(getwd())
  setwd(save.folder)
  printSummary()
}

