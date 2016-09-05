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
  # print("clustering.var")
  # print(clustering.var)
  # print("setdiff(colnames(fcs.files[[1]]), clustering.var)")
  # print(setdiff(colnames(fcs.files[[1]]), clustering.var))
  # print("setdiff(clustering.var, colnames(fcs.files[[1]]))")
  # print(setdiff(clustering.var, colnames(fcs.files[[1]])))
  file.clusters <- ClusterFCS(fcs.files = fcs.files, clustering.var = clustering.var,
                              numcluster = cluster.number, distance.metric = distance.metric)
  results <- BuildFLOWMAP(FLOWMAP.clusters = file.clusters, per = per, min = minimum,
                        max = maximum, distance.metric = distance.metric, cellnum = subsample,
                        clustering.var = clustering.var)
  graph <- results$output.graph
  edgelist.save <- results$edgelist.save
  # print("edgelist.save")
  # print(edgelist.save)
  total.edges <- 0
  # cat("length(edgelist.save)", length(edgelist.save), "\n")
  for (i in 1:length(edgelist.save)) {
    # cat("i is", i, "\n")
    # cat("names(edgelist.save)[i] is", names(edgelist.save)[i], "\n")
    num.edges <- dim(edgelist.save[[i]])[1]
    total.edges <- total.edges + num.edges
  }
  # total.edges <- total.edges + dim(edgelist.save$first)[1]
  # print("total.edges")
  # print(total.edges)
  # graph <- BuildFLOWMAP(FLOWMAP.clusters = file.clusters, per = per, min = minimum,
  #                       max = maximum, distance.metric = distance.metric, cellnum = subsample,
  #                       clustering.var = clustering.var)
  file.name <- paste(basename(folder), "original_edge_choice", sep = "_")
  ConvertToGraphML(output.graph = graph, file.name = file.name)
  graph.xy <- ForceDirectedXY(graph = graph)
  file.name.xy <- paste(basename(folder), "original_edge_choice", "xy", sep = "_")
  final.file.name <- ConvertToGraphML(output.graph = graph.xy, file.name = file.name.xy)
  ConvertToPDF(graphml.file = final.file.name, edge.color = "#FF000000")
  print(getwd())
  setwd(keep.folder)
  printSummary()
  return(graph.xy)
}


MultiFLOWMAP <- function(folder, file.format, var.remove,
                         var.annotate, clustering.var, cluster.number, saveGRAPHML = TRUE,
                         savePDFS = TRUE, subsampleRand = TRUE) {
  fcs.file.names <- GetMultiFCSNames(folder, file.format)
  setwd(save.folder)
  output.folder <- MakeOutFolder(runtype = "multiFLOWMAP")
  setwd(output.folder)
  keep.folder <- getwd()
  print(keep.folder)
  fcs.files <- LoadMultiCleanFCS(fcs.file.names, var.remove, var.annotate,
                                 subsample = subsample, subsample.rand)
  fcs.files.conversion <- ConvertNumericLabel(fcs.files)
  fixed.fcs.files <- fcs.files.conversion$fixed.list.FCS.files
  label.key <- fcs.files.conversion$label.key
  file.clusters <- MultiClusterFCS(fixed.fcs.files, clustering.var = clustering.var, numcluster = cluster.number,
                                   distance.metric = distance.metric)
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
  setwd(keep.folder)
  printSummary()
}

