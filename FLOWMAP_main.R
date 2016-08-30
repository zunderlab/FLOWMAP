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

singleFLOWMAP <- function(FOLDER, FILE_FORMAT, VAR_REMOVE, VAR_ANNOTATE,
                          CLUSTERING_VAR, CLUSTNUM, SUBSAMPLE, distance_metric,
                          minimum, maximum, per, shuffle = FALSE) {
  fcs_file_names <- getFCSNames(FOLDER, FILE_FORMAT)
  output_folder <- makeOutFolder(runtype = "singleFLOWMAP")
  setwd(output_folder)
  save_folder <- getwd()
  print(save_folder)
  fcs_files <- loadCleanFCS(fcs_file_names, VAR_REMOVE, VAR_ANNOTATE, subsample = SUBSAMPLE, subsampleRand = TRUE)
  if (shuffle) {
    for (i in 1:length(fcs_files)) {
    df1 <- fcs_files[[i]]
    df2 <- df1[sample(nrow(df1)), ]
    fcs_files[[i]] <- df2
    }
  }
  file_clusters <- clusterFCS(fcs_files, CLUSTERING_VAR, numcluster = CLUSTNUM,
                              distance_metric = distance_metric)
  graph <- buildFLOWMAP(file_clusters, per = per, min = minimum,
                        max = maximum, distance_metric, cellnum = SUBSAMPLE)
  file_name <- paste(basename(FOLDER), "original_edge_choice", sep = "_")
  convertToGraphML(graph, file_name)
  graph_xy <- forceDirectedXY(graph)
  file_name_xy <- paste(basename(FOLDER), "original_edge_choice", "xy", sep = "_")
  final_file_name <- convertToGraphML(graph_xy, file_name_xy)
  convertToPDF(final_file_name, edge_color = "#FF000000")
  print(getwd())
  setwd(save_folder)
  printSummary()
}


singleFLOWMAP1B <- function(FOLDER, FILE_FORMAT, VAR_REMOVE, VAR_ANNOTATE,
                            CLUSTERING_VAR, CLUSTNUM, saveGRAPHML = TRUE,
                            savePDFS = TRUE, subsampleRand = TRUE,
                            noEdge = FALSE, ...) {
  fcs_file_names <- getFCSNames(FOLDER, FILE_FORMAT)
  output_folder <- makeOutFolder(runtype = "singleFLOWMAP")
  fcs_files <- loadCleanFCS(fcs_file_names, VAR_REMOVE, VAR_ANNOTATE, subsample = SUBSAMPLE, subsampleRand = subsampleRand)
  file_clusters <- clusterFCS(fcs_files, CLUSTERING_VAR, numcluster = CLUSTNUM,
                              distance_metric = distance_metric)
  cellassgn <- file_clusters$cellassgn
  graph <- buildFLOWMAP(file_clusters, per = per, min = minimum,
                        max = maximum, distance_metric, cellnum = SUBSAMPLE)
  x <- "original_edge_choice"
  if (saveGRAPHML) {
    file_name <- paste(basename(FOLDER), x, sep = "_")
    convertToGraphML(graph, output_folder, file_name)
  }
  graph_xy <- forceDirectedXY(graph)
  file_name_xy <- paste(basename(FOLDER), x, "xy", sep = "_")
  convertToGraphML(graph_xy, output_folder, file_name_xy)
  if (savePDFS) {
    in_folder <- output_folder
    file_pattern <- file_name_xy
    if (noEdge) {
      convertToPDF(in_folder, file_pattern, edge_color = "#FF000000") 
    }
    else {
      convertToPDF(in_folder, file_pattern) 
    }
  }
  printSummary(output_folder)
}



singleFLOWMAP2 <- function(FOLDER, FILE_FORMAT, VAR_REMOVE, VAR_ANNOTATE,
                           CLUSTERING_VAR, CLUSTNUM, saveGRAPHML = TRUE,
                           savePDFS = TRUE, subsampleRand = TRUE,
                           findStages = TRUE, exportStages = TRUE,
                           noEdge = FALSE, ...) {
  fcs_file_names <- getFCSNames(FOLDER, FILE_FORMAT)
  output_folder <- makeOutFolder(runtype = "singleFLOWMAP")
  if (exists("VAR_UPSAMPLE")) {
    fcs_files <- loadCleanFCS(fcs_file_names, VAR_REMOVE, VAR_ANNOTATE)
    fcs_files <- upsampleFCS(fcs_files, VAR_UPSAMPLE, percent_upsample_limit = UPLIMIT, subsample = SUBSAMPLE) 
  }
  else {
    fcs_files <- loadCleanFCS(fcs_file_names, VAR_REMOVE, VAR_ANNOTATE, subsample = SUBSAMPLE, subsampleRand)
  }
  file_clusters <- clusterFCS(fcs_files, CLUSTERING_VAR, numcluster = CLUSTNUM)
  cellassgn <- file_clusters$cellassgn
  if (exists("top_npercent")) {
    cellnum <- SUBSAMPLE
    graph <- buildSimpleFLOWMAP(file_clusters, top_npercent, distance_metric, cellnum = cellnum)
    x <- "top_npercent_edge_choice"
  }
  else {
    cellnum <- SUBSAMPLE
    graph <- buildFLOWMAP(file_clusters, per = per, min = minimum,
                          max = maximum, distance_metric, cellnum = cellnum)
    x <- "original_edge_choice"
  }
  if (saveGRAPHML) {
    file_name <- paste(basename(FOLDER), x, sep = "_")
    convertToGraphML(graph, output_folder, file_name)
  }
  # graph_xy <- forceDirectedXY(graph)
  # file_name_xy <- paste(basename(FOLDER), x, "xy", sep = "_")
  # convertToGraphML(graph_xy, output_folder, file_name_xy)
  if (savePDFS) {
    in_folder <- output_folder
    file_pattern <- file_name_xy
    #     visibleChannels <- c("Timepoint", "Treatment")
    if (noEdge) {
      convertToPDF(in_folder, file_pattern, edge_color = "#FF000000") 
    }
    else {
      convertToPDF(in_folder, file_pattern) 
    }
    #     convertToPDF(in_folder, file_pattern, listOfTreatments = listOfTreatments,
    #                  treatInvisible = TRUE, timeInvisible = TRUE, visibleChannels = visibleChannels) 
  }
  printSummary(output_folder)
}


multiFLOWMAP <- function(listOfTreatments, FOLDER, FILE_FORMAT, VAR_REMOVE,
                         VAR_ANNOTATE, CLUSTERING_VAR, CLUSTNUM, saveGRAPHML = TRUE,
                         savePDFS = TRUE, subsampleRand = TRUE, ...) {
  fcs_file_names <- getMultiFCSNames(listOfTreatments, FOLDER, FILE_FORMAT)
  output_folder <- "~/Desktop"
  output_folder <- makeOutFolder(runtype = "multiFLOWMAP")
  if (exists("VAR_UPSAMPLE")) {
    fcs_files <- loadMultiCleanFCS(fcs_file_names, VAR_REMOVE, VAR_ANNOTATE)
    fcs_files <- upsampleMultiFCS(fcs_files, VAR_UPSAMPLE, percent_upsample_limit = UPLIMIT,
                                  subsample = SUBSAMPLE) 
  }
  else {
    fcs_files <- loadMultiCleanFCS(fcs_file_names, VAR_REMOVE, VAR_ANNOTATE,
                                   subsample = SUBSAMPLE, subsampleRand)
  }
  file_clusters <- multiClusterFCS(fcs_files, channel_cluster = CLUSTERING_VAR, numcluster = CLUSTNUM)
  if (exists("top_npercent")) {
    cellnum <- SUBSAMPLE
    graph <- buildSimpleFLOWMAP(file_clusters, top_npercent, distance_metric, cellnum = cellnum)
    x <- "top_npercent_edge_choice"
  }
  else {
    cellnum <- SUBSAMPLE
    graph <- buildMultiFLOWMAP(file_clusters, per = per, min = minimum,
                               max = maximum, distance_metric = distance_metric, cellnum = cellnum)
    x <- "original_edge_choice"
  }
  if (saveGRAPHML) {
    file_name <- paste(basename(FOLDER), x, sep = "_")
    convertToGraphML(graph, output_folder, file_name)
  }
  graph_xy <- forceDirectedXY(graph)
  file_name_xy <- paste(basename(FOLDER), x, "xy", sep = "_")
  convertToGraphML(graph_xy, output_folder, file_name_xy)
  if (savePDFS) {
    in_folder <- output_folder
    file_pattern <- file_name_xy
    #     visibleChannels <- c("Timepoint", "Treatment")
    convertToPDF(in_folder, file_pattern, listOfTreatments = listOfTreatments) 
    #     convertToPDF(in_folder, file_pattern, listOfTreatments = listOfTreatments,
    #                  treatInvisible = TRUE, timeInvisible = TRUE, visibleChannels = visibleChannels) 
  }
  printSummary(output_folder)
}

