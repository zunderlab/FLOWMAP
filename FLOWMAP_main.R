rm(list=ls())

library(igraph)
library(Rclusterpp)
library(flowCore)
library(robustbase)
library(cluster)
library(Biobase)
library(RForceAtlas2)
library(huge)
library(proxy)

setwd("~/Desktop/Work/Research/FLOW-MAP Code/FLOWMAP")
source("FLOWMAP_loadFromFCS.R")
source("FLOWMAP_clusterCells.R")
source("FLOWMAP_buildGraph.R")
source("FLOWMAP_buildMultiGraph.R")
source("FLOWMAP_forceDirected.R")
source("FLOWMAP_cosineDist.R")
source("FLOWMAP_saveGraphs.R")
source("FLOWMAP_exportEvents.R")


# MULTI_FOLDER <- "~/Desktop/Work/Research/Raw FCS Files/20150804TRAILHeLa"
# listOfTreatments <- c("DMSO", "MEKi", "p38i")

FOLDER <- "~/Desktop/Work/Research/Raw FCS Files/20150728TRAILHeLa"

FILE_FORMAT <- "*.fcs"
VAR_ANNOTATE <- list("Pd102Di" = "barcode1", "Pd104Di" = "barcode2", "Pd105Di" = "barcode3",
                     "Pd106Di" = "barcode4", "Pd108Di" = "barcode5", "Pd110Di" = "barcode6",
                     "In113Di" = "pH3", "I127Di" = "IdU", "La139Di" = "Normbeads1",
                     "Pr141Di" = "pBad-S112", "Nd142Di" = "cCaspase3", "Nd143Di" = "p4EBP1",
                     "Nd144Di" = "RSK2", "Nd145Di" = "p-p38", "Nd146Di" = "cCaspase7", 
                     "Sm147Di" = "p-p90RSK", "Sm149Di" = "pNFkB", "Nd150Di" = "S6",
                     "Eu151Di" = "Eubeads", "Sm152Di" = "pAkt", "Eu153Di" = "pMAPKAPK-2",
                     "Gd156Di" = "pBcl-2", "Gd158Di" = "Bid", "Tb159Di" = "Normbeads2",
                     "Dy162Di" = "pBad-S136", "Dy164Di" = "CyclinB1", "Ho165Di" = "pRb",
                     "Er166Di" = "pHSP27", "Er167Di" = "pErk", "Er168Di" = "Ki67",
                     "Tm169Di" = "IkBalpha", "Yb171Di" = "cPARP", "Yb172Di" = "pS6",
                     "Lu175Di" = "pAMPK", "Yb176Di" = "Mcl-1", 
                     "Ir191Di" = "DNA1", "Ir193Di" = "DNA2", "Pt195Di" = "Cisplatin")
VAR_REMOVE <- c("Time", "Event_length", "Cell_length", "beadDist", "barcode",
                "Normbeads1", "Normbeads2", "Eubeads", "DNA1", "DNA2", "barcode1",
                "barcode2", "barcode3", "barcode4", "barcode5", "barcode6",
                "In115Di", "Ce140Di", "Nd148Di", "Gd154Di", "Gd155Di",
                "Gd157Di", "Dy160Di", "Dy161Di", "Dy163Di", "Yb170Di",
                "Yb173Di", "Yb174Di")

# use all variables
CLUSTERING_VAR <- c()
for (n in names(VAR_ANNOTATE)) {
  CLUSTERING_VAR <- c(CLUSTERING_VAR, VAR_ANNOTATE[[n]])
}
CLUSTERING_VAR <- setdiff(CLUSTERING_VAR, VAR_REMOVE)

# VAR_UPSAMPLE <- "cCaspase3"


singleFLOWMAP <- function(FOLDER, FILE_FORMAT, VAR_REMOVE, VAR_ANNOTATE,
                          CLUSTERING_VAR, CLUSTNUM, saveGRAPHML = TRUE,
                          savePDFS = TRUE, subsampleRand = TRUE,
                          findStages = TRUE, exportStages = TRUE, ...) {
  fcs_file_names <- getFCSNames(FOLDER, FILE_FORMAT)
  output_folder <- "~/Desktop"
  name <- gsub(" ", "_", Sys.time(), fixed = TRUE)
  name <- gsub(":", ".", name, fixed = TRUE)
  output_folder <- paste(output_folder, "/", name, "_singleFLOWMAP_run", sep = "")
  dir.create(output_folder)
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
  graph_xy <- forceDirectedXY(graph)
  file_name_xy <- paste(basename(FOLDER), x, "xy", sep = "_")
  convertToGraphML(graph_xy, output_folder, file_name_xy)
  if (savePDFS) {
    in_folder <- output_folder
    file_pattern <- file_name_xy
#     visibleChannels <- c("Timepoint", "Treatment")
    convertToPDF(in_folder, file_pattern) 
#     convertToPDF(in_folder, file_pattern, listOfTreatments = listOfTreatments,
#                  treatInvisible = TRUE, timeInvisible = TRUE, visibleChannels = visibleChannels) 
    printSummary(output_folder)
  }
}


multiFLOWMAP <- function(listOfTreatments, FOLDER, FILE_FORMAT, VAR_REMOVE,
                         VAR_ANNOTATE, CLUSTERING_VAR, CLUSTNUM, saveGRAPHML = TRUE,
                         savePDFS = TRUE, subsampleRand = TRUE, ...) {
  fcs_file_names <- getMultiFCSNames(listOfTreatments, FOLDER, FILE_FORMAT)
  output_folder <- "~/Desktop"
  name <- gsub(" ", "_", Sys.time(), fixed = TRUE)
  name <- gsub(":", ".", name, fixed = TRUE)
  output_folder <- paste(output_folder, "/", name, "_multiFLOWMAP_run", sep = "")
  dir.create(output_folder)
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
    printSummary(output_folder)
  }
}


per <- 1
minimum <- 2
maximum <- 3
distance_metric <- "manhattan"
# distance_metric <- "cosine" or "euclidean"
SUBSAMPLE <- 1000
CLUSTNUM <- 100
numStages <- 10
seedX <- 10
set.seed(seedX)

singleFLOWMAP(FOLDER, FILE_FORMAT, VAR_REMOVE, VAR_ANNOTATE,
              CLUSTERING_VAR, CLUSTNUM, SUBSAMPLE, distance_metric,
              per, minimum, maximum, saveGRAPHML = TRUE,
              savePDFS = TRUE, subsampleRand = TRUE,
              findStages = TRUE, exportStages = TRUE, seedX,
              channel_cluster = CLUSTERING_VAR)

# multiFLOWMAP(listOfTreatments, MULTI_FOLDER, FILE_FORMAT, VAR_REMOVE, VAR_ANNOTATE,
#                CLUSTERING_VAR, CLUSTNUM, SUBSAMPLE, distance_metric = distance_metric,
#                per, minimum, maximum, saveGRAPHML = TRUE,
#                savePDFS = TRUE, subsampleRand = TRUE, seedX)

# in_folder <- "~/Desktop/2015-09-11_00.59.06_multiFLOWMAP_run"
# visibleChannels <- c("Timepoint", "Treatment")
# convertToPDF(in_folder, file_pattern = "DEFAULT", out_folder = in_folder,
#              scale = NULL, normalize = "none", node_size_scale = 2,
#              min_node_size = 12, max_node_size = 24, pdf_width = 100, pdf_height = 100,
#              text_color = "black", PALETTE = "jet", listOfTreatments = listOfTreatments,
#              treatInvisible = TRUE, timeInvisible = TRUE, visibleChannels) 



