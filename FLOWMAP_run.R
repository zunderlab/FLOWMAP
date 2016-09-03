rm(list=ls())

# Copyright: Apache License 2.0

# Authors:
# Melissa Ko (Stanford University)
# Eli Zunder (University of Virginia)

# File description comment,
# including purpose of program,
# inputs, and outputs

# source() and library() statements
# Function definitions
# Executed statements, if applicable (e.g., print, plot)

prefolder <- "/Users/mesako/Desktop/Work/Research/Code/FLOW-MAP/"
# prefolder specifies the folder where all FLOW-MAP
# code is saved on your local computer

# uses prefolder to load in all FLOW-MAP files/functions
source(paste(prefolder, "FLOWMAP_main.R", sep = ""))

folder <- "/Users/mesako/Desktop/Work/Research/Code/FLOW-MAP/Synthetic Data/SingleFLOWMAP"
# folder <- "/Users/mesako/Desktop/Work/Research/Raw FCS Files/20150728TRAILHeLa"
# folder specifies the folder where the FCS files
# to be analyzed are saved on your local computer

file.format <- "*.fcs"
# file.format specifies what kind of files to look
# for in "folder," all functions currently only work
# with FCS files

var.annotate <- list("Marker1" = "Marker1", "Marker2" = "Marker2")
# var.annotate <- list("Pd102Di" = "barcode1", "Pd104Di" = "barcode2", "Pd105Di" = "barcode3",
#                      "Pd106Di" = "barcode4", "Pd108Di" = "barcode5", "Pd110Di" = "barcode6",
#                      "In113Di" = "pH3", "I127Di" = "IdU", "La139Di" = "Normbeads1",
#                      "Pr141Di" = "pBad-S112", "Nd142Di" = "cCaspase3", "Nd143Di" = "p4EBP1",
#                      "Nd144Di" = "RSK2", "Nd145Di" = "p-p38", "Nd146Di" = "cCaspase7", 
#                      "Sm147Di" = "p-p90RSK", "Sm149Di" = "pNFkB", "Nd150Di" = "S6",
#                      "Eu151Di" = "Eubeads", "Sm152Di" = "pAkt", "Eu153Di" = "pMAPKAPK-2",
#                      "Gd156Di" = "pBcl-2", "Gd158Di" = "Bid", "Tb159Di" = "Normbeads2",
#                      "Dy162Di" = "pBad-S136", "Dy164Di" = "CyclinB1", "Ho165Di" = "pRb",
#                      "Er166Di" = "pHSP27", "Er167Di" = "pErk", "Er168Di" = "Ki67",
#                      "Tm169Di" = "IkBalpha", "Yb171Di" = "cPARP", "Yb172Di" = "pS6",
#                      "Lu175Di" = "pAMPK", "Yb176Di" = "Mcl-1", 
#                      "Ir191Di" = "DNA1", "Ir193Di" = "DNA2", "Pt195Di" = "Cisplatin")

# var.annotate specifies what names each parameter
# from the FCS file should be called, for example,
# converting metal/fluorescence channels (e.g. "FITC")
# and renaming them according to marker (e.g. "CD44")

var.remove <- c()
# var.remove <- c("Time", "Event_length", "Cell_length", "beadDist", "barcode",
#                 "Normbeads1", "Normbeads2", "Eubeads", "DNA1", "DNA2", "barcode1",
#                 "barcode2", "barcode3", "barcode4", "barcode5", "barcode6",
#                 "In115Di", "Ce140Di", "Nd148Di", "Gd154Di", "Gd155Di",
#                 "Gd157Di", "Dy160Di", "Dy161Di", "Dy163Di", "Yb170Di",
#                 "Yb173Di", "Yb174Di")
# var.remove specifies any parameters from the FCS files
# that should be entirely excluded from downstream
# analysis, specify based on post-annotation name

per <- 1
# per specifies n for top n percent of edges used
# between all nodes in graph under consideration for
# forming edges between

minimum <- 2
# minimum specifies the minimum number of edges any
# given node in graph will have

maximum <- 3
# maximum specifies the maximum number of edges any
# given node in graph will have

distance.metric <- "manhattan" # other option is "euclidean"
# distance.metric specifies how the distance/weight
# between nodes will be calculated, in order to determine
# which edges are assigned and what is their weight

subsample <- 400
# subsample specifies how many measurements/events/cells
# to take from each FCS file, each file must contain at
# least this many events for analysis to proceed

cluster.number <- 40
# cluster.number specifies how many clusters to identify
# for the subsampled events from each separate FCS file

seed.X <- 5
set.seed(seed.X)
# seed.X specifies the seed for a given run, this should
# lead to reproducible runs of FLOW-MAP and its resulting
# figures for the same seed

clustering.var <- c("marker1", "marker2")
# clustering.var <- c("Bid", "cCaspase3", "cPARP",
#                     "pNFkB", "IkBalpha", "Cisplatin")
# clustering.var specifies which parameters to use for
# calculating distance between nodes, all markers not
# listed will not be considered and will have no influence
# on the shape of the resulting FLOW-MAP, but will still
# be seen as a parameter in the final PDFs

setwd(folder)
SingleFLOWMAP(folder = folder, file.format = file.format, var.remove = var.remove,
              var.annotate = var.annotate, clustering.var = clustering.var,
              cluster.number = cluster.number, subsample = subsample,
              distance.metric = distance.metric, minimum = minimum, maximum = maximum,
              per = per, shuffle = TRUE)


# SingleFLOWMAP function, with correctly provided folders
# and variables above, should run from start to finish,
# producing PDFs and graphml files in a new subfolder within
# the "folder" that contains the FCS files

# fcs.file.names <- GetFCSNames(folder, file.format)
# output.folder <- MakeOutFolder(runtype = "singleFLOWMAP")
# setwd(output.folder)
# save.folder <- getwd()
# print(save.folder)
# fcs.files <- LoadCleanFCS(fcs.file.names, var.remove, var.annotate, subsample = subsample, subsample.rand = TRUE)
# for (i in 1:length(fcs.files)) {
#   df1 <- fcs.files[[i]]
#   df2 <- df1[sample(nrow(df1)), ]
#   fcs.files[[i]] <- df2
# }
# file.clusters <- ClusterFCS(fcs.files, clustering.var, numcluster = cluster.number,
#                             distance.metric = distance.metric)
# graph <- BuildFLOWMAP(file.clusters, per = per, min = minimum,
#                       max = maximum, distance.metric, cellnum = subsample)
# file.name <- paste(basename(folder), "original_edge_choice", sep = "_")
# ConvertToGraphML(graph, file.name)
# graph.xy <- ForceDirectedXY(graph)
# file.name.xy <- paste(basename(folder), "original_edge_choice", "xy", sep = "_")
# final.file.name <- ConvertToGraphML(graph.xy, file.name.xy)
# ConvertToPDF(final.file.name, edge.color = "#FF000000")
# print(getwd())
# setwd(save.folder)
# printSummary()



# multiFLOWMAP(listOfTreatments, MULTI_FOLDER, FILE_FORMAT, VAR_REMOVE, VAR_ANNOTATE,
#                CLUSTERING_VAR, CLUSTNUM, SUBSAMPLE, distance_metric = distance_metric,
#                per, minimum, maximum, saveGRAPHML = TRUE,
#                savePDFS = TRUE, subsampleRand = TRUE, seedX)

# fcs.file.names <- GetMultiFCSNames(folder, file.format)
# output.folder <- MakeOutFolder(runtype = "multiFLOWMAP")
# setwd(output.folder)
# save.folder <- getwd()
# print(save.folder)
# fcs.files <- LoadMultiCleanFCS(fcs.file.names, var.remove, var.annotate,
#                                subsample = subsample, subsample.rand)
# file.clusters <- MultiClusterFCS(fcs.files, channel.cluster = clustering.var, numcluster = cluster.number)
# graph <- BuildMultiFLOWMAP(file.clusters, per = per, min = minimum,
#                            max = maximum, distance.metric = distance.metric, cellnum = subsample)
# file.name <- paste(basename(folder), "original_edge_choice", sep = "_")
# ConvertToGraphML(graph, file.name)
# graph.xy <- forceDirectedXY(graph)
# file.name.xy <- paste(basename(folder), "original_edge_choice", "xy", sep = "_")
# final.file.name <- ConvertToGraphML(graph.xy, file.name.xy)
# ConvertToPDF(final.file.name, edge_color = "#FF000000")
# # visibleChannels <- c("Timepoint", "Treatment")
# ConvertToPDF(in_folder, file_pattern, listOfTreatments = listOfTreatments) 
# # convertToPDF(in_folder, file_pattern, listOfTreatments = listOfTreatments,
# #             treatInvisible = TRUE, timeInvisible = TRUE, visibleChannels = visibleChannels) 
# print(getwd())
# setwd(save.folder)
# printSummary()
