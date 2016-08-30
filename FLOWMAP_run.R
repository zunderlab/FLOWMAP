rm(list=ls())

# Copyright statement comment
# Author comment
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

folder <- "/Users/mesako/Desktop/Work/Research/Code/FLOW-MAP_20160824/Synthetic Data"
# folder specifies the folder where the FCS files
# to be analyzed are saved on your local computer

file.format <- "*.fcs"
# file.format specifies what kind of files to look
# for in "folder," all functions currently only work
# with FCS files

var.annotate <- list("Marker1" = "Marker1", "Marker2" = "Marker2")
# var.annotate specifies what names each parameter
# from the FCS file should be called, for example,
# converting metal/fluorescence channels (e.g. "FITC")
# and renaming them according to marker (e.g. "CD44")

var.remove <- c()
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

subsample <- 5000
# subsample specifies how many measurements/events/cells
# to take from each FCS file, each file must contain at
# least this many events for analysis to proceed

cluster.number <- 500
# cluster.number specifies how many clusters to identify
# for the subsampled events from each separate FCS file

seed.X <- 3
set.seed(seed.X)
# seed.X specifies the seed for a given run, this should
# lead to reproducible runs of FLOW-MAP and its resulting
# figures for the same seed

clustering.var <- c("marker1", "marker2")
# clustering.var specifies which parameters to use for
# calculating distance between nodes, all markers not
# listed will not be considered and will have no influence
# on the shape of the resulting FLOW-MAP, but will still
# be seen as a parameter in the final PDFs

setwd(folder)
singleFLOWMAP(folder, file.format, var.remove, var.annotate,
              clustering.var, cluster.number, subsample, distance.metric,
              minimum, maximum, per, shuffle = TRUE)
# singleFLOWMAP function, with correctly provided folders
# and variables above, should run from start to finish,
# producing PDFs and graphml files in a new subfolder within
# the "folder" that contains the FCS files


# multiFLOWMAP(listOfTreatments, MULTI_FOLDER, FILE_FORMAT, VAR_REMOVE, VAR_ANNOTATE,
#                CLUSTERING_VAR, CLUSTNUM, SUBSAMPLE, distance_metric = distance_metric,
#                per, minimum, maximum, saveGRAPHML = TRUE,
#                savePDFS = TRUE, subsampleRand = TRUE, seedX)
