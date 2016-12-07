# GNU General Public License v3.0

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

files <- "/Users/mesako/Desktop/Work/Research/Code/FLOW-MAP/Synthetic Data/MultiFLOWMAP"
# files <- "/Users/mesako/Desktop/Work/Research/Code/FLOW-MAP/Synthetic Data/SingleFLOWMAP"
# "files" variable could be one of the following:
# a single fcs file path
# a single folder path containing 2+ fcs files
# a vector of fcs file paths
# a single folder path, containing subfolders,
# which each contain 2+ fcs files
# a list name by treatment/conditions, each element is a 
# vector of 2+ fcs file paths

save.folder <- "/Users/mesako/Desktop"
# save.folder specifies the folder where the results should be saved
# a new folder will be made in this directory named with date/time of run

file.format <- "*.fcs"
# file.format specifies what kind of files to look
# for in "folder," all functions currently only work
# with FCS files

var.annotate <- list("marker1" = "marker1", "marker2" = "marker2")
# var.annotate <- list("Pd102Di" = "barcode1", "Pd104Di" = "barcode2", "Pd105Di" = "barcode3",
#                      "Pd106Di" = "barcode4", "Pd108Di" = "barcode5", "Pd110Di" = "barcode6",
#                      "In113Di" = "active_Bax", "La139Di" = "cPARP", "Ce140Di" = "Bak",
#                      "Pr141Di" = "p-p38", "Nd146Di" = "pBcl-2", "Nd148Di" = "pErk",
#                      "Nd150Di" = "pRb", "Sm152Di" = "pAkt", "Eu153Di" = "Bcl-xL",
#                      "Sm154Di" = "Bax", "Gd155Di" = "active_Bak", "Gd156Di" = "CyclinB1",
#                      "Gd157Di" = "Bcl-2", "Gd158Di" = "pSTAT5", "Gd160Di" = "Mcl-1",
#                      "Dy161Di" = "cMyc", "Dy164Di" = "IkBalpha", "Ho165Di" = "Bim",
#                      "Er170Di" = "pBad-S112", "Er167Di" = "CyclinA", "Er168Di" = "pH3",
#                      "Yb172Di" = "Bcl-w", "Yb171Di" = "pZAP70", "Yb174Di" = "p53",
#                      "Lu175Di" = "pS6", "Lu176Di" = "pCREB", "Yb173Di" = "cCaspase3",
#                      "Ir191Di" = "DNA1", "Ir193Di" = "DNA2", "Pt195Di" = "Cisplatin", "Sm149Di" = "APAF")
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
# var.remove <- c("Time", "Event_length", "Cell_length", "beadDist", "barcode", "Nd143Di", "Nd145Di",
#                 "Normbeads1", "Normbeads2", "Eubeads", "DNA1", "DNA2", "barcode1", "Tb159Di",
#                 "barcode2", "barcode3", "Y89Di", "barcode4", "barcode5", "barcode6", "Eu151Di",
#                 "I127Di", "Gd154Di", "Dy160Di", "Tm169Di", "Dy163Di", "Pt194Di", "Nd142Di",
#                 "normbeads1", "eubeads1", "normbeads3", "FileNum", "Yb172Di", "Bi209Di", "Pm147Di",
#                 "bc_separation_dist", "Er166Di", "mahalanobis_dist", "In115Di", "Dy162Di", "Nd144Di")
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

maximum <- 5
# maximum specifies the maximum number of edges any
# given node in graph will have

distance.metric <- "manhattan" # other option is "euclidean"
# distance.metric specifies how the distance/weight
# between nodes will be calculated, in order to determine
# which edges are assigned and what is their weight

subsamples <- 1000
# subsamples <- 500
# subsamples <- c(500, 490, 510, 510, 490, 500)
# subsample specifies how many measurements/events/cells
# to take from each FCS file, each file must contain at
# least this many events for analysis to proceed

cluster.numbers <- 500
# cluster.numbers <- 250
# cluster.numbers <- c(250, 240, 255, 245, 250, 255)
# cluster.number specifies how many clusters to identify
# for the subsampled events from each separate FCS file

seed.X <- 1
set.seed(seed.X)
# seed.X specifies the seed for a given run, this should
# lead to reproducible runs of FLOW-MAP and its resulting
# figures for the same seed

clustering.var <- c("marker1", "marker2")
# clustering.var <- c("cPARP", "Bak", "p-p38", "pBcl-2", "pErk",
#                     "APAF", "pRb", "pAkt", "Bcl-xL", "Bax", "active_Bak",
#                     "CyclinB1", "Bcl-2", "pSTAT5", "Mcl-1", "cMyc", "IkBalpha",
#                     "Bim", "CyclinA", "pH3", "pBad-S112", "pZAP70", "Bcl-w",
#                     "active_Bax", "cCaspase3", "p53", "pS6", "pCREB")
# clustering.var <- c("Bid", "cCaspase3", "cPARP",
#                     "pNFkB", "IkBalpha", "Cisplatin")
# clustering.var specifies which parameters to use for
# calculating distance between nodes, all markers not
# listed will not be considered and will have no influence
# on the shape of the resulting FLOW-MAP, but will still
# be seen as a parameter in the final PDFs

# files <- list.files(folder, pattern = ".fcs", full.names = TRUE)
# one.file <- files[1]
# shuffled.files <- files[sample(x = length(files),
#                                size = length(files),
#                                replace = FALSE)]
# f <- "/Users/mesako/Desktop"

FLOWMAP.results <- FLOWMAP(files = files, file.format = file.format, var.remove = var.remove,
                           var.annotate = var.annotate, clustering.var = clustering.var,
                           cluster.numbers = cluster.numbers, subsamples = subsamples,
                           distance.metric = distance.metric, minimum = minimum, maximum = maximum,
                           per = per, save.folder = save.folder, shuffle = TRUE, name.sort = FALSE)
# FLOWMAP function, with correctly provided folders
# and variables above, should run from start to finish,
# producing PDFs and graphml files in a new subfolder within
# the "folder" that contains the FCS files

