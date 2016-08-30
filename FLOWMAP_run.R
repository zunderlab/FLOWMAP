rm(list=ls())

prefolder <- "/Users/mesako/Desktop/Work/Research/Code/FLOW-MAP/"
# prefolder <- "C:\\Users\\Eli\\Documents\\Lab\\Manuscripts\\2016 - CPiSCB - FLOW-MAP\\Code\\FLOW-MAP\\"
source(paste(prefolder, "FLOWMAP_main.R", sep = ""))

# MULTI_FOLDER <- "~/Desktop/Work/Research/Raw FCS Files/20150804TRAILHeLa"
# listOfTreatments <- c("DMSO", "MEKi", "p38i")
# FOLDER <- "~/Desktop/Work/Research/Raw FCS Files/20150728TRAILHeLa"
# FOLDER <- "/Users/mesako/Desktop/Work/Research/Code/FLOW-MAP_20160824/Synthetic Data"
FOLDER <- "/Users/mesako/Desktop/Work/Research/Raw FCS Files/20160721_MM1S_Timecourse3_Live/WT/J"
# FOLDER <- "/Users/mesako/Desktop/Work/Research/Raw FCS Files/20160721_MM1S_Timecourse3_Live/WT/BO"
# FOLDER <- "C:\\Users\\Eli\\Documents\\Lab\\Manuscripts\\2016 - CPiSCB - FLOW-MAP\\Code\\FLOW-MAP\\Synthetic Data\\"
FILE_FORMAT <- "*.fcs"
VAR_ANNOTATE <- list("In113Di"	= "activeBak", "La139Di"	= "cPARP",
                     "Ce140Di"	= "Bak", "Pr141Di"	= "pp38", "Nd146Di"	= "pBcl2",
                     "Nd148Di"	= "pErk", "Sm149Di"	= "APAF", "Nd150Di"	= "pRb",
                     "Sm152Di"	= "pAkt", "Eu153Di"	= "BclxL", "Sm154Di"	= "Bax",
                     "Gd155Di"	= "activeBak", "Gd156Di"	= "CyclinB", "Gd157Di"	= "Bcl2",
                     "Gd158Di"	= "pSTAT5", "Gd160Di"	= "Mcl1", "Dy161Di"	= "cMyc",
                     "Dy164Di"	= "IkBalpha", "Ho165Di"	= "Bim", "Er167Di"	= "CyclinA",
                     "Er168Di"	= "pH3", "Er170Di"	= "pBadS112", "Yb171Di"	= "pZAP70",
                     "Yb172Di"	= "Bclw", "Yb173Di"	= "cCaspase3", "Yb174Di"	= "p53",
                     "Lu175Di"	= "pS6", "Lu176Di" = "pCREB", "Pt195Di" = "Cisplatin")
VAR_REMOVE <- c("Time", "Event_length", "beadDist", "Bi209Di",
                "Pt194Di", "Ir193Di", "Ir191Di", "Tm169Di",
                "Er166Di", "Dy162Di", "Dy163Di", "Tb159Di",
                "Y89Di", "Pd102Di", "Pd104Di", "Pd105Di", "Pd106Di",
                "Pd108Di", "Pd110Di", "In115Di", "I127Di",
                "Nd142Di", "Nd143Di", "Nd144Di", "Nd145Di",
                "Pm147Di", "Eu151Di", "FileNum")
# VAR_ANNOTATE <- list("Marker1" = "Marker1", "Marker2" = "Marker2")
# VAR_REMOVE <- c()

per <- 1
minimum <- 2
maximum <- 3
distance_metric <- "manhattan" # or "euclidean"
# small
# SUBSAMPLE <- 5000
# CLUSTNUM <- 500
# large
SUBSAMPLE <- 20000
CLUSTNUM <- 1000
seedX <- 3
set.seed(seedX)

# CLUSTERING_VAR <- c("activeBak", "cPARP", "Bak", "pBcl2",
#                     "APAF","BclxL", "Bax", "activeBak", "Bcl2",
#                     "Mcl1", "Bim", "pBadS112", "Bclw", "cCaspase3")
# CLUSTERING_VAR <- c("activeBak", "cPARP", "Bak", "pp38", "pBcl2",
#                     "pErk", "APAF", "pAkt", "BclxL", "Bax",
#                     "activeBak", "Bcl2", "pSTAT5", "Mcl1",
#                     "IkBalpha", "Bim", "pBadS112", "pZAP70",
#                     "Bclw", "cCaspase3", "pS6", "pCREB")
CLUSTERING_VAR <- c("activeBak", "cPARP", "Bak", "pp38", "pBcl2",
                    "pErk", "APAF", "pRb", "pAkt", "BclxL", "Bax",
                    "activeBak", "CyclinB", "Bcl2", "pSTAT5",
                    "Mcl1", "cMyc", "IkBalpha", "Bim", "CyclinA",
                    "pH3", "pBadS112", "pZAP70", "Bclw", "cCaspase3",
                    "p53", "pS6", "pCREB")
# CLUSTERING_VAR <- c("marker1", "marker2")
setwd(FOLDER)

singleFLOWMAP(FOLDER, FILE_FORMAT, VAR_REMOVE, VAR_ANNOTATE,
              CLUSTERING_VAR, CLUSTNUM, SUBSAMPLE, distance_metric,
              minimum, maximum, per, shuffle = TRUE)


# fcs_file_names <- getFCSNames(FOLDER, FILE_FORMAT)
# output_folder <- makeOutFolder(runtype = "singleFLOWMAP")
# setwd(output_folder)
# fcs_files <- loadCleanFCS(fcs_file_names, VAR_REMOVE, VAR_ANNOTATE, subsample = SUBSAMPLE, subsampleRand = TRUE)
# for (i in 1:length(fcs_files)) {
#   df1 <- fcs_files[[i]]
#   df2 <- df1[sample(nrow(df1)), ]
#   fcs_files[[i]] <- df2
# }
# file_clusters <- clusterFCS(fcs_files = fcs_files, channel_cluster = CLUSTERING_VAR, numcluster = CLUSTNUM,
#                             distance_metric = distance_metric)
# graph <- buildFLOWMAP(file_clusters, per = per, min = minimum,
#                       max = maximum, distance_metric, cellnum = SUBSAMPLE)
# file_name <- paste(basename(FOLDER), "original_edge_choice", sep = "_")
# convertToGraphML(graph, file_name)
# graph_xy <- forceDirectedXY(graph)
# file_name_xy <- paste(basename(FOLDER), "original_edge_choice", "xy", sep = "_")
# final_file_name <- convertToGraphML(graph_xy, file_name_xy)
# convertToPDF(final_file_name, edge_color = "#FF000000")
# printSummary()


# multiFLOWMAP(listOfTreatments, MULTI_FOLDER, FILE_FORMAT, VAR_REMOVE, VAR_ANNOTATE,
#                CLUSTERING_VAR, CLUSTNUM, SUBSAMPLE, distance_metric = distance_metric,
#                per, minimum, maximum, saveGRAPHML = TRUE,
#                savePDFS = TRUE, subsampleRand = TRUE, seedX)
