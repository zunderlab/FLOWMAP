library(igraph)
library(proxy)


remodelFLOWMAPClusterList <- function(listOfFLOWMAPclusters) {
  # take FLOWMAP of treats with timeseries and make into one timeseries
  # combine FLOWMAP treatments
  treatments <- names(listOfFLOWMAPclusters)
  timepoints <- length(listOfFLOWMAPclusters[[treatments[1]]]$cluster_medians)
  fullclusters <- data.frame()
  table_breaks <- c()
  table_lengths <- c()
  cluster_medians <- list()
  cluster_counts <- list()
  cellassgn <- list()
#   print(timepoints)
  for (t in 1:timepoints) {
#     print(t)
    temp_medians <- data.frame()
    temp_cellassgn <- data.frame()
    temp_counts <- data.frame()
    for (treat in treatments) {
#       print(treat)
      temp_medians <- rbind(temp_medians, listOfFLOWMAPclusters[[treat]]$cluster_medians[[t]])
      temp_cellassgn <- rbind(temp_cellassgn, listOfFLOWMAPclusters[[treat]]$cellassgn[[t]])
      temp_counts <- rbind(temp_counts, listOfFLOWMAPclusters[[treat]]$cluster_counts[[t]])
#       print(dim(temp_counts))
#       print(head(temp_counts))
#       print(dim(temp_cellassgn))
#       print(head(temp_cellassgn))
    }
    cluster_medians[[t]] <- temp_medians
    cluster_counts[[t]] <- temp_counts
    cellassgn[[t]] <- temp_cellassgn
    table_lengths <- c(table_lengths, dim(temp_medians)[1])
    table_breaks <- c(table_breaks, sum(table_lengths))
    fullclusters <- rbind(fullclusters, temp_medians)
  }
  remodeledFLOWMAPclusters <- FLOWMAPcluster(fullclusters, table_breaks, table_lengths,
                                             cluster_medians, cluster_counts, cellassgn) 
  return(remodeledFLOWMAPclusters)
}


initializeMultiGraph <- function(listOfFLOWMAPclusters) {
  # This section initializes the graph
  total_nodes <- 0
  for (treat in names(listOfFLOWMAPclusters)) {
    FLOWMAPclusters <- listOfFLOWMAPclusters[[treat]]
    total_nodes <- total_nodes + length(FLOWMAPclusters$fullclusters[, 1])
    # create empty graph with the right number of nodes - will fill in the edges later
  }
  initial_graph <- graph.empty(n = total_nodes, directed = FALSE)
  return(initial_graph)
}


buildFirstMultiFLOWMAP <- function(listOfFLOWMAPclusters, per, min, max, distance_metric, ...) {
  n <- 1
  output_graph <- initializeMultiGraph(listOfFLOWMAPclusters)
  # This section creates a flowmap for the first time point
  cat("Building first flowmap\n")
  clusters <- c()
  for (treat in names(listOfFLOWMAPclusters)) {
    table_length <- listOfFLOWMAPclusters[[treat]]$table_lengths[1]
    current_clusters <- listOfFLOWMAPclusters[[treat]]
    clusters <- rbind(clusters, current_clusters$fullclusters[1:table_length, ])
  }
  # get distance matrix from clusters
  clusters <- subset(clusters, select = colnames(clusters)[colnames(clusters) %in% CLUSTERING_VAR])
  if (distance_metric == "cosine") {
    cluster_distances <- cosine_similarity_matrix(clusters)
  }
  else {
    cluster_distances <- dist(clusters, method = distance_metric, diag = TRUE, upper = TRUE)
  }
  cluster_distances_matrix <- as.matrix(cluster_distances)
  # set i-i distances to Inf instead of 0 so they aren't the closest neighbors
  for (i in 1:ncol(cluster_distances_matrix)) {
    cluster_distances_matrix[i, i] <- Inf
  }
  numcluster <- nrow(clusters)
  normalized_densities <- findNormalized(cluster_distances_matrix,
                                         per, min, max, numcluster)
  # build new edgelist with N edges for each cluster based on normalized density
  cat("Building first edgelist\n")
  output_graph <- drawNormalizedEdges(output_graph, cluster_distances_matrix,
                                      normalized_densities, n = n)
  # now add all MST edges that are not yet included in the graph, and annotate all as "MST"
  adjacency_graph <- graph.adjacency(cluster_distances_matrix,
                                     mode = "undirected", weighted = "weight", diag = FALSE)
  mst_graph <- minimum.spanning.tree(adjacency_graph, algorithm = "prim")
  # for each edge of the mst, if it exists in the graph then label MST, if it doesn't exist then add it
  mst_graph_edgelist <- cbind(get.edgelist(mst_graph), E(mst_graph)$weight)
  class(mst_graph_edgelist) <- "numeric"
  # for each edge of the mst, if it exists in the graph then label MST, if it doesn't exist then add it
  for (i in 1:nrow(mst_graph_edgelist)) {
    if (are.connected(output_graph,
                      mst_graph_edgelist[i, 1],
                      mst_graph_edgelist[i, 2])) {
      E(output_graph,P = c(mst_graph_edgelist[i, 1],
                           mst_graph_edgelist[i, 2]))$label <- "MST"
    } else {
      output_graph <- add.edges(output_graph,
                                as.numeric(mst_graph_edgelist[i, 1:2]),
                                weight = mst_graph_edgelist[i, 3],
                                label = "MST", sequence_assignment = 1)
    }
  }
  return(output_graph)
}


buildMultiFLOWMAP <- function(listofFLOWMAPclusters, treatments,
                              per, min, max, distance_metric, cellnum, ...) {
  output_graph <- buildFirstMultiFLOWMAP(listofFLOWMAPclusters, per, min, max, distance_metric = distance_metric)
  remodelFLOWMAPclusters <- remodelFLOWMAPClusterList(listofFLOWMAPclusters)
  # put each treatments clusters together into one timepoint
  table_breaks <- c(0, remodelFLOWMAPclusters$table_breaks)
  # This section builds the flowmap one timepoint at a time
  for (n in 1:(length(remodelFLOWMAPclusters$cluster_medians) - 1)) {
    # offset value is used to correctly index the edges at each sequential step
    offset <- table_breaks[n]
    # go through sequential cluster sets, add edges for each a-a and a-a+1 set, and also label sequential MST
    cat("Build FlowMap from", n, "to", n + 1, "\n")
    # get clusters for time a and a+1
    clusters <- rbind(remodelFLOWMAPclusters$cluster_medians[[n]], remodelFLOWMAPclusters$cluster_medians[[n + 1]])
    clusters <- subset(clusters, select = colnames(clusters)[colnames(clusters) %in% CLUSTERING_VAR])
    numcluster <- nrow(clusters)
    # make adjacency matrix from clusters
    if (distance_metric == "cosine") {
      cluster_distances <- cosine_similarity_matrix(clusters)
    }
    else {
      cluster_distances <- dist(clusters, method = distance_metric, diag = TRUE, upper = TRUE)
    }
    cluster_distances_matrix <- as.matrix(cluster_distances)
    # set i-i distances to Inf instead of 0 so they aren't the closest neighbors
    for (i in 1:ncol(cluster_distances_matrix)) {
      cluster_distances_matrix[i, i] <- Inf
    }
    # This section adds the lowest distance n_n+1 and n+1_n+1 edges to the output graph
    n_n1_table_lengths <- remodelFLOWMAPclusters$table_lengths[n:(n + 1)]
    n_n1_table_breaks <- remodelFLOWMAPclusters$table_breaks[n:(n + 1)]
    normalized_densities <- findNormalized(cluster_distances_matrix, per, min,
                                           max, numcluster, table_lengths = n_n1_table_lengths,
                                           table_breaks = n_n1_table_breaks,
                                           offset = offset)
    # build new edgelist with N edges for each cluster based on normalized density
    cat("Building edgelist\n")
    output_graph <- drawNormalizedEdges(output_graph, cluster_distances_matrix,
                                        normalized_densities, offset = offset, n = n)
    # This section adds the "MST" for n_n+1 and n+1_n+1 nodes
    output_graph <- checkMSTedges(output_graph, cluster_distances_matrix, 
                                  n_n1_table_lengths, offset = offset, n = n)
  }
  # convert graph distances to weights (low distance = high weight and vice versa)
  distances <- E(output_graph)$weight
  #weights <- -distances+max(distances)+min(distances)
  weights <- 1 / distances
  E(output_graph)$weight <- weights
  output_graph <- annotateMultiGraph(output_graph, listofFLOWMAPclusters, cellnum)
  return(output_graph)
}


annotateMultiGraph <- function(output_graph, listOfFLOWMAPclusters, cellnum, ...) {
  # This section annotates the graph
  anno_cat <- c()
  anno <- list()
  # iterate through all times and annotate
  x <- names(listOfFLOWMAPclusters)[1]
  
  for (f in 1:length(listOfFLOWMAPclusters[[x]]$cluster_medians)) {
    for (treat in names(listOfFLOWMAPclusters)) {
      cat("Annotating graph for file", f, "\n")
      # get medians for all parameters and counts for all clusters
      counts <- listOfFLOWMAPclusters[[treat]]$cluster_counts[[f]]$Counts
      anno$count <- counts
      anno$percenttotal <- data.frame(percenttotal = c(counts / cellnum))
      anno$medians  <- listOfFLOWMAPclusters[[treat]]$cluster_medians[[f]]
      # add time information column
      time_matrix <- matrix(f, nrow = length(anno$count))
      colnames(time_matrix) <- c("Timepoint")
      # add treat information column
      treat_matrix <- matrix(treat, nrow = length(anno$count))
      colnames(treat_matrix) <- c("Treatment")
      anno$medians <- cbind(anno$medians, time_matrix, treat_matrix)
      # add median and percent values
      for (a in c("medians", "percenttotal")) {
        # anno_cat is all the annos concatenated, will be
        # used to make "anno.Rsave" file
        anno_cat[[a]] <- rbind(anno_cat[[a]], anno[[a]])
      }
    }
  }
  
  #   for (treat in names(listOfFLOWMAPclusters)) {
  #     for (f in 1:length(listOfFLOWMAPclusters[[treat]]$cluster_medians)) {
  #       cat("Annotating graph for file", f, "\n")
  #       # get medians for all parameters and counts for all clusters
  #       counts <- listOfFLOWMAPclusters[[treat]]$cluster_counts[[f]]$Counts
  #       anno$count <- counts
  #       anno$percenttotal <- data.frame(percenttotal = c(counts / cellnum))
  #       anno$medians  <- listOfFLOWMAPclusters[[treat]]$cluster_medians[[f]]
  #       # add time information column
  #       time_matrix <- matrix(f, nrow = length(anno$count))
  #       colnames(time_matrix) <- c("Timepoint")
  #       # add treat information column
  #       treat_matrix <- matrix(treat, nrow = length(anno$count))
  #       colnames(treat_matrix) <- c("Treatment")
  #       anno$medians <- cbind(anno$medians, time_matrix)
  #       # add median and percent values
  #       for (a in c("medians", "percenttotal")) {
  #         # anno_cat is all the annos concatenated, will be
  #         # used to make "anno.Rsave" file
  #         anno_cat[[a]] <- rbind(anno_cat[[a]], anno[[a]])
  #       }
  #     }
  #   }
  
  # combine anno_cat matrices
  output_anno <- cbind(anno_cat[[1]], anno_cat[[2]])
  for (c in colnames(output_anno)) {
    output_graph <- set.vertex.attribute(output_graph, c,
                                         index = as.numeric(rownames(output_anno)),
                                         value = output_anno[, c])
  }
  # add name attribute
  V(output_graph)$name <- 1:length(V(output_graph))
  return(output_graph)
}






# 
# fcs_files <- list()
# 
# for (treat in INDIV_DIR_NAMES) {''
#   dir_name <- paste(FOLDER,"/",treat,sep="")
#   # get FCS files
#   fcs_files[[treat]] = list.files(path=dir_name,pattern=FILE_FORMAT,recursive=FALSE,full.names=TRUE)
#   # sort to organize by hour
#   fcs_files[[treat]] <- sort(fcs_files[[treat]])
# }
# 
# clusters_full <- data.frame() # list of lists????
# table_breaks   <- c() # this will be used later in the annotation section
# table_lengths <- c()
# cluster_tables <- list()
# cluster_medians <- list()
# cluster_counts <- list()
# 
# singlecellclust <- as.data.frame(matrix(0, ncol = 3, nrow = 0))
# colnames(singlecellclust) <- c("Treatment","ClusterID","SingleCellIndices")
# clustcount <- 0
# 
# cellassgn_matrix <- data.frame(Treat = numeric(0),Timepoint = numeric(0),Cluster = numeric(0),Members = numeric(0))
# 
# for (treat in INDIV_DIR_NAMES) {
#   for (t in 1:length(cellassgn[[treat]])) {
#     for (i in 1:(max(cellassgn[[treat]][[t]]))) {
#       cellmembers <- which(cellassgn[[treat]][[t]] == i)
#       cellmembers <- paste(cellmembers,collapse=", ")
#       tmp <- cbind(treat,t,i,cellmembers)
#       cellassgn_matrix <- rbind(cellassgn_matrix,tmp)
#     }  
#   }
# }
# 
# write.csv(cellassgn_matrix, file = paste(FOLDER,"/practice_alltreat_cell_cluster_assignments.csv",sep=""))
# 
# # create empty graph with the right number of nodes - will fill in the edges later
# total_nodes <- length(clusters_full[,1])
# output_graph <- graph.empty(n=total_nodes,directed=FALSE)
# 
# sequence_assignments <- c()
# keep = 0
# numfiles <- length(cluster_tables[[1]])
# 
# for (treat in INDIV_DIR_NAMES) {
#   for (i in 1:numfiles) {
#     sequence_assignments <- c(sequence_assignments,rep(i,table_lengths[i + keep * numfiles]))
#   }
#   keep = keep + 1
# }
# 
# E(output_graph)$sequence_assignment <- sequence_assignments
# 
# 
# cat("building first flowmap\n")
# # get distance matrix from clusters
# 
# # for pairwise comparison
# #       clusters <- rbind(clusters_full[(firstind[1]:firstind[2]),],
# #                         clusters_full[(firstind[3]:firstind[4]),])
# 
# initialclusters <- clusters_full[which(sequence_assignments==1),]
# 
# # remove medians for markers not involved with clustering/distance
# clustfordist <- subset(initialclusters, select = colnames(initialclusters)[colnames(initialclusters) %in% CLUSTERING_VAR])
# 
# cluster_distances <- dist(clustfordist, method=DISTANCE_METHOD)
# cluster_distances_matrix <- as.matrix(cluster_distances)
# # set i-i distances to Inf instead of 0 so they aren't the closest neighbors
# for (i in 1:ncol(cluster_distances_matrix)) {
#   cluster_distances_matrix[i,i] <- Inf
# }
# # make fully connected graph from adjacency list
# # note: "weight" here is really distance, calling it weight for the mst function later needs this
# edgelist_with_distances <- cbind(as.vector(row(cluster_distances_matrix)),
#                                  as.vector(col(cluster_distances_matrix)),as.vector(cluster_distances_matrix))
# # take only the upper triangle of the matrix
# edgelist_with_distances <- edgelist_with_distances[upper.tri(matrix(data=1:length(cluster_distances_matrix),
#                                                                     nrow=nrow(cluster_distances_matrix),ncol=ncol(cluster_distances_matrix))),]
# # convert strings to numeric
# class(edgelist_with_distances) <- "numeric"
# # sort edges by distance (column 3)
# edgelist_with_distances <- edgelist_with_distances[order(edgelist_with_distances[,3]),]
# # take only the top X-% edges by distance based on PERCENT_TOTAL
# trim_edgelist_with_distances <- edgelist_with_distances[1:floor(length(edgelist_with_distances[,1])*per/100),]
# # calculate "density" for each cluster and normalize
# densities_no_zeros <- table(trim_edgelist_with_distances[,1:2])
# # add in zeros for clusters with no edges (table function leaves these out)
# 
# indtoclust <- list()
# 
# for (i in 1:length(rownames(initialclusters))) {
#   indtoclust[[as.character(i)]] <- rownames(initialclusters)[i]
# }
# 
# for (elem in names(densities_no_zeros)) {
#   elemind <- match(elem,names(densities_no_zeros))
#   names(densities_no_zeros)[elemind] <- indtoclust[elem][[1]]
# }
# 
# densities <- rep(0,nrow(initialclusters))
# names(densities) <- rownames(initialclusters)
# densities[names(densities_no_zeros)] <- densities_no_zeros
# normalized_densities <- round(densities/max(densities)*(max - min) + min)
# 
# # build new edgelist with N edges for each cluster based on normalized density
# cat("building first edgelist\n")
# final_edgelist_with_distances <- c()
# for (n in 1:length(normalized_densities)) {
#   clustID <- names(normalized_densities)[n]        
#   tmp_final_edgelist_with_distances <- cbind(as.numeric(clustID),order(cluster_distances_matrix[,clustID]),
#                                              unname(sort(cluster_distances_matrix[,clustID])))[1:normalized_densities[clustID],]
#   for (i in 1:length(tmp_final_edgelist_with_distances[,2])) {
#     tmp_final_edgelist_with_distances[,2][i] <- as.numeric(indtoclust[tmp_final_edgelist_with_distances[,2][i]][[1]])
#   }
#   
#   final_edgelist_with_distances <- rbind(final_edgelist_with_distances,tmp_final_edgelist_with_distances)
# }
# 
# # remove duplicate edges from edgelist
# final_edgelist_with_distances <- unique(cbind(t(apply(final_edgelist_with_distances[,1:2],1,sort)),final_edgelist_with_distances[,3]))
# 
# 
# # add edges to graph
# # note: weight here is really distance, will be converted to weight after the graph is completed
# output_graph <- add.edges(output_graph,edges=as.vector(t(final_edgelist_with_distances[,1:2])),
#                           weight=final_edgelist_with_distances[,3],label="EXTRA",sequence_assignment=1)
# 
# # now add all MST edges that are not yet included in the graph, and annotate all as "MST"
# adjacency_graph <- graph.adjacency(cluster_distances_matrix, mode="undirected", weighted="weight", diag=FALSE)
# mst_graph <- minimum.spanning.tree(adjacency_graph, algorithm="prim")
# # for each edge of the mst, if it exists in the graph then label MST, if it doesn't exist then add it
# mst_graph_edgelist <- cbind(get.edgelist(mst_graph),E(mst_graph)$weight)
# class(mst_graph_edgelist) <- "numeric"
# # for each edge of the mst, if it exists in the graph then label MST, if it doesn't exist then add it
# for (i in 1:nrow(mst_graph_edgelist)) {
#   if (are.connected(output_graph,mst_graph_edgelist[i,1],mst_graph_edgelist[i,2])) {
#     E(output_graph,P=c(mst_graph_edgelist[i,1],mst_graph_edgelist[i,2]))$label <- "MST"
#   } else {
#     output_graph <- add.edges(output_graph,as.numeric(mst_graph_edgelist[i,1:2]),weight=mst_graph_edgelist[i,3],label="MST",sequence_assignment=1)
#   }
# }
# 
# rm(indtoclust)
# 
# ################################################################################
# # This section builds the flowmap one timepoint at a time
# ################################################################################
# 
# # offset value is used to correctly index the edges at each sequential step
# #offset=0
# 
# # go through sequential cluster sets, add edges for each a-a and a-a+1 set, and also label sequential MST
# for (n in 1:(numfiles-1)) {
#   cat(paste("Build FlowMap from",n,"to",n+1,"\n"))
#   # get clusters for time a and a+1
#   
#   clusters <- rbind(clusters_full[which(sequence_assignments==n),],clusters_full[which(sequence_assignments==(n+1)),])
#   
#   clustfordist <- subset(clusters, select = colnames(clusters)[colnames(clusters) %in% CLUSTERING_VAR])
#   
#   #clusters <- rbind(read.table(cluster_table_files[n], header=TRUE, check.names=FALSE),read.table(cluster_table_files[n+1], header=TRUE, check.names=FALSE))
#   # make adjacency matrix from clusters
#   cluster_distances <- dist(clustfordist, method=DISTANCE_METHOD)
#   cluster_distances_matrix <- as.matrix(cluster_distances)
#   # set i-i distances to Inf instead of 0 so they aren't the closest neighbors
#   for (i in 1:ncol(cluster_distances_matrix)) {
#     cluster_distances_matrix[i,i] <- Inf
#   }
#   
#   ###################################################################################
#   # This section adds the lowest distance n_n+1 and n+1_n+1 edges to the output graph
#   ###################################################################################
#   
#   # take the bottom half of the matrix, giving (n to n+1) and (n+1 to n+1) distances
#   # as separate matrices
#   
#   lengthofn <- length(rownames(clusters_full[which(sequence_assignments==n),]))
#   
#   n_n1 <- cluster_distances_matrix[(lengthofn + 1):dim(cluster_distances_matrix)[1],
#                                    1:lengthofn]
#   n1_n1 <- cluster_distances_matrix[(lengthofn + 1):dim(cluster_distances_matrix)[1],
#                                     (lengthofn + 1):dim(cluster_distances_matrix)[1]]
#   
#   #         n_n1 <- cluster_distances_matrix[(table_lengths[n]+1):(table_lengths[n]+table_lengths[n+1]),
#   #                                          1:table_lengths[n]]
#   #         n1_n1 <- cluster_distances_matrix[(table_lengths[n]+1):(table_lengths[n]+table_lengths[n+1]),
#   #                                           (table_lengths[n]+1):(table_lengths[n]+table_lengths[n+1])]
#   # convert the distance matrices to lists of coordinates with distances as the 3rd row
#   edgelist_with_distances_n_n1 <- cbind(rep(colnames(n_n1),each=length(rownames(n_n1))),
#                                         rep(rownames(n_n1),times=length(colnames(n_n1))),as.vector(n_n1))
#   edgelist_with_distances_n1_n1 <- cbind(rep(colnames(n1_n1),each=length(rownames(n1_n1))),
#                                          rep(rownames(n1_n1),times=length(colnames(n1_n1))),as.vector(n1_n1))
#   
#   #         edgelist_with_distances_n_n1 <- cbind(as.vector(row(n_n1)+table_breaks[n]),
#   #                                               as.vector(col(n_n1)+table_breaks[n]-table_lengths[n]),as.vector(n_n1))
#   #         edgelist_with_distances_n1_n1 <- cbind(as.vector(row(n1_n1)+table_breaks[n]),
#   #                                                as.vector(col(n1_n1)+table_breaks[n]),as.vector(n1_n1))
#   # only take upper triangle from n+1 to n+1 matrix
#   
#   edgelist_with_distances_n1_n1 <- edgelist_with_distances_n1_n1[upper.tri(matrix(data=1:length(n1_n1),
#                                                                                   nrow=nrow(n1_n1),ncol=ncol(n1_n1))),]
#   
#   # edgelist_with_distances_n1_n1 <- edgelist_with_distances_n1_n1[edgelist_with_distances_n1_n1[,1] != edgelist_with_distances_n1_n1[,2],,drop=FALSE]
#   
#   #combine n_n1 and n1_n1 distances into a single "edgelist" with distances
#   edgelist_with_distances <- rbind(edgelist_with_distances_n_n1,edgelist_with_distances_n1_n1)
#   # convert strings to numeric
#   class(edgelist_with_distances) <- "numeric"
#   # sort edges by distance (column 3)
#   edgelist_with_distances <- edgelist_with_distances[order(edgelist_with_distances[,3]),]
#   # take only the top X-% edges by distance based on PERCENT_TOTAL
#   trim_edgelist_with_distances <- edgelist_with_distances[1:floor(length(edgelist_with_distances[,1])*per/100),]
#   # add edges to graph
#   # note: weight here is really distance, will be converted to weight after the graph is completed
#   
#   # calculate "density" for each cluster and normalize
#   densities_no_zeros <- table(trim_edgelist_with_distances[,1:2])
#   # add in zeros for clusters with no edges (table function leaves these out)
#   
#   indtoclust <- list()
#   
#   for (i in 1:length(rownames(clusters))) {
#     indtoclust[[as.character(i)]] <- rownames(clusters)[i]
#   }
#   
#   #         for (elem in names(densities_no_zeros)) {
#   #           elemind <- match(elem,names(densities_no_zeros))
#   #           names(densities_no_zeros)[elemind] <- indtoclust[elem][[1]]
#   #         }
#   #         
#   densities <- rep(0,nrow(clusters))
#   names(densities) <- rownames(clusters)        
#   densities[names(densities_no_zeros)] <- densities_no_zeros
#   normalized_densities <- round(densities/max(densities)*(max - min) + min)        
#   
#   cat("building edgelist\n")
#   final_edgelist_with_distances <- c()
#   for (n in 1:length(normalized_densities)) {
#     clustID <- names(normalized_densities)[n]        
#     tmp_final_edgelist_with_distances <- cbind(as.numeric(clustID),order(cluster_distances_matrix[,clustID]),
#                                                unname(sort(cluster_distances_matrix[,clustID])))[1:normalized_densities[clustID],]
#     for (i in 1:length(tmp_final_edgelist_with_distances[,2])) {
#       tmp_final_edgelist_with_distances[,2][i] <- as.numeric(indtoclust[tmp_final_edgelist_with_distances[,2][i]][[1]])
#     }
#     
#     final_edgelist_with_distances <- rbind(final_edgelist_with_distances,tmp_final_edgelist_with_distances)
#   }
#   
#   # remove duplicate edges from edgelist
#   final_edgelist_with_distances <- unique(cbind(t(apply(final_edgelist_with_distances[,1:2],1,sort)),final_edgelist_with_distances[,3]))
#   
#   output_graph <- add.edges(output_graph,edges=as.vector(t(final_edgelist_with_distances[,1:2])),
#                             weight=final_edgelist_with_distances[,3],label="EXTRA",sequence_assignment=n+1)
#   
#   ###################################################################################
#   # This section adds the "MST" for n_n+1 and n+1_n+1 nodes, a little complicated. . .
#   ###################################################################################
#   
#   # create dummy graph contains n and n+1 vertices, but n_n section is fully connected so
#   # it will act as a single connected component while n_n+1 and n+1_n+1 sections are not
#   # connected, so each n+1 vertex will be its own disconnected component
#   adjacency_nn <- cluster_distances_matrix[1:lengthofn,1:lengthofn]
#   el_nn_with_dist <- cbind(as.vector(row(adjacency_nn)),as.vector(col(adjacency_nn)),as.vector(adjacency_nn))
#   dummy_graph <- graph.empty(n=nrow(cluster_distances_matrix), directed=FALSE)
#   dummy_graph <- add.edges(dummy_graph,edges=as.vector(t(el_nn_with_dist[,1:2])))
#   
#   #                 first <- rep(as.numeric(rownames(adjacency_nn)),times=length(rownames(adjacency_nn)))
#   #                 second <- rep(as.numeric(rownames(adjacency_nn)),each=length(rownames(adjacency_nn)))
#   #                 el_nn_with_dist <- cbind(first,second,as.vector(adjacency_nn))
#   #         dummy_graph <- graph.empty(n=nrow(cluster_distances_matrix), directed=FALSE)
#   #         dummy_graph <- graph.empty(n=nrow(adjacency_nn), directed=FALSE)
#   
#   #           dummy_graph <- graph.empty(n=max(as.numeric(rownames(cluster_distances_matrix))), directed=FALSE)
#   #           dummy_graph <- add.edges(dummy_graph,edges=as.vector(t(el_nn_with_dist[,1:2])))
#   
#   clu <- clusters(dummy_graph)
#   # get the cluster that each node belongs to
#   members <- clu$membership
#   # initialize matrix to hold distances between graph components
#   linkweights <- matrix(nrow=clu$no,ncol=clu$no)
#   
#   # loop through all components vs. all components to fill in distances
#   for (i in 1:clu$no) {
#     members_i <- which(members==i)
#     for (j in 1:clu$no) {
#       members_j <- which(members==j)
#       # get the minimum distance between components i and j
#       # (note members has zero index)
#       tmpmatrix <- cluster_distances_matrix[members_i,members_j]
#       
#       #             tmpmatrix <- cluster_distances_matrix[as.character(members_i),as.character(members_j)]
#       linkweights[i,j] <- min(tmpmatrix)
#     }
#   }
#   
#   # set i-i distances to Inf instead of 0 so they aren't picked
#   # as the closest neighbors
#   
#   # linkweights
#   
#   for (i in 1:ncol(linkweights)) {
#     linkweights[i,i] <- Inf
#   }
#   
#   # make the minimum spanning tree for graph components
#   components_adjacency <- graph.adjacency(linkweights,mode="undirected",weighted=TRUE,diag=TRUE)
#   #         components_adjacency <- graph.adjacency(adjacency_nn,mode="undirected",weighted=TRUE,diag=TRUE)
#   components_mst <- minimum.spanning.tree(components_adjacency)
#   components_mst_el <- get.edgelist(components_mst,names=FALSE)
#   
#   for (i in 1:nrow(components_mst_el)) {
#     # get index of shortest connection between components
#     # tmp <- which.min(adjacency[which(members==components_mst_el[i,1]),which(members==components_mst_el[i,2])])
#     members_x <- which(members==components_mst_el[i,1])
#     for (x in 1:length(members_x)) {
#       members_x[x] <- indtoclust[x][[1]]
#     }
#     #convertx <- indtoclust[members_x][[1]]
#     members_y <- which(members==components_mst_el[i,2])
#     for (y in 1:length(members_x)) {
#       members_y[y] <- indtoclust[y][[1]]
#     }
#     #converty <- indtoclust[members_y][[1]]
#     tmpmatrix <- as.matrix(cluster_distances_matrix[members_x,members_y])
#     if(length(members_x)==1){
#       tmpmatrix <- t(tmpmatrix)
#     }
#     tmp_min <- min(tmpmatrix) # get the largest weight = shortest distance
#     tmp_index <- which(tmpmatrix==tmp_min,arr.ind=TRUE)
#     
#     # add new edge to graph with the shortest connection vertices and the weight from linkweights
#     if (are.connected(output_graph,members_x[tmp_index[1]],members_y[tmp_index[2]])) {
#       # cat(paste("\ntrue = ",i))
#       E(output_graph,P=c(members_x[tmp_index[1]],members_y[tmp_index[2]]))$label <- "MST"
#     } else {
#       output_graph <- add.edges(output_graph,c(members_x[tmp_index[1]],
#                                                members_y[tmp_index[2]]),weight=tmp_min,label="MST",sequence_assignment=n+1)
#     }
#   }
# }
# # convert graph distances to weights (low distance = high weight and vice versa)
# distances <- E(output_graph)$weight
# #weights <- -distances+max(distances)+min(distances)
# weights <- 1/distances
# E(output_graph)$weight <- weights
# 
# ################################################################################
# # This section annotates the graph
# ################################################################################
# 
# # iterate through all times and annotate
# 
# anno_cat <- c()
# 
# for (treat in INDIV_DIR_NAMES) {
#   treat_anno_cat <- c()
#   treat_anno <- list()
#   
#   for (f in 1:length(fcs_files[[treat]])) {
#     
#     # get fcs file corresponding to time[i]
#     fcs_file = fcs_files[[treat]][f]
#     
#     currentfile <- tail(strsplit(fcs_file,"/")[[1]],n=1)
#     cat("Annotating graph for",currentfile,"\n")
#     
#     # get medians (and count) for all parameters
#     treat_anno$count  <- cluster_counts[[treat]][[f]]
#     treat_anno$percenttotal <- data.frame(percenttotal = c(cluster_counts[[treat]][[f]]$counts/SUBSAMPLE))
#     treat_anno$medians  <- cluster_medians[[treat]][[f]]
#     
#     # add time[f] information column
#     time_matrix <- matrix(f,nrow=length(treat_anno$count$counts))
#     colnames(time_matrix) <- c("Timepoint")
#     name_column <- matrix(data=treat,nrow=nrow(treat_anno$medians),ncol=1,byrow=FALSE,dimnames=list(NULL,"Treatment"))
#     treat_anno$medians <- cbind(treat_anno$medians,time_matrix,name_column)
#     
#     # add median and percent values
#     for (a in c("medians", "percenttotal")) {
#       treat_anno_cat[[a]] <- rbind(treat_anno_cat[[a]],treat_anno[[a]])
#     }
#   }
#   
#   anno_cat <- rbind(anno_cat, treat_anno_cat)
#   rm(treat_anno_cat)
#   rm(treat_anno)
#   # combine anno_cat matrices
# }
# 
# output_anno <- c()
# 
# for (t in 1:dim(anno_cat)[1]) {
#   tmp <- cbind(anno_cat[,1][[t]], anno_cat[,2][[t]])
#   output_anno <- rbind(output_anno,tmp)
# }
# 
# for (c in colnames(output_anno)) {
#   output_graph <- set.vertex.attribute(output_graph,c,index=as.numeric(rownames(output_anno)), value=output_anno[,c])
# }


