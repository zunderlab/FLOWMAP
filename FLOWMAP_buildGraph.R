
findNormalized <- function(cluster_distances_matrix, per, min,
                           max, numcluster, table_lengths = FALSE,
                           table_breaks, offset) {
  # make fully connected graph from adjacency list
  # note: "weight" here is really distance, calling it weight for the mst function later needs this
  # cat("table_lengths are", table_lengths, "\n")
  if (table_lengths[1] == FALSE) {
    # print("a")
    edgelist_with_distances <- cbind(as.vector(row(cluster_distances_matrix)),
                                     as.vector(col(cluster_distances_matrix)),
                                     as.vector(cluster_distances_matrix))
    # take only the upper triangle of the matrix
    # print("head of edgelist_with_distances")
    # print(head(edgelist_with_distances))
    # cat("dim of edgelist_with_distances", dim(edgelist_with_distances), "\n")
    edgelist_with_distances <- edgelist_with_distances[upper.tri(matrix(data = 1:length(cluster_distances_matrix),
                                                                        nrow = nrow(cluster_distances_matrix),
                                                                        ncol = ncol(cluster_distances_matrix))),]
    # print("head of edgelist_with_distances")
    # print(head(edgelist_with_distances))
    # cat("dim of edgelist_with_distances", dim(edgelist_with_distances), "\n")
    # convert strings to numeric
    class(edgelist_with_distances) <- "numeric"
    # print("head of edgelist_with_distances")
    # print(head(edgelist_with_distances))
    # cat("dim of edgelist_with_distances", dim(edgelist_with_distances), "\n")
    # sort edges by distance (column 3)
    edgelist_with_distances <- edgelist_with_distances[order(edgelist_with_distances[, 3]), ]
    # print("head of edgelist_with_distances")
    # print(head(edgelist_with_distances))
    # cat("dim of edgelist_with_distances", dim(edgelist_with_distances), "\n")
    val <- max(floor(length(edgelist_with_distances[, 1]) * per / 100), 1)
    trim_edgelist_with_distances <- edgelist_with_distances[(1:val), ]
    # print("trim_edgelist_with_distances")
    # print(trim_edgelist_with_distances)
    # cat("dim of trim_edgelist_with_distances", dim(trim_edgelist_with_distances), "\n")
    # print(head(trim_edgelist_with_distances))
    # calculate "density" for each cluster and normalize
    trim_edgelist_with_distances <- as.data.frame(trim_edgelist_with_distances)
    # print("trim_edgelist_with_distances")
    # print(trim_edgelist_with_distances)
    # cat("dim of trim_edgelist_with_distances", dim(trim_edgelist_with_distances), "\n")
    # trim_edgelist_with_distances <- t(trim_edgelist_with_distances)
    # cat("dim of trim_edgelist_with_distances", dim(trim_edgelist_with_distances), "\n")
    # print("trim_edgelist_with_distances[, 1:2]")
    # print(head(trim_edgelist_with_distances[, 1:2]))
    densities_no_zeros <- table(trim_edgelist_with_distances[, 1:2])
    # print("densities_no_zeros")
    # print(densities_no_zeros)
    # cat("length of densities_no_zeros", length(densities_no_zeros), "\n")
    # cat("dim of densities_no_zeros", dim(densities_no_zeros), "\n")
    # add in zeros for clusters with no edges (table function leaves these out)
    densities <- rep(0, numcluster)
    # print("densities")
    # print(densities)
    names(densities) <- (1:numcluster)
    # cat("dim of densities_no_zeros", dim(densities_no_zeros), "\n")
    # cat("names of densities_no_zeros", names(densities_no_zeros), "\n")
    # print("table(densities_no_zeros)")
    # print(table(densities_no_zeros))
    
    # which(boop != 0, arr.ind = TRUE)
    # as.numeric(rownames(which(boop != 0, arr.ind = TRUE)))
    # table(as.numeric(rownames(which(boop != 0, arr.ind = TRUE))))
    table_ind <- which(densities_no_zeros != 0, arr.ind = TRUE)
    densities_no_zeros <- table(as.numeric(rownames(table_ind)))
    densities[as.numeric(names(densities_no_zeros))] <- densities_no_zeros
    # densities[as.numeric(names(densities_no_zeros))] <- densities_no_zeros
    # print("densities")
    # print(densities)
    normalized_densities <- round(densities/max(densities) * (max - min) + min)
    # print("normalized_densities")
    # print(normalized_densities)
  }
  else {
    # print("b")
    n <- 1
    # take the bottom half of the matrix, giving (n to n+1) and (n+1 to n+1) distances
    # as separate matrices
    n_n1 <- cluster_distances_matrix[(table_lengths[n] + 1):(table_lengths[n] + table_lengths[n + 1]),
                                     1:table_lengths[n]]
    n1_n1 <- cluster_distances_matrix[(table_lengths[n] + 1):(table_lengths[n] + table_lengths[n + 1]),
                                      (table_lengths[n] + 1):(table_lengths[n] + table_lengths[n + 1])]
    # convert the distance matrices to lists of coordinates with distances as the 3rd row
    edgelist_with_distances_n_n1 <- cbind(as.vector(row(n_n1) + table_breaks[n]),
                                          as.vector(col(n_n1) + table_breaks[n] - table_lengths[n]),
                                          as.vector(n_n1))
    # print("head of edgelist_with_distances_n_n1")
    # print(head(edgelist_with_distances_n_n1))
    # cat("dim of edgelist_with_distances_n_n1", dim(edgelist_with_distances_n_n1), "\n")
    edgelist_with_distances_n1_n1 <- cbind(as.vector(row(n1_n1) + table_breaks[n]),
                                           as.vector(col(n1_n1) + table_breaks[n]),
                                           as.vector(n1_n1))
    # print("head of edgelist_with_distances_n1_n1")
    # print(head(edgelist_with_distances_n1_n1))
    # cat("dim of edgelist_with_distances_n1_n1", dim(edgelist_with_distances_n1_n1), "\n")    # only take upper triangle from n+1 to n+1 matrix
    edgelist_with_distances_n1_n1 <- edgelist_with_distances_n1_n1[upper.tri(matrix(data = 1:length(n1_n1),
                                                                                    nrow = nrow(n1_n1),
                                                                                    ncol = ncol(n1_n1))),]
    # print("edgelist_with_distances_n1_n1")
    # print(edgelist_with_distances_n1_n1)
    # cat("dim of edgelist_with_distances_n1_n1", dim(edgelist_with_distances_n1_n1), "\n")    # only take upper triangle from n+1 to n+1 matrix
    #combine n_n1 and n1_n1 distances into a single "edgelist" with distances
    edgelist_with_distances <- rbind(edgelist_with_distances_n_n1, edgelist_with_distances_n1_n1)
    # print("head of edgelist_with_distances")
    # print(head(edgelist_with_distances))
    # cat("dim of edgelist_with_distances", dim(edgelist_with_distances), "\n")
    # convert strings to numeric
    class(edgelist_with_distances) <- "numeric"
    # sort edges by distance (column 3)
    edgelist_with_distances <- edgelist_with_distances[order(edgelist_with_distances[, 3]), ]
    # print("head of edgelist_with_distances")
    # print(head(edgelist_with_distances))
    # cat("dim of edgelist_with_distances", dim(edgelist_with_distances), "\n")
    val <- max(floor(length(edgelist_with_distances[, 1]) * per / 100), 1)
    # take only the top X-% edges by distance based on PERCENT_TOTAL
    trim_edgelist_with_distances <- edgelist_with_distances[(1:val), ]
    # print("trim_edgelist_with_distances")
    # print(trim_edgelist_with_distances)
    # cat("dim of trim_edgelist_with_distances", dim(trim_edgelist_with_distances), "\n")
    trim_edgelist_with_distances <- as.data.frame(trim_edgelist_with_distances)
    # print("trim_edgelist_with_distances")
    # print(trim_edgelist_with_distances)
    # cat("dim of trim_edgelist_with_distances", dim(trim_edgelist_with_distances), "\n")
    # trim_edgelist_with_distances <- t(trim_edgelist_with_distances)
    # cat("dim of trim_edgelist_with_distances", dim(trim_edgelist_with_distances), "\n")
    # calculate "density" for each cluster and normalize
    # print("trim_edgelist_with_distances")
    # print(trim_edgelist_with_distances)
    # print("trim_edgelist_with_distances[, 1:2]")
    # print(head(trim_edgelist_with_distances[, 1:2]))
    densities_no_zeros <- table(trim_edgelist_with_distances[, 1:2])
    # print("densities_no_zeros")
    # print(densities_no_zeros)
    # print("table(densities_no_zeros)")
    # print(table(densities_no_zeros))
    # cat("length of densities_no_zeros", length(densities_no_zeros), "\n")
    # cat("dim of densities_no_zeros", dim(densities_no_zeros), "\n")
    # add in zeros for clusters with no edges (table function leaves these out)
    densities <- rep(0, numcluster)
    # print("densities")
    # print(densities)
    names(densities) <- (offset + 1):(offset + table_lengths[n] + table_lengths[n + 1])
    # cat("names of densities_no_zeros", names(densities_no_zeros), "\n")
    table_ind <- which(densities_no_zeros != 0, arr.ind = TRUE)
    densities_no_zeros <- table(as.numeric(rownames(table_ind)))
    # print("densities_no_zeros")
    # print(densities_no_zeros)
    densities[names(densities_no_zeros)] <- densities_no_zeros
    # densities[names(densities_no_zeros)] <- densities_no_zeros
    # densities[names(densities_no_zeros)] <- densities_no_zeros
    # print("densities")
    # print(densities)
    normalized_densities <- round(densities / max(densities) * (max - min) + min)
    # print("normalized_densities")
    # print(normalized_densities)
  }
  return(normalized_densities)
}


drawNormalizedEdges <- function(output_graph, cluster_distances_matrix,
                                normalized_densities, offset = FALSE, n) {
  final_edgelist_with_distances <- c()
  # print("normalized_densities")
  # print(normalized_densities)
  # cat("dim of cluster_distances_matrix", dim(cluster_distances_matrix), "\n")
  if (!offset) {
    # print("woof")
    # cat("length of normalized_densities", length(normalized_densities), "\n")
    for (i in 1:length(normalized_densities)) {
      tmp_final_edgelist_with_distances <- cbind(i, order(cluster_distances_matrix[, i]),
                                                 sort(cluster_distances_matrix[, i]))[1:normalized_densities[i], ]
      final_edgelist_with_distances <- rbind(final_edgelist_with_distances,
                                             tmp_final_edgelist_with_distances)
    }
  }
  else {
    # print("meow")
    # cat("length of normalized_densities", length(normalized_densities), "\n")
    for (i in 1:length(normalized_densities)) {
      tmp_final_edgelist_with_distances <- cbind((i + offset), (order(cluster_distances_matrix[, i]) + offset),
                                                 sort(cluster_distances_matrix[, i]))[1:normalized_densities[i], ]
      final_edgelist_with_distances <- rbind(final_edgelist_with_distances,
                                             tmp_final_edgelist_with_distances)
    }
  }
  # remove duplicate edges from edgelist
  final_edgelist_with_distances <- unique(cbind(t(apply(final_edgelist_with_distances[, 1:2], 1, sort)),
                                                final_edgelist_with_distances[, 3]))
  # add edges to graph
  # note: weight here is really distance, will be converted to weight after the graph is completed
  output_graph <- add.edges(output_graph, edges = as.vector(t(final_edgelist_with_distances[, 1:2])),
                            weight = final_edgelist_with_distances[, 3],
                            label = "EXTRA", sequence_assignment = n)
  return(output_graph)
}


checkMSTedges <- function(output_graph, cluster_distances_matrix,
                          table_lengths, offset, n) {
  # create dummy graph contains n and n+1 vertices, but n_n section is fully connected so
  # it will act as a single connected component while n_n+1 and n+1_n+1 sections are not
  # connected, so each n+1 vertex will be its own disconnected component
  adjacency_nn <- cluster_distances_matrix[1:table_lengths[1], 1:table_lengths[1]]
  el_nn_with_dist <- cbind(as.vector(row(adjacency_nn)),
                           as.vector(col(adjacency_nn)), 
                           as.vector(adjacency_nn))
  dummy_graph <- graph.empty(n = nrow(cluster_distances_matrix), directed = FALSE)
  dummy_graph <- add.edges(dummy_graph, edges = as.vector(t(el_nn_with_dist[, 1:2])))
  clu <- igraph::clusters(dummy_graph)
  # get the cluter that each node belongs to
  members <- clu$membership
  # initialize matrix to hold distances between graph components
  linkweights <- matrix(nrow = clu$no, ncol = clu$no)
  # loop through all components vs. all components to fill in distances
  for (i in 1:clu$no) {
    members_i <- which(members == i)
    for (j in 1:clu$no) {
      members_j <- which(members == j)
      # get the minimum distance between components i and j
      # (note members has zero index)
      tmpmatrix <- cluster_distances_matrix[members_i, members_j]
      linkweights[i, j] <- min(tmpmatrix)
    }
  }
  # set i-i distances to Inf instead of 0 so they aren't picked
  # as the closest neighbors
  for (i in 1:ncol(linkweights)) {
    linkweights[i, i] <- Inf
  }
  # make the minimum spanning tree for graph components
  components_adjacency <- graph.adjacency(linkweights, mode = "undirected",
                                          weighted = TRUE, diag = TRUE)
  components_mst <- minimum.spanning.tree(components_adjacency)
  components_mst_el <- get.edgelist(components_mst, names = FALSE)
  # go from MST edges to the actual edges between individual nodes
  for (i in 1:nrow(components_mst_el)) {
    # get index of shortest connection between components
    members_x <- which(members == components_mst_el[i, 1])
    members_y <- which(members == components_mst_el[i, 2])
    tmpmatrix <- as.matrix(cluster_distances_matrix[members_x, members_y])
    if(length(members_x) == 1){
      tmpmatrix <- t(tmpmatrix)
    }
    tmp_min <- min(tmpmatrix) # get the largest weight = shortest distance
    tmp_index <- which(tmpmatrix == tmp_min, arr.ind = TRUE)
    # add new edge to graph with the shortest connection vertices and the weight from linkweights
    if (are.connected(output_graph, members_x[tmp_index[1]] + offset, 
                      members_y[tmp_index[2]] + offset)) {
      E(output_graph,P = c(members_x[tmp_index[1]] + offset,
                           members_y[tmp_index[2]] + offset))$label <- "MST"
    } else {
      output_graph <- add.edges(output_graph, c(members_x[tmp_index[1]] + offset,
                                                members_y[tmp_index[2]] + offset),
                                weight = tmp_min, label = "MST",
                                sequence_assignment = n)
    }
  }
  return(output_graph)
}


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
  for (t in 1:timepoints) {
    temp_medians <- data.frame()
    temp_cellassgn <- data.frame()
    temp_counts <- data.frame()
    for (treat in treatments) {
      temp_medians <- rbind(temp_medians, listOfFLOWMAPclusters[[treat]]$cluster_medians[[t]])
      temp_cellassgn <- rbind(temp_cellassgn, listOfFLOWMAPclusters[[treat]]$cellassgn[[t]])
      temp_counts <- rbind(temp_counts, listOfFLOWMAPclusters[[treat]]$cluster_counts[[t]])
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


initializeGraph <- function(FLOWMAPclusters) {
  # This section initializes the graph
  total_nodes <- length(FLOWMAPclusters$fullclusters[, 1])
  # create empty graph with the right number of nodes - will fill in the edges later
  initial_graph <- graph.empty(n = total_nodes, directed = FALSE)
  return(initial_graph)
}

buildFirstFLOWMAP <- function(FLOWMAPclusters, per, min, max, distance_metric, clustering.var) {
  n <- 1
  output_graph <- initializeGraph(FLOWMAPclusters)
  # This section creates a flowmap for the first time point
  cat("Building first flowmap\n")
  table_lengths <- FLOWMAPclusters$table_lengths
  # get distance matrix from clusters
  clusters <- FLOWMAPclusters$fullclusters[1:table_lengths[1], ]
  clusters <- subset(clusters, select = clustering.var)
  cluster_distances <- dist(clusters, method = distance_metric, diag = TRUE, upper = TRUE)
  cluster_distances_matrix <- as.matrix(cluster_distances)
  # set i-i distances to Inf instead of 0 so they aren't the closest neighbors
  for (i in 1:ncol(cluster_distances_matrix)) {
    cluster_distances_matrix[i, i] <- Inf
  }
  numcluster <- nrow(clusters)
  normalized_densities <- findNormalized(cluster_distances_matrix,
                                         per, min, max, numcluster)
  # cat("length of normalized_densities is", length(normalized_densities), "\n")
  # print("normalized_densities")
  # print(normalized_densities)
  # build new edgelist with N edges for each cluster based on normalized density
  # cat("Building first edgelist\n")
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
    if (are.connected(output_graph, mst_graph_edgelist[i, 1], mst_graph_edgelist[i, 2])) {
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


buildFLOWMAP <- function(FLOWMAPclusters, per, min, max, distance_metric, cellnum, clustering.var) {
  output_graph <- buildFirstFLOWMAP(FLOWMAPclusters, per, min, max, distance_metric)
  table_breaks <- c(0, FLOWMAPclusters$table_breaks)
  # This section builds the flowmap one timepoint at a time
  for (n in 1:(length(FLOWMAPclusters$cluster_medians) - 1)) {
    # offset value is used to correctly index the edges at each sequential step
    offset <- table_breaks[n]
    # go through sequential cluster sets, add edges for each a-a and a-a+1 set, and also label sequential MST
    cat("Build FlowMap from", n, "to", n + 1, "\n")
    # get clusters for time a and a+1
    clusters <- rbind(FLOWMAPclusters$cluster_medians[[n]], FLOWMAPclusters$cluster_medians[[n + 1]])
    clusters <- subset(clusters, select = clustering.var)
    numcluster <- nrow(clusters)
    # make adjacency matrix from clusters
    # if (distance_metric == "cosine") {
    #   cluster_distances <- cosine_similarity_matrix(clusters)
    # }
    # else {
    cluster_distances <- dist(clusters, method = distance_metric, diag = TRUE, upper = TRUE)
    # }
    cluster_distances_matrix <- as.matrix(cluster_distances)
    # set i-i distances to Inf instead of 0 so they aren't the closest neighbors
    for (i in 1:ncol(cluster_distances_matrix)) {
      cluster_distances_matrix[i, i] <- Inf
    }
    # This section adds the lowest distance n_n+1 and n+1_n+1 edges to the output graph
    n_n1_table_lengths <- FLOWMAPclusters$table_lengths[n:(n + 1)]
    n_n1_table_breaks <- FLOWMAPclusters$table_breaks[n:(n + 1)]
    # print(dim(cluster_distances_matrix))
    normalized_densities <- findNormalized(cluster_distances_matrix, per, min,
                                           max, numcluster, table_lengths = n_n1_table_lengths,
                                           table_breaks = n_n1_table_breaks,
                                           offset = offset)
    # cat("length of normalized_densities is", length(normalized_densities), "\n")
    # print(head(normalized_densities))
    # print(tail(normalized_densities))
    # build new edgelist with N edges for each cluster based on normalized density
    # cat("Building edgelist\n")
    # cat("dim of cluster_distances_matrix", dim(cluster_distances_matrix), "\n")
    output_graph <- drawNormalizedEdges(output_graph, cluster_distances_matrix,
                                        normalized_densities, offset = offset, n = n)
    # This section adds the "MST" for n_n+1 and n+1_n+1 nodes
    output_graph <- checkMSTedges(output_graph, cluster_distances_matrix, 
                                  n_n1_table_lengths, offset = offset, n = n)
  }
  # convert graph distances to weights (low distance = high weight and vice versa)
  distances <- E(output_graph)$weight
  #weights <- -distances+max(distances)+min(distances)
  weights <- (1 / distances)
  E(output_graph)$weight <- weights
  output_graph <- annotateGraph(output_graph, FLOWMAPclusters, cellnum = cellnum)
  return(output_graph)
}


buildSimpleFLOWMAP <- function(FLOWMAPclusters, top_npercent, distance_metric, cellnum, clustering.var) {
  # SimpleFLOWMAP does not enforce minimum or maximum edges, just 
  # picks the top n% of edges for tn to tn and tn to tn+1
  output_graph <- initializeGraph(FLOWMAPclusters)
  table_breaks <- c(0, FLOWMAPclusters$table_breaks)
  # This section builds the flowmap one timepoint at a time
  for (n in 1:(length(FLOWMAPclusters$cluster_medians) - 1)) {
    # offset value is used to correctly index the edges at each sequential step
    offset <- table_breaks[n]
    # go through sequential cluster sets, add edges for each a-a and a-a+1 set, and also label sequential MST
    cat("Build FlowMap from", n, "to", n + 1, "\n")
    # get clusters for time a and a+1
    clusters <- rbind(FLOWMAPclusters$cluster_medians[[n]], FLOWMAPclusters$cluster_medians[[n + 1]])
    clusters <- subset(clusters, select = clustering.var)
    numcluster <- nrow(clusters)
    # make adjacency matrix from clusters
    # if (distance_metric == "cosine") {
    #   cluster_distances <- cosine_similarity_matrix(clusters)
    # }
    # else {
    cluster_distances <- dist(clusters, method = distance_metric, diag = TRUE, upper = TRUE)
    # }
    cluster_distances_matrix <- as.matrix(cluster_distances)
    # set i-i distances to Inf instead of 0 so they aren't the closest neighbors
    for (i in 1:ncol(cluster_distances_matrix)) {
      cluster_distances_matrix[i, i] <- Inf
    }
    # This section adds the lowest distance n_n+1 and n+1_n+1 edges to the output graph
    n_n1_table_lengths <- FLOWMAPclusters$table_lengths[n:(n + 1)]
    n_n1_table_breaks <- FLOWMAPclusters$table_breaks[n:(n + 1)]
    
    # build edgelist for tn to tn and tn to tn+1 with top n% of edges
    cat("Building edgelist\n")
    output_graph <- drawTopEdges(output_graph, cluster_distances_matrix,
                                 top_npercent, offset = offset, n = n)
    #     print(is.igraph(output_graph))
    #     # add MST to tn and tn and tn to tn+1 so no nodes are lost (no edges)
    #     output_graph <- checkMSTedges(output_graph, cluster_distances_matrix,
    #                                   n_n1_table_lengths, offset = offset, n = n)
    #     print(is.igraph(output_graph))
  }
  clusters <- FLOWMAPclusters$fullclusters
  clusters <- subset(clusters, select = clustering.var)
  numcluster <- nrow(clusters)
  # if (distance_metric == "cosine") {
  #   cluster_distances <- cosine_similarity_matrix(clusters)
  # }
  # else {
  cluster_distances <- dist(clusters, method = distance_metric, diag = TRUE, upper = TRUE)
  # }
  cluster_distances_matrix <- as.matrix(cluster_distances)
  
  #   print(is.igraph(output_graph))
  # convert graph distances to weights (low distance = high weight and vice versa)
  distances <- E(output_graph)$weight
  #weights <- -distances+max(distances)+min(distances)
  weights <- 1/distances
  E(output_graph)$weight <- weights
  output_graph <- annotateGraph(output_graph, FLOWMAPclusters, cellnum)
  return(output_graph)
}


annotateGraph <- function(output_graph, FLOWMAPclusters, cellnum, ...) {
  # This section annotates the graph
  anno_cat <- c()
  anno <- list()
  # iterate through all times and annotate
  for (f in 1:length(FLOWMAPclusters$cluster_medians)) {
    cat("Annotating graph for file", f, "\n")
    # get medians for all parameters and counts for all clusters
    counts <- FLOWMAPclusters$cluster_counts[[f]]$Counts
    anno$count <- counts
    anno$percenttotal <- data.frame(percenttotal = c(counts / cellnum))
    anno$medians  <- FLOWMAPclusters$cluster_medians[[f]]
    # add time information column
    time_matrix <- matrix(f, nrow = length(anno$count))
    colnames(time_matrix) <- c("Timepoint")
    anno$medians <- cbind(anno$medians, time_matrix)
    # add median and percent values
    for (a in c("medians", "percenttotal")) {
      # anno_cat is all the annos concatenated, will be
      # used to make "anno.Rsave" file
      anno_cat[[a]] <- rbind(anno_cat[[a]], anno[[a]])
    }
  }
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


