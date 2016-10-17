
ConvertIndex <- function(inds, orig.data.frame, keep.names = TRUE) {
  # inds go down through the rows of the first column
  # before starting at the first row of the second column
  # and going back down, etc.
  col.inds <- ceiling(inds / nrow(orig.data.frame))
  row.inds <- (inds %% nrow(orig.data.frame))
  row.inds[which(row.inds == 0)] <- nrow(orig.data.frame)
  if (keep.names == TRUE) {
    row.inds <- as.numeric(rownames(orig.data.frame)[row.inds])
    col.inds <- as.numeric(colnames(orig.data.frame)[col.inds])
  }
  arr.inds <- cbind(row.inds, col.inds)
  return(arr.inds)
}

MergeValues <- function(arr.inds, orig.data.frame) {
  merged <- cbind(arr.inds, values = sort(as.matrix(orig.data.frame)))
  return(merged)
}

removeDuplicateValues <- function(arr.inds.with.values) {
  x <- arr.inds.with.values[!duplicated(arr.inds.with.values[, "values"]),]
  return(x)
}

RemoveWithinNEdges <- function(arr.inds.with.values, inds.in.n) {
  d <- arr.inds.with.values
  d <- d[!(!is.na(match(d[, "row.inds"], inds.in.n)) &
             !is.na(match(d[, "col.inds"], inds.in.n))), ]
  inds.no.n_n <- d
  return(inds.no.n_n)
}

FindNormalized <- function(cluster.distances.matrix, per, min,
                           max, numcluster, table.lengths,
                           table.breaks, offset) {
  # make fully connected graph from adjacency list
  # note: "weight" here is really distance, calling it weight for the mst function later needs this
  arr.inds <- ConvertIndex(inds = order(cluster.distances.matrix), orig.data.frame = cluster.distances.matrix, keep.names = TRUE)
  edgelist.with.distances <- MergeValues(arr.inds = arr.inds, orig.data.frame = cluster.distances.matrix)
  edgelist.with.distances <- removeDuplicateValues(edgelist.with.distances)
  # remove last edge with Inf weight
  edgelist.with.distances <- edgelist.with.distances[1:(nrow(edgelist.with.distances) - 1), ]
  # print("nrow(edgelist.with.distances)")
  # print(nrow(edgelist.with.distances))
  # print("head(edgelist.with.distances)")
  # print(head(edgelist.with.distances))
  if (table.lengths[1] != FALSE) {
    inds.in.n <- (offset + 1):(offset + table.lengths[1])
    # cat("inds.in.n is", inds.in.n, "\n")
    edgelist.with.distances <- RemoveWithinNEdges(edgelist.with.distances,
                                                  inds.in.n = inds.in.n)
  }
  # print("nrow(edgelist.with.distances)")
  # print(nrow(edgelist.with.distances))
  # print("head(edgelist.with.distances)")
  # print(head(edgelist.with.distances))
  # val <- max(floor(length(edgelist.with.distances[, 1]) * per / 100), 1)
  num.edges <- length(edgelist.with.distances[, 1])
  val <- max(floor(num.edges * per / 100), 1)
  # val <- max(floor(choose(numcluster, 2) * per / 100), 1)
  # print("choose(numcluster, 2)")
  # print(choose(numcluster, 2))
  # cat("val is", val, "\n")
  trim.edgelist.with.distances <- edgelist.with.distances[(1:val), ]
  # calculate "density" for each cluster and normalize
  if (val == 1) {
    trim.edgelist.with.distances <- t(trim.edgelist.with.distances)
  }
  # print("trim.edgelist.with.distances")
  # print(trim.edgelist.with.distances)
  densities.no.zeros <- table(trim.edgelist.with.distances[, 1:2])
  # print("densities.no.zeros")
  # print(densities.no.zeros)
  # add in zeros for clusters with no edges (table function leaves these out)
  densities <- rep(0, numcluster)
  if (table.lengths[1] == FALSE) {
    names(densities) <- (1:numcluster)
  } else {
    names(densities) <- (offset + 1):(offset + table.lengths[1] + table.lengths[2])
  }
  # table.ind <- which(densities.no.zeros != 0, arr.ind = TRUE)
  # print("table.ind")
  # print(table.ind)
  # densities.no.zeros <- table(as.numeric(rownames(table.ind)))
  # print("densities.no.zeros")
  # print(densities.no.zeros)
  densities[names(densities.no.zeros)] <- densities.no.zeros
  # print("densities")
  # print(densities)
  normalized.densities <- round(densities/max(densities) * (max - min) + min)
  # print("normalized.densities")
  # print(normalized.densities)
  return(list(normalized.densities = normalized.densities,
              edgelist.with.distances = edgelist.with.distances))
  
  # edgelist_with_distances_n_n1 <- cbind(as.vector(row(n_n1)+table_breaks[n]),
  #                                       as.vector(col(n_n1)+table_breaks[n]-table_lengths[n]),as.vector(n_n1))
  # edgelist_with_distances_n1_n1 <- cbind(as.vector(row(n1_n1)+table_breaks[n]),
  #                                        as.vector(col(n1_n1)+table_breaks[n]),as.vector(n1_n1))
  # # only take upper triangle from n+1 to n+1 matrix
  # edgelist_with_distances_n1_n1 <- edgelist_with_distances_n1_n1[upper.tri(matrix(data=1:length(n1_n1),
  #                                                                                 nrow=nrow(n1_n1),ncol=ncol(n1_n1))),]
  # #combine n_n1 and n1_n1 distances into a single "edgelist" with distances
  # edgelist_with_distances <- rbind(edgelist_with_distances_n_n1,edgelist_with_distances_n1_n1)
  # # convert strings to numeric
  # class(edgelist_with_distances) <- "numeric"
  # # sort edges by distance (column 3)
  # edgelist_with_distances <- edgelist_with_distances[order(edgelist_with_distances[,3]),]
  # # print("edgelist_with_distances")
  # # print(edgelist_with_distances)
  # # # take only the top X-% edges by distance based on PERCENT_TOTAL
  # trim_edgelist_with_distances <- edgelist_with_distances[1:floor(length(edgelist_with_distances[,1])*PERCENT_EDGES/100),]
  # print("trim_edgelist_with_distances")
  # print(trim_edgelist_with_distances)
  # # add edges to graph
  # # note: weight here is really distance, will be converted to weight after the graph is completed
  # 
  # # calculate "density" for each cluster and normalize
  # densities_no_zeros <- table(trim_edgelist_with_distances[,1:2])
  # print("densities_no_zeros")
  # print(densities_no_zeros)
  # # add in zeros for clusters with no edges (table function leaves these out)
  # densities <- rep(0,nrow(clusters))
  # names(densities) <- (offset+1):(offset+table_lengths[n]+table_lengths[n+1])
  # densities[names(densities_no_zeros)] <- densities_no_zeros
  # normalized_densities <- round(densities/max(densities)*(MAX_EDGES - MIN_EDGES) + MIN_EDGES)
  # print("normalized_densities")
  # print(normalized_densities)
  # 
  # # build new edgelist with N edges for each cluster based on normalized density
  # cat("building edgelist\n")
  # final_edgelist_with_distances <- c()
  # for (i in 1:length(normalized_densities)) {
  #   tmp_final_edgelist_with_distances <- cbind(offset+i,order(cluster_distances_matrix[,i])+offset,
  #                                              sort(cluster_distances_matrix[,i]))[1:normalized_densities[i],]
  #   final_edgelist_with_distances <- rbind(final_edgelist_with_distances,tmp_final_edgelist_with_distances)
  # }
  
}



DrawNormalizedEdges <- function(output.graph, cluster.distances.matrix,
                                normalized.densities, n, offset = FALSE) {
  final.edgelist.with.distances <- c()
  # print("names(normalized.densities)")
  # print(names(normalized.densities))
  # print("cluster.distances.matrix")
  # print(cluster.distances.matrix)
  # print("names(normalized.densities)")
  # print(names(normalized.densities))
  for (i in names(normalized.densities)) {
    # print("i")
    # print(i)
    matches.in.order <- order(cluster.distances.matrix[, i])
    # print("order(cluster.distances.matrix[, i])")
    # print(order(cluster.distances.matrix[, i]))
    # print("rownames(cluster.distances.matrix)[matches.in.order]")
    # print(rownames(cluster.distances.matrix)[matches.in.order])
    tmp.edgelist <- cbind(as.numeric(i), as.numeric(rownames(cluster.distances.matrix)[matches.in.order]),
                                               sort(cluster.distances.matrix[, i]))[1:normalized.densities[i], ]
    final.edgelist.with.distances <- rbind(final.edgelist.with.distances, tmp.edgelist)
  }
  colnames(final.edgelist.with.distances) <- c("row.inds", "col.inds", "values")
  final.edgelist.with.distances <- removeDuplicateValues(final.edgelist.with.distances)
  # note: weight here is really distance, will be converted to weight after the graph is completed
  # print("dim(final.edgelist.with.distances)")
  # print(dim(final.edgelist.with.distances))
  # print("head(final.edgelist.with.distances)")
  # print(head(final.edgelist.with.distances))
  # print("final.edgelist.with.distances")
  # print(final.edgelist.with.distances)
  vertices.edges <- as.vector(t(as.matrix(final.edgelist.with.distances[, 1:2])))
  output.graph <- add.edges(output.graph, edges = vertices.edges,
                            weight = final.edgelist.with.distances[, 3],
                            label = "EXTRA", sequence_assignment = n)
  # return(output.graph)
  return(list(output.graph = output.graph,
              final.edgelist.with.distances = final.edgelist.with.distances))
}


# DrawNormalizedEdges <- function(output.graph, edgelist.with.distances,
#                                 normalized.densities, n, offset = FALSE) {
#   
#   final.edgelist.with.distances <- c()
#   for (i in names(normalized.densities)) {
#     which.ind <- c(which(edgelist.with.distances[, 1] == as.numeric(i)),
#                    which(edgelist.with.distances[, 2] == as.numeric(i)))
#     tmp.edgelist <- edgelist.with.distances[which.ind, ]
#     tmp.edgelist <- tmp.edgelist[order(tmp.edgelist[, "values"]), ]
#     tmp.edgelist <- tmp.edgelist[1:unname(normalized.densities[i]), ]
#     final.edgelist.with.distances <- rbind(final.edgelist.with.distances,
#                                            tmp.edgelist)
#   }
#   final.edgelist.with.distances <- removeDuplicateValues(final.edgelist.with.distances)
#   # note: weight here is really distance, will be converted to weight after the graph is completed
#   print("final.edgelist.with.distances")
#   print(final.edgelist.with.distances)
#   vertices.edges <- as.vector(t(as.matrix(final.edgelist.with.distances[, 1:2])))
#   output.graph <- add.edges(output.graph, edges = vertices.edges,
#                             weight = final.edgelist.with.distances[, 3],
#                             label = "EXTRA", sequence_assignment = n)
#   # return(output.graph)
#   return(list(output.graph = output.graph,
#               final.edgelist.with.distances = final.edgelist.with.distances))
# }



CheckMSTEdges <- function(output.graph, cluster.distances.matrix,
                          table.lengths, offset, n) {
  # create dummy graph contains n and n+1 vertices, but n_n section is fully connected so
  # it will act as a single connected component while n_n+1 and n+1_n+1 sections are not
  # connected, so each n+1 vertex will be its own disconnected component
  # print("cluster.distances.matrix")
  # print(cluster.distances.matrix)
  # cat("dim(cluster.distances.matrix)", dim(cluster.distances.matrix), "\n")
  # cat("colnames(cluster.distances.matrix)", colnames(cluster.distances.matrix), "\n")
  # cat("rownames(cluster.distances.matrix)", rownames(cluster.distances.matrix), "\n")
  adjacency.nn <- cluster.distances.matrix[1:table.lengths[1], 1:table.lengths[1]]
  el.nn.with.dist <- cbind(as.vector(row(adjacency.nn)),
                           as.vector(col(adjacency.nn)), 
                           as.vector(adjacency.nn))
  dummy.graph <- graph.empty(n = nrow(cluster.distances.matrix), directed = FALSE)
  dummy.graph <- add.edges(dummy.graph, edges = as.vector(t(el.nn.with.dist[, 1:2])))
  clu <- igraph::clusters(dummy.graph)
  # get the cluter that each node belongs to
  members <- clu$membership
  # initialize matrix to hold distances between graph components
  link.weights <- matrix(nrow = clu$no, ncol = clu$no)
  # loop through all components vs. all components to fill in distances
  for (i in 1:clu$no) {
    members.i <- which(members == i)
    for (j in 1:clu$no) {
      members.j <- which(members == j)
      # get the minimum distance between components i and j
      # (note members has zero index)
      tmp.matrix <- cluster.distances.matrix[members.i, members.j]
      link.weights[i, j] <- min(tmp.matrix)
    }
  }
  # set i-i distances to Inf instead of 0 so they aren't picked
  # as the closest neighbors
  for (i in 1:ncol(link.weights)) {
    link.weights[i, i] <- Inf
  }
  # make the minimum spanning tree for graph components
  components.adjacency <- graph.adjacency(link.weights, mode = "undirected",
                                          weighted = TRUE, diag = TRUE)
  components.mst <- minimum.spanning.tree(components.adjacency)
  components.mst.el <- get.edgelist(components.mst, names = FALSE)
  # go from MST edges to the actual edges between individual nodes
  # print("components.mst.el")
  # print(components.mst.el)
  for (i in 1:nrow(components.mst.el)) {
    # get index of shortest connection between components
    members.x <- which(members == components.mst.el[i, 1])
    members.y <- which(members == components.mst.el[i, 2])
    tmp.matrix <- as.matrix(cluster.distances.matrix[members.x, members.y])
    if(length(members.x) == 1){
      tmp.matrix <- t(tmp.matrix)
    }
    tmp.min <- min(tmp.matrix) # get the largest weight = shortest distance
    tmp.index <- which(tmp.matrix == tmp.min, arr.ind = TRUE)
    # print("tmp.index")
    # print(tmp.index)
    # cat("tmp.index[1]", tmp.index[1], "\n")
    # cat("tmp.index[2]", tmp.index[2], "\n")
    # add new edge to graph with the shortest connection vertices and the weight from linkweights
    # cat("members.x are", members.y, "\n")
    # cat("members.y are", members.y, "\n")
    v1 <- members.x[tmp.index[1]] + offset
    v2 <- members.y[tmp.index[2]] + offset
    # cat("v1 is", v1, "and v2 is", v2, "\n")
    if (are.connected(output.graph, v1, v2)) {
      E(output.graph, P = c(members.x[tmp.index[1]] + offset,
                            members.y[tmp.index[2]] + offset))$label <- "MST"
    } else {
      output.graph <- add.edges(output.graph, c(members.x[tmp.index[1]] + offset,
                                                members.y[tmp.index[2]] + offset),
                                weight = tmp.min, label = "MST",
                                sequence_assignment = n)
    }
  }
  return(output.graph)
}


InitializeGraph <- function(FLOWMAP.clusters) {
  # This section initializes the graph
  total.nodes <- length(FLOWMAP.clusters$full.clusters[, 1])
  # create empty graph with the right number of nodes - will fill in the edges later
  initial.graph <- graph.empty(n = total.nodes, directed = FALSE)
  return(initial.graph)
}

BuildFirstFLOWMAP <- function(FLOWMAP.clusters, per, min, max, distance.metric,
                              clustering.var) {
  n <- 0
  output.graph <- InitializeGraph(FLOWMAP.clusters = FLOWMAP.clusters)
  # This section creates a flowmap for the first time point
  cat("Building first FLOWMAP:\n")
  table.lengths <- FLOWMAP.clusters$table.lengths
  # cat("table.lengths are", table.lengths, "\n")
  # get distance matrix from clusters
  clusters <- FLOWMAP.clusters$full.clusters[1:table.lengths[1], ]
  clusters <- subset(clusters, select = clustering.var)
  # cat("colnames(clusters) are", colnames(clusters), "\n")
  # cat("dim(clusters) is", dim(clusters), "\n")
  cluster.distances <- dist(clusters, method = distance.metric, diag = TRUE, upper = TRUE)
  # cat("dim(cluster.distances) is", dim(cluster.distances), "\n")
  cluster.distances.matrix <- as.matrix(cluster.distances)
  cluster.distances.matrix <- as.data.frame(cluster.distances.matrix)
  # cat("rownames(cluster.distances.matrix) are", rownames(cluster.distances.matrix), "\n")
  # cat("colnames(cluster.distances.matrix) are", colnames(cluster.distances.matrix), "\n")
  
  # set i-i distances to Inf instead of 0 so they aren't the closest neighbors
  for (i in 1:ncol(cluster.distances.matrix)) {
    cluster.distances.matrix[i, i] <- Inf
  }
  numcluster <- nrow(clusters)
  # cat("numcluster is", numcluster, "\n")
  normalized.results <- FindNormalized(cluster.distances.matrix = cluster.distances.matrix,
                                       per = per, min = min, max = max, numcluster = numcluster,
                                       table.lengths = FALSE)
  normalized.densities <- normalized.results$normalized.densities
  # print("normalized.densities")
  # print(normalized.densities)
  edgelist.with.distances <- normalized.results$edgelist.with.distances
  # print("edgelist.with.distances")
  # print(edgelist.with.distances)
  # print("head(edgelist.with.distances)")
  # print(head(edgelist.with.distances))
  # print("dim(edgelist.with.distances)")
  # print(dim(edgelist.with.distances))
  
  # build new edgelist with N edges for each cluster based on normalized density
  results <- DrawNormalizedEdges(output.graph = output.graph,
                                 cluster.distances.matrix = cluster.distances.matrix,
                                 normalized.densities = normalized.densities,
                                 n = n, offset = FALSE)
  # results <- DrawNormalizedEdges(output.graph = output.graph,
  #                                     edgelist.with.distances = edgelist.with.distances,
  #                                     normalized.densities = normalized.densities,
  #                                     n = n, offset = FALSE)
  output.graph <- results$output.graph
  # output.graph <- DrawNormalizedEdges(output.graph = output.graph,
  #                                     edgelist.with.distances = edgelist.with.distances,
  #                                     normalized.densities = normalized.densities,
  #                                     n = n, offset = FALSE)
  # print("E(output.graph) post-normalized, pre-MST")
  # print(E(output.graph))
  edgelist.save <- results$final.edgelist.with.distances
  
  # now add all MST edges that are not yet included in the graph, and annotate all as "MST"
  adjacency.graph <- graph.adjacency(as.matrix(cluster.distances.matrix),
                                     mode = "undirected", weighted = "weight", diag = FALSE)
  mst.graph <- minimum.spanning.tree(adjacency.graph, algorithm = "prim")
  # for each edge of the mst, if it exists in the graph then label MST, if it doesn't exist then add it
  mst.graph.edgelist <- cbind(get.edgelist(mst.graph), E(mst.graph)$weight)
  class(mst.graph.edgelist) <- "numeric"
  # for each edge of the mst, if it exists in the graph then label MST, if it doesn't exist then add it
  for (i in 1:nrow(mst.graph.edgelist)) {
    if (are.connected(output.graph, mst.graph.edgelist[i, 1], mst.graph.edgelist[i, 2])) {
      E(output.graph, P = c(mst.graph.edgelist[i, 1],
                            mst.graph.edgelist[i, 2]))$label <- "MST"
    } else {
      output.graph <- add.edges(output.graph,
                                as.numeric(mst.graph.edgelist[i, 1:2]),
                                weight = mst.graph.edgelist[i, 3],
                                label = "MST", sequence_assignment = n)
    }
  }
  
  # print("E(output.graph) post-MST")
  # print(E(output.graph))
  
  # # now add all MST edges that are not yet included in the graph, and annotate all as "MST"
  # adjacency.graph <- graph.adjacency(as.matrix(cluster.distances.matrix),
  #                                    mode = "undirected", weighted = "weight", diag = FALSE)
  # mst.graph <- minimum.spanning.tree(adjacency.graph, algorithm = "prim")
  # # for each edge of the mst, if it exists in the graph then label MST, if it doesn't exist then add it
  # mst.graph.edgelist <- cbind(get.edgelist(mst.graph), E(mst.graph)$weight)
  # class(mst.graph.edgelist) <- "numeric"
  # # for each edge of the mst, if it exists in the graph then label MST, if it doesn't exist then add it
  # for (i in 1:nrow(mst.graph.edgelist)) {
  #   if (are.connected(output.graph, mst.graph.edgelist[i, 1], mst.graph.edgelist[i, 2])) {
  #     E(output.graph, P = c(mst.graph.edgelist[i, 1],
  #                           mst.graph.edgelist[i, 2]))$label <- "MST"
  #   } else {
  #     output.graph <- add.edges(output.graph,
  #                               as.numeric(mst.graph.edgelist[i, 1:2]),
  #                               weight = mst.graph.edgelist[i, 3],
  #                               label = "MST", sequence_assignment = 1)
  #   }
  # }
  
  
  # return(output.graph)
  return(list(output.graph = output.graph,
              edgelist.save = edgelist.save))
}


BuildFLOWMAP <- function(FLOWMAP.clusters, per, min, max,
                         distance.metric, cellnum, clustering.var) {
  edgelist.save <- list()
  first.results <- BuildFirstFLOWMAP(FLOWMAP.clusters = FLOWMAP.clusters,
                                    per = per, min = min, max = max,
                                    distance.metric = distance.metric,
                                    clustering.var = clustering.var)
  output.graph <- first.results$output.graph
  edgelist.save[["first"]] <- first.results$edgelist.save
  # output.graph <- BuildFirstFLOWMAP(FLOWMAP.clusters = FLOWMAP.clusters,
  #                                   per = per, min = min, max = max,
  #                                   distance.metric = distance.metric,
  #                                   clustering.var = clustering.var)
  table.breaks <- c(0, FLOWMAP.clusters$table.breaks)
  # cat("table.breaks are", table.breaks, "\n")
  # This section builds the flowmap one timepoint at a time
  for (n in 1:(length(FLOWMAP.clusters$cluster.medians) - 1)) {
    # offset value is used to correctly index the edges at each sequential step
    offset <- table.breaks[n]
    # cat("offset is", offset, "\n")
    # go through sequential cluster sets, add edges for each a-a and a-a+1 set, and also label sequential MST
    cat("Building FLOWMAP from", n, "to", n + 1, "\n")
    # get clusters for time a and a+1
    clusters <- rbind(FLOWMAP.clusters$cluster.medians[[n]], FLOWMAP.clusters$cluster.medians[[n + 1]])
    clusters <- subset(clusters, select = clustering.var)
    # cat("colnames(clusters) are", colnames(clusters), "\n")
    numcluster <- nrow(clusters)
    # cat("numcluster is", numcluster, "\n")
    # cat("dim(clusters) is", dim(clusters), "\n")
    # make adjacency matrix from clusters
    cluster.distances <- dist(clusters, method = distance.metric, diag = TRUE, upper = TRUE)
    # cat("dim(cluster.distances) is", dim(cluster.distances), "\n")
    cluster.distances.matrix <- as.matrix(cluster.distances)
    cluster.distances.matrix <- as.data.frame(cluster.distances.matrix)
    rownames(cluster.distances.matrix) <- (offset + 1):table.breaks[n + 2]
    colnames(cluster.distances.matrix) <- (offset + 1):table.breaks[n + 2]
    # cat("rownames(cluster.distances.matrix) are", rownames(cluster.distances.matrix), "\n")
    # cat("colnames(cluster.distances.matrix) are", colnames(cluster.distances.matrix), "\n")
    # cluster.distances.matrix <- as.matrix(cluster.distances)
    # set i-i distances to Inf instead of 0 so they aren't the closest neighbors
    for (i in 1:ncol(cluster.distances.matrix)) {
      # cat("i is", i, "\n")
      cluster.distances.matrix[i, i] <- Inf
    }
    # This section adds the lowest distance n_n+1 and n+1_n+1 edges to the output graph
    n_n1.table.lengths <- FLOWMAP.clusters$table.lengths[n:(n + 1)]
    n_n1.table.breaks <- table.breaks[n:(n + 1)]
    # cat("n_n1.table.lengths are", n_n1.table.lengths, "\n")
    # cat("n_n1.table.breaks are", n_n1.table.breaks, "\n")
    normalized.results <- FindNormalized(cluster.distances.matrix = cluster.distances.matrix,
                                         per = per, min = min, max = max, numcluster = numcluster,
                                         table.lengths = n_n1.table.lengths,
                                         table.breaks = n_n1.table.breaks,
                                         offset = offset)
    normalized.densities <- normalized.results$normalized.densities
    # print("normalized.densities")
    # print(normalized.densities)
    edgelist.with.distances <- normalized.results$edgelist.with.distances
    # print("edgelist.with.distances")
    # print(edgelist.with.distances)
    # print("head(edgelist.with.distances)")
    # print(head(edgelist.with.distances))
    # print("dim(edgelist.with.distances)")
    # print(dim(edgelist.with.distances))
    
    # build new edgelist with N edges for each cluster based on normalized density
    results <- DrawNormalizedEdges(output.graph = output.graph,
                                   cluster.distances.matrix = cluster.distances.matrix,
                                   normalized.densities = normalized.densities,
                                   n = n, offset = offset)
    # results <- DrawNormalizedEdges(output.graph = output.graph,
    #                                     edgelist.with.distances = edgelist.with.distances,
    #                                     normalized.densities = normalized.densities,
    #                                     n = n, offset = offset)
    output.graph <- results$output.graph
    edgelist.save[[n]] <- results$final.edgelist.with.distances
    # output.graph <- DrawNormalizedEdges(output.graph = output.graph,
    #                                     edgelist.with.distances = edgelist.with.distances,
    #                                     normalized.densities = normalized.densities,
    #                                     n = n, offset = offset)
    # print("E(output.graph) post-normalized, pre-MST")
    # print(E(output.graph))

    # This section adds the "MST" for n_n+1 and n+1_n+1 nodes
    output.graph <- CheckMSTEdges(output.graph = output.graph,
                                  cluster.distances.matrix = cluster.distances.matrix, 
                                  table.lengths = n_n1.table.lengths, n = n, offset = offset)
    # print("E(output.graph) post-MST")
    # print(E(output.graph))
  }
  
  # convert graph distances to weights (low distance = high weight and vice versa)
  distances <- E(output.graph)$weight
  # print("E(output.graph)$weight")
  # print(E(output.graph)$weight)
  # weights <- -distances + max(distances) + min(distances)
  weights <- (1 / distances)
  E(output.graph)$weight <- weights
  # print("E(output.graph)$weight")
  # print(E(output.graph)$weight)
  # print("final E(output.graph)")
  # print(E(output.graph))
  output.graph <- AnnotateGraph(output.graph = output.graph,
                                FLOWMAP.clusters = FLOWMAP.clusters,
                                cellnum = cellnum)
  return(list(output.graph = output.graph,
              edgelist.save = edgelist.save))
}


AnnotateGraph <- function(output.graph, FLOWMAP.clusters, cellnum) {
  # This section annotates the graph
  anno.cat <- c()
  anno <- list()
  # print("FLOWMAP.clusters$full.clusters")
  # print(FLOWMAP.clusters$full.clusters)
  # iterate through all times and annotate
  for (f in 1:length(FLOWMAP.clusters$cluster.medians)) {
    cat("Annotating graph for file", f, "\n")
    # get medians for all parameters and counts for all clusters
    counts <- FLOWMAP.clusters$cluster.counts[[f]]$Counts
    anno$count <- counts
    anno$percent.total <- data.frame(percent.total = c(counts / cellnum))
    anno$medians  <- FLOWMAP.clusters$cluster.medians[[f]]
    # add time information column
    time.matrix <- matrix(f, nrow = length(anno$count))
    colnames(time.matrix) <- c("Timepoint")
    anno$medians <- cbind(anno$medians, time.matrix)
    # add median and percent values
    # print("anno$count")
    # print(anno$count)
    # print("anno$medians")
    # print(anno$medians)
    for (a in c("medians", "percent.total")) {
      # anno.cat is all the annos concatenated, will be
      # used to make "anno.Rsave" file
      anno.cat[[a]] <- rbind(anno.cat[[a]], anno[[a]])
    }
    # print("anno.cat")
    # print(anno.cat)
  }
  # combine anno.cat matrices
  output.anno <- cbind(anno.cat[[1]], anno.cat[[2]])
  # print("output.anno")
  # print(output.anno)
  for (c in colnames(output.anno)) {
    # cat("c is", c, "\n")
    output.graph <- set.vertex.attribute(output.graph, c,
                                         index = as.numeric(rownames(output.anno)),
                                         value = output.anno[, c])
    # print("as.numeric(rownames(output.anno))")
    # print(as.numeric(rownames(output.anno)))
    # print("output.anno[, c]")
    # print(output.anno[, c])
  }
  # add name attribute
  # print("V(output.graph)$name")
  # print(V(output.graph)$name)
  V(output.graph)$name <- 1:length(V(output.graph))
  # print("V(output.graph)$name")
  # print(V(output.graph)$name)
  return(output.graph)
}

