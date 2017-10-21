
InitializeMultiGraph <- function(list.of.FLOWMAP.clusters) {
  # This section initializes the graph
  total.nodes <- dim(list.of.FLOWMAP.clusters$full.clusters)[1]
  initial.graph <- graph.empty(n = total.nodes, directed = FALSE)
  return(initial.graph)
}

BuildFirstMultiFLOWMAP <- function(list.of.FLOWMAP.clusters, per, min, max, distance.metric,
                                   clustering.var) {
  output.graph <- InitializeMultiGraph(list.of.FLOWMAP.clusters)
  cat("Building first FLOWMAP\n")
  # This section creates a flowmap for the first time point
  clusters <- list.of.FLOWMAP.clusters$cluster.medians[[1]]
  # get distance matrix from clusters
  clusters <- subset(clusters, select = clustering.var)
  cluster.distances <- dist(clusters, method = distance.metric, diag = TRUE, upper = TRUE)
  cluster.distances.matrix <- as.matrix(cluster.distances)
  # set i-i distances to Inf instead of 0 so they aren't the closest neighbors
  for (i in 1:ncol(cluster.distances.matrix)) {
    cluster.distances.matrix[i, i] <- Inf
  }
  numcluster <- nrow(clusters)
  results <- FindNormalized(cluster.distances.matrix,
                            per, min, max, numcluster,
                            table.lengths = FALSE)
  normalized.densities <- results$normalized.densities
  # build new edgelist with N edges for each cluster based on normalized density
  cat("Building first edgelist\n")
  results <- DrawNormalizedEdges(output.graph, cluster.distances.matrix,
                                 normalized.densities, n = 1)
  output.graph <- results$output.graph
  # now add all MST edges that are not yet included in the graph, and annotate all as "MST"
  adjacency.graph <- graph.adjacency(cluster.distances.matrix,
                                     mode = "undirected", weighted = "weight", diag = FALSE)
  mst.graph <- minimum.spanning.tree(adjacency.graph, algorithm = "prim")
  # for each edge of the mst, if it exists in the graph then label MST, if it doesn't exist then add it
  mst.graph.edgelist <- cbind(get.edgelist(mst.graph), E(mst.graph)$weight)
  class(mst.graph.edgelist) <- "numeric"
  # for each edge of the mst, if it exists in the graph then label MST, if it doesn't exist then add it
  for (i in 1:nrow(mst.graph.edgelist)) {
    if (are.connected(output.graph,
                      mst.graph.edgelist[i, 1],
                      mst.graph.edgelist[i, 2])) {
      E(output.graph, P = c(mst.graph.edgelist[i, 1],
                            mst.graph.edgelist[i, 2]))$label <- "MST"
    } else {
      output.graph <- add.edges(output.graph,
                                as.numeric(mst.graph.edgelist[i, 1:2]),
                                weight = mst.graph.edgelist[i, 3],
                                label = "MST", sequence.assignment = 1)
    }
  }
  return(output.graph)
}

BuildMultiFLOWMAP <- function(list.of.FLOWMAP.clusters, per, min,
                              max, distance.metric, label.key, clustering.var) {
  remodel.FLOWMAP.clusters <- RemodelFLOWMAPClusterList(list.of.FLOWMAP.clusters)
  output.graph <- BuildFirstMultiFLOWMAP(remodel.FLOWMAP.clusters, per,
                                         min, max, distance.metric = distance.metric,
                                         clustering.var)
  # put each conditions clusters together into one timepoint
  table.breaks <- c(0, remodel.FLOWMAP.clusters$table.breaks)
  # This section builds the flowmap one timepoint at a time
  for (n in 1:(length(remodel.FLOWMAP.clusters$cluster.medians) - 1)) {
    # offset value is used to correctly index the edges at each sequential step
    offset <- table.breaks[n]
    # go through sequential cluster sets, add edges for each a-a and a-a+1 set, and also label sequential MST
    cat("Build FLOWMAP from", n, "to", n + 1, "\n")
    # get clusters for time a and a+1
    clusters <- rbind(remodel.FLOWMAP.clusters$cluster.medians[[n]], remodel.FLOWMAP.clusters$cluster.medians[[n + 1]])
    clusters <- subset(clusters, select = clustering.var)
    numcluster <- nrow(clusters)
    # make adjacency matrix from clusters
    cluster.distances <- dist(clusters, method = distance.metric, diag = TRUE, upper = TRUE)
    cluster.distances.matrix <- as.matrix(cluster.distances)
    rownames(cluster.distances.matrix) <- (offset + 1):table.breaks[n + 2]
    colnames(cluster.distances.matrix) <- (offset + 1):table.breaks[n + 2]
    # set i-i distances to Inf instead of 0 so they aren't the closest neighbors
    for (i in 1:ncol(cluster.distances.matrix)) {
      cluster.distances.matrix[i, i] <- Inf
    }
    # This section adds the lowest distance n_n+1 and n+1_n+1 edges to the output graph
    n_n1.table.lengths <- remodel.FLOWMAP.clusters$table.lengths[n:(n + 1)]
    n_n1.table.breaks <- remodel.FLOWMAP.clusters$table.breaks[n:(n + 1)]
    results <- FindNormalized(cluster.distances.matrix, per, min,
                              max, numcluster, table.lengths = n_n1.table.lengths,
                              table.breaks = n_n1.table.breaks,
                              offset = offset)
    normalized.densities <- results$normalized.densities
    # build new edgelist with N edges for each cluster based on normalized density
    cat("Building edgelist\n")
    results <- DrawNormalizedEdges(output.graph, cluster.distances.matrix,
                                   normalized.densities, offset = offset, n = n)
    output.graph <- results$output.graph
    # This section adds the "MST" for n_n+1 and n+1_n+1 nodes
    output.graph <- CheckMSTEdges(output.graph, cluster.distances.matrix, 
                                  n_n1.table.lengths, offset = offset, n = n)
  }
  # convert graph distances to weights (low distance = high weight and vice versa)
  distances <- E(output.graph)$weight
  #### TEMPORARY FIX FOR IDENTICAL CELLS WITH DIST = 0, WEIGHT = INF
  fix.identical.dist <- which(distances == 0)
  if (length(fix.identical.dist) != 0) {
    distances.no.identical <- distances[-fix.identical.dist]
    distances[fix.identical.dist] <- min(distances.no.identical)
  }
  ####
  weights <- 1 / distances
  E(output.graph)$weight <- weights
  output.graph <- AnnotateMultiGraph(output.graph, remodel.FLOWMAP.clusters, label.key)
  return(output.graph)
}

AnnotateMultiGraph <- function(output.graph, list.of.FLOWMAP.clusters,
                               label.key) {
  # This section annotates the graph
  anno <- list()
  # iterate through all times and annotate
  times <- 1:length(list.of.FLOWMAP.clusters$cluster.medians)
  anno$medians <- data.frame()
  anno$count <- data.frame()
  anno$percent.total <- data.frame()
  for (t in times) {
    cat("Annotating graph for file", t, "\n")
    # get medians for all parameters and counts for all clusters
    anno$medians <- rbind(anno$medians, list.of.FLOWMAP.clusters$cluster.medians[[t]])
    anno$medians[, "Condition"]
    anno$count <- rbind(anno$count, list.of.FLOWMAP.clusters$cluster.counts[[t]])
    total.cell <- sum(list.of.FLOWMAP.clusters$cluster.counts[[t]])
    percent.total <- data.frame(list.of.FLOWMAP.clusters$cluster.counts[[t]] / total.cell)
    colnames(percent.total) <- "percent.total"
    anno$percent.total <- rbind(anno$percent.total,
                                percent.total)
  }
  output.anno <- cbind(anno$medians, anno$count, anno$percent.total)
  print("label.key")
  print(label.key)
  output.anno <- ConvertCharacterLabel(output.anno, label.key)
  for (c in colnames(output.anno)) {
    output.graph <- set.vertex.attribute(output.graph, c,
                                         index = as.numeric(1:dim(output.anno)[1]),
                                         value = output.anno[, c])
  }
  # add name attribute
  V(output.graph)$name <- 1:length(V(output.graph))
  return(output.graph)
}
