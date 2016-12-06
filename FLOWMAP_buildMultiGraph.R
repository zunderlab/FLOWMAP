
RemodelFLOWMAPClusterList <- function(list.of.FLOWMAP.clusters) {
  # take FLOWMAP of conditions with timeseries and make into one timeseries
  # combine FLOWMAP conditions
  conditions <- names(list.of.FLOWMAP.clusters)
  timepoints <- length(list.of.FLOWMAP.clusters[[conditions[1]]]$cluster.medians)
  full.clusters <- data.frame()
  table.breaks <- c()
  table.lengths <- c()
  cluster.medians <- list()
  cluster.counts <- list()
  cell.assgn <- list()
  for (t in 1:timepoints) {
    temp.medians <- data.frame()
    temp.cell.assgn <- data.frame()
    temp.counts <- data.frame()
    for (condition in conditions) {
      temp.medians <- rbind(temp.medians, list.of.FLOWMAP.clusters[[condition]]$cluster.medians[[t]])
      temp.cell.assgn <- rbind(temp.cell.assgn, list.of.FLOWMAP.clusters[[condition]]$cell.assgn[[t]])
      temp.counts <- rbind(temp.counts, list.of.FLOWMAP.clusters[[condition]]$cluster.counts[[t]])
    }
    cluster.medians[[t]] <- temp.medians
    cluster.counts[[t]] <- temp.counts
    cell.assgn[[t]] <- temp.cell.assgn
    table.lengths <- c(table.lengths, dim(temp.medians)[1])
    table.breaks <- c(table.breaks, sum(table.lengths))
    full.clusters <- rbind(full.clusters, temp.medians)
  }
  remodeled.FLOWMAP.clusters <- FLOWMAPcluster(full.clusters, table.breaks, table.lengths,
                                               cluster.medians, cluster.counts, cell.assgn) 
  return(remodeled.FLOWMAP.clusters)
}


InitializeMultiGraph <- function(list.of.FLOWMAP.clusters) {
  # This section initializes the graph
  total.nodes <- 0
  for (condition in names(list.of.FLOWMAP.clusters)) {
    FLOWMAP.clusters <- list.of.FLOWMAP.clusters[[condition]]
    total.nodes <- total.nodes + length(FLOWMAP.clusters$full.clusters[, 1])
    # create empty graph with the right number of nodes - will fill in the edges later
  }
  initial.graph <- graph.empty(n = total.nodes, directed = FALSE)
  return(initial.graph)
}


BuildFirstMultiFLOWMAP <- function(list.of.FLOWMAP.clusters, per, min, max, distance.metric) {
  n <- 1
  output.graph <- InitializeMultiGraph(list.of.FLOWMAP.clusters)
  # This section creates a flowmap for the first time point
  cat("Building first FLOWMAP\n")
  clusters <- c()
  for (condition in names(list.of.FLOWMAP.clusters)) {
    table.length <- list.of.FLOWMAP.clusters[[condition]]$table.lengths[1]
    current.clusters <- list.of.FLOWMAP.clusters[[condition]]
    clusters <- rbind(clusters, current.clusters$full.clusters[1:table.length, ])
  }
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
                                 normalized.densities, n = n)
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
      E(output.graph,P = c(mst.graph.edgelist[i, 1],
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


BuildMultiFLOWMAP <- function(list.of.FLOWMAP.clusters, conditions,
                              per, min, max, distance.metric, cellnum) {
  output.graph <- BuildFirstMultiFLOWMAP(list.of.FLOWMAP.clusters, per,
                                         min, max, distance.metric = distance.metric)
  remodel.FLOWMAP.clusters <- RemodelFLOWMAPClusterList(list.of.FLOWMAP.clusters)
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
  # print("sum(distances == 0)")
  # print(sum(distances == 0))
  
  #### TEMPORARY FIX FOR IDENTICAL CELLS WITH DIST = 0, WEIGHT = INF
  fix.identical.dist <- which(distances == 0)
  distances.no.identical <- distances[-fix.identical.dist]
  # print("min(distances.no.identical)")
  # print(min(distances.no.identical))
  distances[fix.identical.dist] <- min(distances.no.identical)
  # print("distances")
  # print(distances)
  ####
  
  weights <- 1 / distances
  E(output.graph)$weight <- weights
  
  output.graph <- AnnotateMultiGraph(output.graph, list.of.FLOWMAP.clusters, cellnum)
  return(output.graph)
}


AnnotateMultiGraph <- function(output.graph, list.of.FLOWMAP.clusters, cellnum) {
  # This section annotates the graph
  anno.cat <- c()
  anno <- list()
  # iterate through all times and annotate
  x <- names(list.of.FLOWMAP.clusters)[1]
  for (f in 1:length(list.of.FLOWMAP.clusters[[x]]$cluster.medians)) {
    for (condition in names(list.of.FLOWMAP.clusters)) {
      cat("Annotating graph for file", f, "\n")
      # get medians for all parameters and counts for all clusters
      counts <- list.of.FLOWMAP.clusters[[condition]]$cluster.counts[[f]]$Counts
      anno$count <- counts
      anno$percent.total <- data.frame(percent.total = c(counts / cellnum))
      anno$medians  <- list.of.FLOWMAP.clusters[[condition]]$cluster.medians[[f]]
      # add time information column
      time.matrix <- matrix(f, nrow = length(anno$count))
      colnames(time.matrix) <- c("Timepoint")
      # add condition information column
      condition.matrix <- matrix(condition, nrow = length(anno$count))
      colnames(condition.matrix) <- c("Condition")
      anno$medians <- cbind(anno$medians, time.matrix, condition.matrix)
      # add median and percent values
      for (a in c("medians", "percent.total")) {
        # anno_cat is all the annos concatenated, will be
        # used to make "anno.Rsave" file
        anno.cat[[a]] <- rbind(anno.cat[[a]], anno[[a]])
      }
    }
  }
  # combine anno_cat matrices
  output.anno <- cbind(anno.cat[[1]], anno.cat[[2]])
  for (c in colnames(output.anno)) {
    output.graph <- set.vertex.attribute(output.graph, c,
                                         index = as.numeric(rownames(output.anno)),
                                         value = output.anno[, c])
  }
  # add name attribute
  V(output.graph)$name <- 1:length(V(output.graph))
  return(output.graph)
}



