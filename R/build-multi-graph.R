
InitializeMultiGraph <- function(list.of.FLOWMAP.clusters) {
  # This section initializes the graph
  total.nodes <- dim(list.of.FLOWMAP.clusters$full.clusters)[1]
  initial.graph <- graph.empty(n = total.nodes, directed = FALSE)
  return(initial.graph)
}

BuildFirstMultiFLOWMAPkNN <- function(list.of.FLOWMAP.clusters, k, min, max, distance.metric,
                                   clustering.var) {
  output.graph <- InitializeMultiGraph(list.of.FLOWMAP.clusters)
  cat("Building first FLOWMAP\n")
  
  ## Read in info from FLOWMAP.clusters
  table.lengths <- list.of.FLOWMAP.clusters$table.lengths
  table.breaks <- c(0, list.of.FLOWMAP.clusters$table.breaks)

  clusters <- list.of.FLOWMAP.clusters$cluster.medians[[1]]
  clusters <- subset(clusters, select = clustering.var)
  
  output.graph <- graph.empty()
  
  results <- BaseBuildKNN(clusters=clusters,table.breaks=table.breaks, offset=0, n=0,
                          output.graph=output.graph, k=k, min=min, max=max,
                          distance.metric=distance.metric)
  return(results)
}

BuildMultiFLOWMAPkNN <- function(remodel.FLOWMAP.clusters, k, min,
                              max, distance.metric, label.key, clustering.var) {
 
  global.remodel.FLOWMAP.clusters <<- remodel.FLOWMAP.clusters
  
  first.results <- BuildFirstMultiFLOWMAPkNN(remodel.FLOWMAP.clusters, k,
                                         min, max, distance.metric = distance.metric,
                                         clustering.var)
  output.graph <- first.results$output.graph
  knn.indexes <- list()
  knn.distances <- list()
  knn.indexes[[1]] <- first.results$indexes  
  knn.distances[[1]] <- first.results$distances
  
  ## Read in info for indexing timepoint clusters from FLOWMAP.clusters
  ### (put each conditions clusters together into one timepoint)
  table.lengths <- remodel.FLOWMAP.clusters$table.lengths
  table.breaks <- c(0, remodel.FLOWMAP.clusters$table.breaks)
  
  ## This section builds the flowmap one timepoint at a time for subsequent timepoints
  
  for (n in 1:(length(remodel.FLOWMAP.clusters$cluster.medians) - 1)) {
    ## Offset value is used to correctly index the edges at each sequential step
    offset <- table.breaks[n] + 1
    cat(paste0("Starting next round, offset = ", as.character(offset), "\n"))
    
    ## Go through sequential cluster sets, add edges for each a-a and a-a+1 set
    cat(paste0("Building FLOWMAP from", n, "to", n + 1, "\n"))
    ## Get clusters for time a and a+1
    clusters <- rbind(remodel.FLOWMAP.clusters$cluster.medians[[n]], remodel.FLOWMAP.clusters$cluster.medians[[n + 1]])
    clusters <- subset(clusters, select = clustering.var)
    
    results <- BaseBuildKNN(clusters=clusters, table.breaks=table.breaks, offset=offset,
                            output.graph=output.graph,n=n, k=k, min=min, max=max,
                            distance.metric=distance.metric)
    
    output.graph <- results$output.graph
    #First fill in new connections from cells in previous timepoint
    knn.indexes[[n]] <- data.frame(cbind(knn.indexes[[n]],results$indexes[1:table.lengths[n],])) #to make cols for adding connections back from second timepoint
    knn.distances[[n]] <- data.frame(cbind(knn.distances[[n]],results$distances[1:table.lengths[n],]))
    #Then add for cells in new timepoint
    #knn.indexes[[n+1]] <- results$indexes[1:(nrow(results$indexes)/2),] 
    knn.indexes[[n+1]] <- results$indexes[(table.lengths[n]+1):nrow(results$indexes),]
    knn.distances[[n+1]] <- results$distances[(table.lengths[n]+1):nrow(results$distances),]
  }
  # PREP KNN FOR OUTPUT #######################################
  #make final timepoint data frame as well
  knn.indexes[[n+1]] <- data.frame(knn.indexes[[n+1]]) #to make cols for adding connections back from second timepoint
  knn.distances[[n+1]] <- data.frame(knn.distances[[n+1]])
  ## Remove duplicate connections====
  global.pre.simplify <<- knn.indexes
  keep_ind = lapply(knn.indexes, function(y) apply(y,1,function(x){
    to.keep <- c()
    for (i in unique(x)) {
      to.keep <- rbind(to.keep, head(which(x==i),1))
    }
    to.keep
  }))
  global.keep.ind <<- keep_ind
  #knn.ind.sim <- c()
  #knn.dist.sim <- c()
  
  knn.ind.sim <- lapply(1:(length(keep_ind)-1), function(x) {
    knn.indexes.keep <- c()
    for (i in 1:length(keep_ind[[x]])){
      knn.indexes.keep <- dplyr::bind_rows(knn.indexes.keep,knn.indexes[[x]][i,keep_ind[[x]][[i]]])
    }
    print("kept")
    knn.indexes.keep
  })
  knn.ind.sim <- dplyr::bind_rows(knn.ind.sim,knn.indexes[[length(knn.indexes)]]) #join into one data frame
  
  knn.dist.sim <- lapply(1:(length(keep_ind)-1), function(x) {
    knn.distances.keep <- c()
    for (i in 1:length(keep_ind[[x]])){
      knn.distances.keep <- dplyr::bind_rows(knn.distances.keep,knn.distances[[x]][i,keep_ind[[x]][[i]]])
    }
    knn.distances.keep
  })
  knn.dist.sim <- dplyr::bind_rows(knn.dist.sim,knn.distances[[length(knn.distances)]]) #join into one data frame
  print("post-simplify")
  
  knn.out <- list("indexes" = knn.ind.sim, "distances" = knn.dist.sim)
  global.knn.out.pre <<- knn.out
  # knn.out$indexes <- tidyr::unnest(knn.out$indexes, cols = colnames(knn.out$indexes))
  # knn.out$distances <- tidyr::unnest(knn.out$distances, cols = colnames(knn.out$distances))
  # for(i in 1:nrow(knn.out$indexes)) knn.out$indexes[i,] <- knn.out$indexes[i,order(unlist(knn.out$distances[i,]))]
  # for(i in 1:nrow(knn.out$distances)) knn.out$distances[i,] <- knn.out$distances[i,][order(knn.out$distances[i,])]
  # for(i in 1:ncol(knn.out$indexes)) knn.out$indexes[,i] <- unlist(knn.out$indexes[,i])
  # for(i in 1:ncol(knn.out$distances)) knn.out$distances[,i] <- unlist(knn.out$distances[,i])
  print("processing knn output")
  for(i in 1:ncol(knn.out$indexes)) knn.out$indexes[sapply(knn.out$indexes[,i], is.null),i] <- 0
  for(i in 1:ncol(knn.out$distances)) knn.out$distances[sapply(knn.out$distances[,i], is.null),i] <- 10000
  knn.out$indexes  <-  as.data.frame(matrix(unlist(knn.out$indexes), nrow=length(unlist(knn.out$indexes[1])))) 
  knn.out$distances  <-  as.data.frame(matrix(unlist(knn.out$distances), nrow=length(unlist(knn.out$distances[1]))))
  dist_order = t(apply(knn.out$distances, 1, order))
  for(i in 1:nrow(knn.out$indexes)) knn.out$indexes[i,] <- knn.out$indexes[i,dist_order[i,]] #order indexes by increasing distance (necessary to incorporate second timepoint connections into overall order)
  for(i in 1:nrow(knn.out$distances)) knn.out$distances[i,] <- knn.out$distances[i,][dist_order[i,]] #increasing order distances
  knn.out$indexes  <-  as.data.frame(matrix(unlist(knn.out$indexes), nrow=length(unlist(knn.out$indexes[1])))) 
  knn.out$distances  <-  as.data.frame(matrix(unlist(knn.out$distances), nrow=length(unlist(knn.out$distances[1]))))
  print("finished processing knn output")
  global.knn.out <<- knn.out
  
  # PREP GRAPH FOR OUTPUT #######################################
  global.pre.output.graph <<- output.graph
  output.graph.final <- graph.empty()
  vertices.edges <- get.edgelist(output.graph)
  vertices.edges <- as.matrix(data.frame(lapply(data.frame(vertices.edges), function(x) as.numeric(as.character(x)))))
  output.graph.final <- graph_from_edgelist(vertices.edges, directed = FALSE)
  E(output.graph.final)$weight <- unlist(lapply(E(output.graph)$weight, function(x) 1/as.numeric(x)))
  output.graph.final <- AnnotateMultiGraph(output.graph, remodel.FLOWMAP.clusters, label.key)
  ## Remove duplicates and self-edges
  output.graph.final <- simplify(output.graph.final)
  
  #TODO remove knn outputs, no longer needed for umap
  return(list("output.graph" = output.graph.final,
              "knn.out" = knn.out))
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
  global.output.anno.pre <<- output.anno
  print("label.key")
  print(label.key)
  output.anno <- ConvertCharacterLabel(output.anno, label.key)
  global.output.anno.post <<- output.anno
  for (c in colnames(output.anno)) {
    #if (c == "Condition") {
    #  output.graph <- set.vertex.attribute(output.graph, c,
    #                                       index = as.numeric(1:dim(output.anno)[1]),
    #                                       value = output.anno[, c]) #as.character()
    #} else {
    output.graph <- set.vertex.attribute(output.graph, c,
                                           index = as.numeric(1:dim(output.anno)[1]),
                                           value = output.anno[, c])
    #}
  }
  # add name attribute
  V(output.graph)$name <- 1:length(V(output.graph))
  global.output.graph <<- output.graph
  return(output.graph)
}
