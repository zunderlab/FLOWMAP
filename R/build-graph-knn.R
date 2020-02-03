# #############################################################################
# ##FUNCTIONS FOR KNN-BASED DENSITY GRAPH BUILDING====
# #############################################################################

# kNN density function ====
KnnDensity <- function(k, min, max, n, nn.ids.df, nn.dists.df,
                       numcluster = numcluster, #table.lengths
                       table.breaks, offset) {
  all.densities.list <- list()
  for(i in 1:nrow(nn.ids.df)) {
    ## Alternatives for calculating density
    ## Dist to kth nearest neighbor:
    density.by.knn <- abs(nn.dists.df[i,k])
    ## Transformed metric from X.shift paper:
    #density.by.Xshift <- (1/(n*(longest.dist^d))) * (sum(c(1:k)^d)/sum(compile.dists))^d
    ## Sum of distances from 1:kth nearest neighbor:
    #density.by.knn <- sum(abs(nn.dists.df[i,1:k]))
    ## Mean of distances from 1:kth nearest neighbor:
    #density.by.knn <- sum(abs(nn.dists.df[i,1:k]))/k
    all.densities.list[[paste(i,".nn", sep='')]] <- density.by.knn
  }
  densities.df <- rbind(matrix(unlist(all.densities.list), byrow=T))
  normalized.densities <- round(densities.df / max(densities.df) * (max - min) + min)
  if (offset == 0) {
    names(normalized.densities) <- (1:numcluster)
  } else {
    names(normalized.densities) <- offset:table.breaks[n + 2]
  }
  return(normalized.densities)
}

#Check connectedness for first flowmap ====
FirstConnectSubgraphs <- function(output.graph, edge.list, offset,
                                  table.breaks, n, distance.metric, clusters) {
  output.graph <- set.vertex.attribute(output.graph,'name',index=V(output.graph),as.character(1:vcount(output.graph)))
  ## Identify subgraphs of density-based nearest neighbor graph
  subgraphs.ls <- decompose.graph(output.graph)
  ## Change row names of cluster medians df to have correct indexing
  rownames(clusters) <- c(offset:table.breaks[n + 2])
  ## Make edge.list a df and change colnames
  edge.list <- as.data.frame(edgelist.save)
  colnames(edge.list) <- c("Vertex1", "Vertex2", "Distance")
  ## Set lists and ind for while loop
  y <- 0
  graphs.ls <- list()
  subgraphs.ls.el <- list()
  while (length(subgraphs.ls) >= 2) {
    y <- y + 1
    print(y)
    ## Compute inter-subgraph edges
    if (length(get.edgelist(subgraphs.ls[[2]])) != 0) {
      print("in loop")
      i <- 1
      j <- 1
      to.add.df <- data.frame()
      for (p in 1:length(subgraphs.ls)) {
        si.graph <- subgraphs.ls[[p]]
        ## Prepare df to hold 1 nn from this subraph to other subgraphs
        inter.sub.nns.df <- data.frame()
        for (m in 1:length(subgraphs.ls)){
          # if (m > length(subgraphs.ls)){
          #   break
          # }
          if (m == p){
            next
          }
          sj.graph <- subgraphs.ls[[m]]
          if (distance.metric == 'manhattan') {
            print("Manhattan distance no longer supported. Using euclinean distance.")
            inter.sub.nns <- RANN::nn2(data = clusters[V(sj.graph)$name,],
                                       query=clusters[V(si.graph)$name,],
                                       k=1, searchtype="priority", eps=0.1)
          } else if (distance.metric == 'euclidean') {
            inter.sub.nns <- RANN::nn2(data = clusters[V(sj.graph)$name,],
                                       query=clusters[V(si.graph)$name,],
                                       k=1, searchtype="priority", eps=0.1)
          }
          ## Add nn output to temp df, then bind to full list
          temp.inter.sub.nns.df <- data.frame("Vertex1" = V(si.graph)$name)
          temp.inter.sub.nns.df <- cbind(temp.inter.sub.nns.df,
                                         as.data.frame(V(sj.graph)$name[as.numeric(inter.sub.nns$nn.idx)]),
                                         as.data.frame(inter.sub.nns$nn.dists))
          colnames(temp.inter.sub.nns.df) <- c("Vertex1", "Vertex2", "Distance")
          inter.sub.nns.df <-  rbind(inter.sub.nns.df, temp.inter.sub.nns.df)
        }
        ## Put list of nearest neighbors in order with shortest dist at top
        inter.sub.nns.df <- inter.sub.nns.df[order(inter.sub.nns.df$Distance), ]
        to.add.df <- rbind(to.add.df, inter.sub.nns.df[1,])
      }
      output.graph.update <- graph.empty()
      og.el <- get.edgelist(output.graph)
      to.add.mat <- as.matrix(data.frame(lapply(data.frame(to.add.df),
                                                function(x) as.character(x))))
      vertices.edges <- rbind(og.el[,1:2], to.add.mat[,1:2])
      ogw.df <- data.frame()
      output.graph.update <- graph_from_edgelist(vertices.edges, directed = FALSE)
      ogw.df <- data.frame(E(output.graph)$weight)
      ta.df <- data.frame(to.add.df[,3])
      colnames(ogw.df)[1] <- "weight"
      colnames(ta.df)[1] <- "weight"
      ogw.ta.df <- data.frame(rbind(ogw.df, ta.df))
      E(output.graph.update)$weight <- ogw.ta.df[,1]
      output.graph.update <- set.vertex.attribute(output.graph.update,'name',
                                                  index=V(output.graph.update),
                                                  as.character(offset:table.breaks[n + 2]))
      subgraphs.ls <- decompose.graph(output.graph.update)
      subgraphs.el.ls <- list()
      for (ind in 1:length(subgraphs.ls)) {
        subgraphs.el.ls[[ind]] <- get.edgelist(subgraphs.ls[[ind]])
      }
      graphs.ls[[y]] <- output.graph.update
      subgraphs.ls.el[[y]] <- subgraphs.el.ls
      output.graph <- output.graph.update
    }#if
  }#while

  ## Return edge.list and output.graph as results
  colnames(edge.list) <- c("row.inds", "col.inds", "values")
  results <- list('output.graph' = output.graph, 'edgelist.with.distances' = edge.list)
  return(results)

}

# kNN verion of "CheckMSTEdges"
ConnectSubgraphs <- function(output.graph, edge.list, offset,
                             table.breaks, n, distance.metric, clusters) {
  print("in connected")
  ## Make graph from edge.list
  time.prox.graph <- graph.empty()
  time.prox.graph <- graph_from_edgelist(edge.list[,1:2], directed = FALSE)
  E(time.prox.graph)$weight <- edge.list[,3]
  time.prox.graph <- set.vertex.attribute(time.prox.graph,'name',index=V(time.prox.graph),as.character(offset:table.breaks[n + 2]))

  ## Identify subgraphs of density-based nearest neighbor graph
  subgraphs.ls <- decompose.graph(time.prox.graph)

  ## Set index for while loop
  y <- 0

  while (length(subgraphs.ls) >= 2) {
    y <- y + 1
    print(y)
    ## Compute inter-subgraph edges
    if (length(get.edgelist(subgraphs.ls[[2]])) != 0) { #decompose.graph() returns list of two, second one empty, when graph is fully connected
      print("in loop")
      ## Change row names of cluster medians df to have correct indexing
      rownames(clusters) <- c(offset:table.breaks[n + 2])

      ## Make edge.list a df and change colnames
      edge.list <- as.data.frame(edge.list)
      colnames(edge.list) <- c("Vertex1", "Vertex2", "Distance")

      ## Compute inter-subgraph edges
      i <- 1
      j <- 1
      to.add.df <- data.frame()
      for (p in 1:length(subgraphs.ls)) {
        si.graph <- subgraphs.ls[[p]]
        print(p)
        ## Prepare df to hold 1 nn from this subraph to other subgraphs
        inter.sub.nns.df <- data.frame()
        for (m in 1:length(subgraphs.ls)){
          if (m > length(subgraphs.ls)){
            break
          }
          if (m == p){
            next
          }
          sj.graph <- subgraphs.ls[[m]]
          if (distance.metric == 'manhattan') {
            print("Manhattan distance no longer supported. Using euclinean distance.")
            inter.sub.nns <- RANN::nn2(data = clusters[V(sj.graph)$name,],
                                       query=clusters[V(si.graph)$name,],
                                       k=1, searchtype="priority", eps=0.1)
          } else if (distance.metric == 'euclidean') {
            inter.sub.nns <- RANN::nn2(data = clusters[V(sj.graph)$name,],
                                       query=clusters[V(si.graph)$name,],
                                       k=1, searchtype="priority", eps=0.1)
          }
          ## Add nn output to temp df, then bind to full list
          temp.inter.sub.nns.df <- data.frame("Vertex1" = V(si.graph)$name)
          temp.inter.sub.nns.df <- cbind(temp.inter.sub.nns.df,
                                         as.data.frame(V(sj.graph)$name[as.numeric(inter.sub.nns$nn.idx)]),
                                         as.data.frame(inter.sub.nns$nn.dists))
          colnames(temp.inter.sub.nns.df) <- c("Vertex1", "Vertex2", "Distance")
          inter.sub.nns.df <-  rbind(inter.sub.nns.df, temp.inter.sub.nns.df)
        }
        ## Put list of nearest neighbors in order with shortest dist at top
        inter.sub.nns.df <- inter.sub.nns.df[order(inter.sub.nns.df$Distance), ]
        to.add.df <- rbind(to.add.df, inter.sub.nns.df[1,])

      }
      ## Add new edges and rebuild full graph
      #TODO make this and following three sections into a function, repeated below
      output.graph.update <- graph.empty()
      to.add.mat <- as.matrix(data.frame(lapply(data.frame(to.add.df),
                                                function(x) as.character(x))))
      edge.list.mat <- as.matrix(edge.list)
      vertices.edges <- rbind(edge.list.mat[,1:2], to.add.mat[,1:2])
      output.graph.update <- graph_from_edgelist(vertices.edges, directed = FALSE)

      ## Add edge weights to graph
      ogw.df <- data.frame()
      ogw.df <- data.frame(edge.list[,3])
      ta.df <- data.frame(to.add.df[,3])
      colnames(ogw.df)[1] <- "weight"
      colnames(ta.df)[1] <- "weight"
      ogw.ta.df <- data.frame(rbind(ogw.df, ta.df))
      E(output.graph.update)$weight <- ogw.ta.df[,1]
      ## Add vertex names
      output.graph.update <- set.vertex.attribute(output.graph.update,'name', #TODO WHY IS THIS DIFFERENT FROM VERTEXT NAMES AFTER WHILE LOOP????????? *************
                                                  index=V(output.graph.update),
                                                  as.character(offset:table.breaks[n + 2]))
      ## Make new edgelist
      edge.list <- data.frame(vertices.edges)
      edge.list$weights <- ogw.ta.df[,1]

      ## Test for connectedness of new graph
      subgraphs.ls <- decompose.graph(output.graph.update)
      subgraphs.el.ls <- list()
      for (ind in 1:length(subgraphs.ls)) {
        subgraphs.el.ls[[ind]] <- get.edgelist(subgraphs.ls[[ind]])
      }
    }#if
    else { break }
  }#while
  output.graph.update <- graph.empty() #TODO <-make this and following three sections into a function, repeated above <OR JUST REMOVE THIS SECTION?>
  og.el <- get.edgelist(output.graph)
  to.add.mat <- as.matrix(data.frame(lapply(data.frame(to.add.df),
                                            function(x) as.character(x))))
  vertices.edges <- rbind(og.el[,1:2], to.add.mat[,1:2])
  output.graph.update <- graph_from_edgelist(vertices.edges, directed = FALSE)

  ## Add edge weights to graph
  ogw.df <- data.frame()
  ogw.df <- data.frame(E(output.graph)$weight)
  ta.df <- data.frame(to.add.df[,3])
  colnames(ogw.df)[1] <- "weight"
  colnames(ta.df)[1] <- "weight"
  ogw.ta.df <- data.frame(rbind(ogw.df, ta.df))
  E(output.graph.update)$weight <- ogw.ta.df[,1]

  ##Add vertex names
  V(output.graph.update)$name <- 1:length(V(output.graph.update))


  ##return edge.list and output.graph as results
  colnames(edge.list) <- c("row.inds", "col.inds", "values")
  results <- list('output.graph' = output.graph.update, 'edgelist.with.distances' = edge.list)
  return(results)

}

DrawNormalizedEdgesKnn <- function(output.graph, nn.ids.df, nn.dists.df,
                                   normalized.densities, n, table.breaks,
                                   offset) {
  final.edgelist.with.distances <- c()
  for (i in names(normalized.densities)) {
    ## New not relying on cluster dist matrix
    tmp.edgelist <- cbind(as.numeric(i), as.numeric(nn.ids.df[i,]),as.numeric(nn.dists.df[i,]))[1:normalized.densities[i], ]
    final.edgelist.with.distances <- rbind(final.edgelist.with.distances, tmp.edgelist)
  }
  colnames(final.edgelist.with.distances) <- c("row.inds", "col.inds", "values")

  ## Make df with vertex names, make graph from edge.list df
  if (offset == 0) {
    output.graph.update <- graph.empty()
    vertices.edges <- final.edgelist.with.distances[,1:2]
    vertices.edges <- as.matrix(vertices.edges)
    output.graph.update <- graph_from_edgelist(vertices.edges, directed = FALSE)
    E(output.graph.update)$weight <- final.edgelist.with.distances[,3]
  } else {
    output.graph.update <- graph.empty()
    vertices.edges <- get.edgelist(output.graph)
    vertices.edges <- rbind(vertices.edges, final.edgelist.with.distances[,1:2])
    vertices.edges <- as.matrix(data.frame(lapply(data.frame(vertices.edges), function(x) as.numeric(as.character(x)))))
    output.graph.update <- graph_from_edgelist(vertices.edges, directed = FALSE)
    ogw.df <- data.frame()
    ogw.df <- data.frame(E(output.graph)$weight)
    fewd.df <- data.frame(final.edgelist.with.distances[,3])
    colnames(ogw.df)[1] <- "weight"
    colnames(fewd.df)[1] <- "weight"
    ogw.fewd.df <- data.frame(rbind(ogw.df, fewd.df))
    E(output.graph.update)$weight <- ogw.fewd.df[,1]
    V(output.graph.update)$name <- 1:length(V(output.graph.update))
  }
  results <- list('output.graph' = output.graph.update,
                  'final.edgelist.with.distances' = final.edgelist.with.distances)
  return(results)
}

BaseBuildKNN <- function(clusters, table.breaks, offset, n,
                         output.graph, k, min, max, distance.metric) {

  ## Find k+1 nearest neighbors
  cat(paste0("Clusters length:", as.character(dim(clusters)), "\n"))
  if (distance.metric == 'manhattan') {
    cat("Manhattan distance no longer supported. Using euclinean distance.\n")
    nns <- RANN::nn2(data=clusters, k=k+1, searchtype="priority", eps=0.1)
  } else if (distance.metric == 'euclidean') {
    nns <- RANN::nn2(data=clusters, k=k+1, searchtype="priority", eps=0.1)
  }
  temp_nnids.df <- as.data.frame(nns$nn.idx)
  temp_nndists.df <- as.data.frame(nns$nn.dists)
  nn.ids.df <- temp_nnids.df[,2:length(temp_nnids.df)]
  nn.dists.df <- temp_nndists.df[,2:length(temp_nndists.df)]

  ## Correctly index the edges at each sequential step
  if (offset > 0) {
    cat(paste0("In BaseBuildKNN, offset = ", as.character(offset), "\n"))
    cat(paste0("In BaseBuildKNN, table.breaks[n + 2] = ", as.character(table.breaks[n + 2]), "\n"))
    row.names(nn.ids.df) <- offset:table.breaks[n + 2]
    nn.ids.df[] <- lapply(nn.ids.df, function(x) x+table.breaks[n])
    row.names(nn.dists.df) <- offset:table.breaks[n + 2]
  }
  numcluster <- nrow(clusters)

  ## Calculate density at each vertex based on kNN
  normalized.densities <- KnnDensity(k=k, min, max, n=n,
                                     nn.ids.df = nn.ids.df,
                                     nn.dists.df = nn.dists.df,
                                     numcluster = numcluster,
                                     table.breaks = table.breaks,
                                     offset = offset)

  ## Build new edgelist with N edges for each cluster based on normalized density
  results <- DrawNormalizedEdgesKnn(output.graph = output.graph,
                                    nn.ids.df = nn.ids.df,
                                    nn.dists.df = nn.dists.df,
                                    normalized.densities = normalized.densities,
                                    n = n, table.breaks = table.breaks, offset = offset)
  output.graph <- results$output.graph
  edgelist.save <- results$final.edgelist.with.distances

  ## Test whether graph is connected, connect by kNN if not
  if (!is_connected(output.graph)) {
    print("Graph has disconnected components!")
    if (offset == 0) {
      results <- FirstConnectSubgraphs(output.graph = output.graph,
                                       edge.list = edgelist.save,
                                       offset = offset,
                                       table.breaks = table.breaks,
                                       n = n,
                                       distance.metric = distance.metric,
                                       clusters = clusters)
    } else if (offset > 0) {
      results <- ConnectSubgraphs(output.graph = output.graph,
                                  edge.list = edgelist.save[[as.character(n)]],
                                  offset = offset,
                                  table.breaks = table.breaks,
                                  n = n,
                                  distance.metric = distance.metric,
                                  clusters = clusters)
    }
    output.graph <- results$output.graph
    edgelist.save <- results$edgelist.with.distances

  }

  return(list("output.graph" = output.graph,
              "edgelist.save" = edgelist.save,
              "indexes" = nn.ids.df, #TODO remove knn outputs, no longer needed for umap
              "distances" = nn.dists.df))
}

BuildFirstFLOWMAPkNN <- function(FLOWMAP.clusters, k, min, max, distance.metric,
                                 clustering.var) {

  ## This section creates a flowmap for the first time point
  cat("Building first FLOWMAP:\n")

  ## Read in info from FLOWMAP.clusters
  table.lengths <- FLOWMAP.clusters$table.lengths
  table.breaks <- c(0, FLOWMAP.clusters$table.breaks)
  clusters <- FLOWMAP.clusters$full.clusters[1:table.lengths[1], ]
  clusters <- subset(clusters, select = clustering.var)
  output.graph <- graph.empty()

  results <- BaseBuildKNN(clusters=clusters,table.breaks=table.breaks, offset=0, n=0,
                          output.graph=output.graph, k=k, min=min, max=max,
                          distance.metric=distance.metric)

  return(results)
}

BuildFLOWMAPkNN <- function(FLOWMAP.clusters, k, min, max,
                            distance.metric, clustering.var) {

  ## Build graph from first timepoint
  first.results <- BuildFirstFLOWMAPkNN(FLOWMAP.clusters = FLOWMAP.clusters,
                                        k = k, min = min, max = max,
                                        distance.metric = distance.metric,
                                        clustering.var = clustering.var)
  output.graph <- first.results$output.graph
  edgelist.save <- list()
  edgelist.save[["0"]] <- first.results$edgelist.save
  #TODO remove all this nn index and dist stuff
  knn.indexes <- data.frame(cbind(first.results$indexes,first.results$indexes)) #to make cols for adding connections back from second timepoint
  knn.distances <- data.frame(cbind(first.results$distances,first.results$distances))

  ## Read in info for indexing timepoint clusters from FLOWMAP.clusters
  table.lengths <- FLOWMAP.clusters$table.lengths
  table.breaks <- c(0, FLOWMAP.clusters$table.breaks)

  ## This section builds the flowmap one timepoint at a time for subsequent timepoints

  for (n in 1:(length(FLOWMAP.clusters$cluster.medians) - 1)) {
    ## Offset value is used to correctly index the edges at each sequential step
    offset <- table.breaks[n] + 1
    cat(paste0("Starting next round, offset = ", as.character(offset), "\n"))

    ## Go through sequential cluster sets, add edges for each a-a and a-a+1 set
    cat(paste0("Building FLOWMAP from", n, "to", n + 1, "\n"))
    ## Get clusters for time a and a+1
    clusters <- rbind(FLOWMAP.clusters$cluster.medians[[n]], FLOWMAP.clusters$cluster.medians[[n + 1]])
    clusters <- subset(clusters, select = clustering.var)

    results <- BaseBuildKNN(clusters=clusters, table.breaks=table.breaks, offset=offset,
                            output.graph=output.graph,n=n, k=k, min=min, max=max,
                            distance.metric=distance.metric)

    output.graph <- results$output.graph
    edgelist.save[[n]] <- results$edgelist.save

    #TODO remove all this nn index and dist stuff
    knn.indexes[offset:table.breaks[n+1],11:20] <- results$indexes[offset:table.breaks[n+1],] #data.frame(cbind(knn.indexes, ))
    knn.indexes <- data.frame(dplyr::bind_rows(knn.indexes, results$indexes[(table.breaks[n+1]+1):table.breaks[n+2],]))
    knn.distances[offset:table.breaks[n+1],11:20] <- results$distances[offset:table.breaks[n+1],] #data.frame(cbind(knn.distances, ))
    knn.distances <- data.frame(dplyr::bind_rows(knn.distances, results$distances[(table.breaks[n+1]+1):table.breaks[n+2],]))
  }
  output.graph.final <- graph.empty()
  vertices.edges <- get.edgelist(output.graph)
  vertices.edges <- as.matrix(data.frame(lapply(data.frame(vertices.edges), function(x) as.numeric(as.character(x)))))
  output.graph.final <- graph_from_edgelist(vertices.edges, directed = FALSE)
  E(output.graph.final)$weight <- unlist(lapply(E(output.graph)$weight, function(x) 1/as.numeric(x)))
  output.graph.final <- AnnotateGraph(output.graph = output.graph.final,
                                      FLOWMAP.clusters = FLOWMAP.clusters)
  ## Remove duplicates and self-edges
  output.graph.final <- simplify(output.graph.final)

  return(list("output.graph" = output.graph.final,
              "edgelist.save" = edgelist.save,
              "knn.indexes" = knn.indexes, #TODO remove knn outputs, no longer needed for umap
              "knn.distances" = knn.distances))
}
