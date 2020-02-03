# #############################################################################
# ##FUNCTIONS FOR KNN-BASED DENSITY GRAPH BUILDING====
# #############################################################################

## kNN density function ====
KnnDensity <- function(k, min, max, n, nn.ids.df, nn.dists.df,
                       numcluster = numcluster, #table.lengths
                       table.breaks, offset) {
  all.densities.list <- list()
  for(i in 1:nrow(nn.ids.df)) {
    ## Alternatives for calculating density
    ##Dist to kth nearest neighbor
    density.by.knn <- abs(nn.dists.df[i,k])
    ##Transformed metric from X.shift paper
    #density.by.Xshift <- (1/(n*(longest.dist^d))) * (sum(c(1:k)^d)/sum(compile.dists))^d
    ##Sum of distances from 1:kth nearest neighbor
    #density.by.knn <- sum(abs(nn.dists.df[i,1:k]))
    ##Mean of distances from 1:kth nearest neighbor
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
  #global.densities.ls[[n+1]] <<- normalized.densities
  return(normalized.densities)
}

##Check connectedness for first flowmap ====
FirstConnectSubgraphs <- function(output.graph, edge.list, offset,
                                  table.breaks, n, distance.metric, clusters) {
  #global.clusters.first <<- clusters
  output.graph <- set.vertex.attribute(output.graph,'name',index=V(output.graph),as.character(1:vcount(output.graph)))
  ##Identify subgraphs of density-based nearest neighbor graph
  subgraphs.ls <- decompose.graph(output.graph)
  #global.subgraph.ls.first <<- subgraphs.ls
  ##change row names of cluster medians df to have correct indexing
  rownames(clusters) <- c(offset:table.breaks[n + 2])
  ## make edge.list a df and change colnames
  edge.list <- as.data.frame(edgelist.save)
  colnames(edge.list) <- c("Vertex1", "Vertex2", "Distance")
  #Set lists and ind for while loop
  y <- 0
  graphs.ls <- list()
  subgraphs.ls.el <- list()
  while (length(subgraphs.ls) >= 2) {
    y <- y + 1
    print(y)
    ##Compute inter-subgraph edges
    if (length(get.edgelist(subgraphs.ls[[2]])) != 0) {
      print("in loop")
      i <- 1
      j <- 1
      to.add.df <- data.frame()
      for (p in 1:length(subgraphs.ls)) {
        si.graph <- subgraphs.ls[[p]]
        ##Prepare df to hold 1 nn from this subraph to other subgraphs
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
          ##add nn output to temp df, then bind to full list
          temp.inter.sub.nns.df <- data.frame("Vertex1" = V(si.graph)$name)
          temp.inter.sub.nns.df <- cbind(temp.inter.sub.nns.df,
                                         as.data.frame(V(sj.graph)$name[as.numeric(inter.sub.nns$nn.idx)]),
                                         as.data.frame(inter.sub.nns$nn.dists))
          colnames(temp.inter.sub.nns.df) <- c("Vertex1", "Vertex2", "Distance")
          inter.sub.nns.df <-  rbind(inter.sub.nns.df, temp.inter.sub.nns.df)
        }
        ##put list of nearest neighbors in order with shortest dist at top
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

  #global.subraphs.ls.ls <<- subgraphs.ls.el
  #global.graphs.ls <<- graphs.ls

  ##return edge.list and output.graph as results
  colnames(edge.list) <- c("row.inds", "col.inds", "values")
  results <- list('output.graph' = output.graph, 'edgelist.with.distances' = edge.list)
  return(results)

}

## kNN verion of "CheckMSTEdges"
ConnectSubgraphs <- function(output.graph, edge.list, offset,
                             table.breaks, n, distance.metric, clusters) {
  # global.cs.clusters <<- clusters
  # global.cs.output.graph <<- output.graph
  # global.cs.el <<- edge.list
  # global.cs.offset <<- offset
  # global.cs.table.breaks <<- table.breaks
  # global.cs.n <<- n
  # global.cs.distance.metric <<- distance.metric
  print("in connected")
  ##make df with vertex names, make graph from edge.list df
  # v.df <- data.frame(names = (offset + 1):table.breaks[n + 2])
  # time.prox.graph <- graph.data.frame(edge.list[,1:2], vertices = v.df)
  # E(time.prox.graph)$weight <- edge.list[,3]
  # testing.time.prox.graph.pre <- time.prox.graph
  # testing.is.conn.pre <- is.connected(time.prox.graph)

  time.prox.graph <- graph.empty()
  # vertices.edges <- get.edgelist(output.graph)
  # vertices.edges <- rbind(vertices.edges, final.edgelist.with.distances[,1:2])
  # vertices.edges <- as.matrix(data.frame(lapply(data.frame(vertices.edges),
  #                                     function(x) as.numeric(as.character(x)))))
  time.prox.graph <- graph_from_edgelist(edge.list[,1:2], directed = FALSE)
  E(time.prox.graph)$weight <- edge.list[,3]
  time.prox.graph <- set.vertex.attribute(time.prox.graph,'name',index=V(time.prox.graph),as.character(offset:table.breaks[n + 2]))
  # E(output.graph.update)$weight <- rbind(unlist(lapply(E(output.graph)$weight,
  #                                  function(x) 1/as.numeric(x))), final.edgelist.with.distances[,3])

  ##Identify subgraphs of density-based nearest neighbor graph
  subgraphs.ls <- decompose.graph(time.prox.graph)
  #global.subgraphs.ls <<- subgraphs.ls
  #Set lists and ind for while loop
  y <- 0

  while (length(subgraphs.ls) >= 2) {
    y <- y + 1
    print(y)
    ##Compute inter-subgraph edges
    if (length(get.edgelist(subgraphs.ls[[2]])) != 0) {
      print("in loop")
      # if(length(subgraphs.ls) <= 1) {
      #   print("huh?")
      #   results <- list('output.graph' = output.graph, 'final.edgelist.with.distances' = edge.list)
      #   return(results)
      # }
      # ##Test if true subgraphs ====
      # else if (length(get.edgelist(subgraphs.ls[[2]])) >= 1) {
      ##change row names of cluster medians df to have correct indexing
      rownames(clusters) <- c(offset:table.breaks[n + 2])

      ## make edge.list a df and change colnames
      edge.list <- as.data.frame(edge.list)
      colnames(edge.list) <- c("Vertex1", "Vertex2", "Distance")

      ##Compute inter-subgraph edges
      i <- 1
      j <- 1
      to.add.df <- data.frame()
      for (p in 1:length(subgraphs.ls)) {
        si.graph <- subgraphs.ls[[p]]
        print(p)
        ##Prepare df to hold 1 nn from this subraph to other subgraphs
        inter.sub.nns.df <- data.frame()
        for (m in 1:length(subgraphs.ls)){
          if (m > length(subgraphs.ls)){
            break
          }
          if (m == p){
            next
          }
          #print(m)
          sj.graph <- subgraphs.ls[[m]]
          if (distance.metric == 'manhattan') {
            #global.sj.graph.names <<- V(sj.graph)$name
            #global.cluster.sj.subset <<- list(data = clusters[V(sj.graph)$name,],
            #query=clusters[V(si.graph)$name,])
            print("Manhattan distance no longer supported. Using euclinean distance.")
            inter.sub.nns <- RANN::nn2(data = clusters[V(sj.graph)$name,],
                                       query=clusters[V(si.graph)$name,],
                                       k=1, searchtype="priority", eps=0.1)
          } else if (distance.metric == 'euclidean') {
            inter.sub.nns <- RANN::nn2(data = clusters[V(sj.graph)$name,],
                                       query=clusters[V(si.graph)$name,],
                                       k=1, searchtype="priority", eps=0.1)
          }
          ##add nn output to temp df, then bind to full list
          temp.inter.sub.nns.df <- data.frame("Vertex1" = V(si.graph)$name)
          temp.inter.sub.nns.df <- cbind(temp.inter.sub.nns.df,
                                         as.data.frame(V(sj.graph)$name[as.numeric(inter.sub.nns$nn.idx)]),
                                         as.data.frame(inter.sub.nns$nn.dists))
          colnames(temp.inter.sub.nns.df) <- c("Vertex1", "Vertex2", "Distance")
          inter.sub.nns.df <-  rbind(inter.sub.nns.df, temp.inter.sub.nns.df)
        }
        ##put list of nearest neighbors in order with shortest dist at top
        inter.sub.nns.df <- inter.sub.nns.df[order(inter.sub.nns.df$Distance), ]
        to.add.df <- rbind(to.add.df, inter.sub.nns.df[1,])

      }
      ##Now read in full edge list from graph so you can rebuild full graph
      # edge.list <- data.frame(get.edgelist(output.graph, names = TRUE))
      # edge.list$Distance <- E(output.graph)$weight
      # colnames(edge.list) <- c("Vertex1", "Vertex2", "Distance")
      # edge.list[,1] <- as.double(edge.list[,1])
      # edge.list[,2] <- as.double(edge.list[,2])
      # edge.list <- data.frame(rbind(edge.list, to.add.df))

      #Build full graph with new edges
      # v.df <- data.frame(names =1:table.breaks[n + 2])
      # try.fresh.graph <- graph.data.frame(edge.list[,1:2], vertices = v.df)
      # E(try.fresh.graph)$weight <- edge.list[,3]

      output.graph.update <- graph.empty()
      to.add.mat <- as.matrix(data.frame(lapply(data.frame(to.add.df),
                                                function(x) as.character(x))))
      edge.list.mat <- as.matrix(edge.list)
      vertices.edges <- rbind(edge.list.mat[,1:2], to.add.mat[,1:2])
      #global.vertices.edges <<- vertices.edges
      # vertices.edges <- as.matrix(data.frame(lapply(data.frame(vertices.edges),
      #                                 function(x) as.numeric(as.character(x)))))
      output.graph.update <- graph_from_edgelist(vertices.edges, directed = FALSE)
      # E(output.graph.update)$weight <- rbind(unlist(lapply(E(output.graph)$weight,
      #                                   function(x) 1/as.numeric(x))), to.add.df[,3])
      ##Add edge weights to graph
      ogw.df <- data.frame()
      ogw.df <- data.frame(edge.list[,3])
      ta.df <- data.frame(to.add.df[,3])
      colnames(ogw.df)[1] <- "weight"
      colnames(ta.df)[1] <- "weight"
      ogw.ta.df <- data.frame(rbind(ogw.df, ta.df))
      E(output.graph.update)$weight <- ogw.ta.df[,1]
      ##Add vertex names
      output.graph.update <- set.vertex.attribute(output.graph.update,'name',
                                                  index=V(output.graph.update),
                                                  as.character(offset:table.breaks[n + 2]))
      ##Make new edgelist
      edge.list <- data.frame(vertices.edges)
      edge.list$weights <- ogw.ta.df[,1]

      ##Test for connectedness of new graph
      subgraphs.ls <- decompose.graph(output.graph.update)
      subgraphs.el.ls <- list()
      for (ind in 1:length(subgraphs.ls)) {
        subgraphs.el.ls[[ind]] <- get.edgelist(subgraphs.ls[[ind]])
      }
      #global.graphs.ls[[y]] <<- output.graph.update
      #global.subgraphs.ls.el[[y]] <<- subgraphs.el.ls
    }#if
    else { break }
  }#while
  output.graph.update <- graph.empty()
  og.el <- get.edgelist(output.graph)
  to.add.mat <- as.matrix(data.frame(lapply(data.frame(to.add.df),
                                            function(x) as.character(x))))
  vertices.edges <- rbind(og.el[,1:2], to.add.mat[,1:2])
  #global.vertices.edges <<- vertices.edges
  # vertices.edges <- as.matrix(data.frame(lapply(data.frame(vertices.edges),
  #                                 function(x) as.numeric(as.character(x)))))
  output.graph.update <- graph_from_edgelist(vertices.edges, directed = FALSE)
  # E(output.graph.update)$weight <- rbind(unlist(lapply(E(output.graph)$weight,
  #                                   function(x) 1/as.numeric(x))), to.add.df[,3])
  ##Add edge weights to graph
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
  global.normalized.densities <<- normalized.densities
  #print(paste0("Names of normalized.densities: ", names(normalized.densities)))
  for (i in names(normalized.densities)) {
    #print(i)
    ##new not relying on cluster dist matrix
    tmp.edgelist <- cbind(as.numeric(i), as.numeric(nn.ids.df[i,]),as.numeric(nn.dists.df[i,]))[1:normalized.densities[i], ]
    final.edgelist.with.distances <- rbind(final.edgelist.with.distances, tmp.edgelist)
  }
  colnames(final.edgelist.with.distances) <- c("row.inds", "col.inds", "values")
  #global.draw.norm.edgelist <<- final.edgelist.with.distances
  ##make df with vertex names, make graph from edge.list df
  if (offset == 0) {
    # v.df <- data.frame(names = 1:table.breaks[n + 2])
    # output.graph <- graph.data.frame(final.edgelist.with.distances[,1:2], vertices = v.df)
    # E(output.graph)$weight <- final.edgelist.with.distances[,3]
    output.graph.update <- graph.empty()
    vertices.edges <- final.edgelist.with.distances[,1:2]
    vertices.edges <- as.matrix(vertices.edges)
    output.graph.update <- graph_from_edgelist(vertices.edges, directed = FALSE)
    E(output.graph.update)$weight <- final.edgelist.with.distances[,3]
  } else {
    # curr.edgelist <- cbind(get.edgelist(output.graph),E(output.graph)$weight)
    # curr.edgelist <- rbind(curr.edgelist, final.edgelist.with.distances)
    # v.df <- data.frame(names = 1:table.breaks[n + 2])
    # output.graph <- graph.data.frame(curr.edgelist[,1:2], vertices = v.df)
    # E(output.graph)$weight <- curr.edgelist[,3]
    output.graph.update <- graph.empty()
    vertices.edges <- get.edgelist(output.graph)
    vertices.edges <- rbind(vertices.edges, final.edgelist.with.distances[,1:2])
    vertices.edges <- as.matrix(data.frame(lapply(data.frame(vertices.edges), function(x) as.numeric(as.character(x)))))
    output.graph.update <- graph_from_edgelist(vertices.edges, directed = FALSE)
    #E(output.graph.update)$weight <- rbind(unlist(lapply(E(output.graph)$weight, function(x) 1/as.numeric(x))), final.edgelist.with.distances[,3])
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

  ##Find k+1 nearest neighbors
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

  #correctly index the edges at each sequential step
  if (offset > 0) {
    cat(paste0("In BaseBuildKNN, offset = ", as.character(offset), "\n"))
    cat(paste0("In BaseBuildKNN, table.breaks[n + 2] = ", as.character(table.breaks[n + 2]), "\n"))
    row.names(nn.ids.df) <- offset:table.breaks[n + 2]
    #print(row.names(nn.ids.df))
    nn.ids.df[] <- lapply(nn.ids.df, function(x) x+table.breaks[n])
    #print(row.names(nn.ids.df))
    row.names(nn.dists.df) <- offset:table.breaks[n + 2]
  }
  global.inner.knn.ls[[as.character(offset)]] <<- list(indexes = nn.ids.df,
                                                       distances = nn.dists.df)
  numcluster <- nrow(clusters)

  ##Calculate density at each vertex based on kNN
  normalized.densities <- KnnDensity(k=k, min, max, n=n,
                                     nn.ids.df = nn.ids.df,
                                     nn.dists.df = nn.dists.df,
                                     numcluster = numcluster,
                                     #table.lengths = table.lengths,
                                     table.breaks = table.breaks,
                                     offset = offset)
  global.post.kd.knn.ls[[as.character(offset)]] <<- list(indexes = nn.ids.df,
                                                         distances = nn.dists.df)

  ##Build new edgelist with N edges for each cluster based on normalized density
  results <- DrawNormalizedEdgesKnn(output.graph = output.graph,
                                    nn.ids.df = nn.ids.df,
                                    nn.dists.df = nn.dists.df,
                                    normalized.densities = normalized.densities,
                                    n = n, table.breaks = table.breaks, offset = offset)
  output.graph <- results$output.graph
  edgelist.save <- results$final.edgelist.with.distances
  global.post.nd.knn.ls[[as.character(offset)]] <<- list(indexes = nn.ids.df,
                                                         distances = nn.dists.df)

  ##Test whether graph is connected, connect by kNN if not
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
    edgelist.save <- results$edgelist.with.distances#do we need to even store edgelist at all?

  }
  global.post.connect.knn.ls[[as.character(offset)]] <<- list(indexes = nn.ids.df,
                                                              distances = nn.dists.df)
  return(list("output.graph" = output.graph,
              "edgelist.save" = edgelist.save,
              "indexes" = nn.ids.df,
              "distances" = nn.dists.df))
}

BuildFirstFLOWMAPkNN <- function(FLOWMAP.clusters, k, min, max, distance.metric,
                                 clustering.var) {

  # This section creates a flowmap for the first time point
  cat("Building first FLOWMAP:\n")

  ##Read in info from FLOWMAP.clusters
  table.lengths <- FLOWMAP.clusters$table.lengths
  table.breaks <- c(0, FLOWMAP.clusters$table.breaks)
  clusters <- FLOWMAP.clusters$full.clusters[1:table.lengths[1], ]
  clusters <- subset(clusters, select = clustering.var)
  output.graph <- graph.empty()

  results <- BaseBuildKNN(clusters=clusters,table.breaks=table.breaks, offset=0, n=0,
                          output.graph=output.graph, k=k, min=min, max=max,
                          distance.metric=distance.metric)
  global.first.results <<- results
  return(results)
}

BuildFLOWMAPkNN <- function(FLOWMAP.clusters, k, min, max,
                            distance.metric, clustering.var) {

  global.inner.knn.ls <<- list()
  global.post.kd.knn.ls <<- list()
  global.post.nd.knn.ls <<- list()
  global.post.connect.knn.ls <<- list()
  global.knn.ls <<- list()
  #build graph from first timepoint
  first.results <- BuildFirstFLOWMAPkNN(FLOWMAP.clusters = FLOWMAP.clusters,
                                        k = k, min = min, max = max,
                                        distance.metric = distance.metric,
                                        clustering.var = clustering.var)
  output.graph <- first.results$output.graph
  edgelist.save <- list()
  edgelist.save[["0"]] <- first.results$edgelist.save
  knn.indexes <- data.frame(cbind(first.results$indexes,first.results$indexes)) #to make cols for adding connections back from second timepoint
  global.knn.indexes.first <<- first.results$indexes
  knn.distances <- data.frame(cbind(first.results$distances,first.results$distances))

  global.knn.ls[["0"]] <<- list("knn.indexes"=knn.indexes,
                                "knn.distances"=knn.distances)
  ##Read in info for indexing timepoint clusters from FLOWMAP.clusters
  table.lengths <- FLOWMAP.clusters$table.lengths
  table.breaks <- c(0, FLOWMAP.clusters$table.breaks)

  # This section builds the flowmap one timepoint at a time for subsequent timepoints

  for (n in 1:(length(FLOWMAP.clusters$cluster.medians) - 1)) {
    # offset value is used to correctly index the edges at each sequential step
    offset <- table.breaks[n] + 1
    cat(paste0("Starting next round, offset = ", as.character(offset), "\n"))
    # go through sequential cluster sets, add edges for each a-a and a-a+1 set
    cat(paste0("Building FLOWMAP from", n, "to", n + 1, "\n"))

    #n_n1.table.lengths <- FLOWMAP.clusters$table.lengths[n:(n + 1)]
    #n_n1.table.breaks <- table.breaks[n:(n + 1)]

    ##get clusters for time a and a+1
    clusters <- rbind(FLOWMAP.clusters$cluster.medians[[n]], FLOWMAP.clusters$cluster.medians[[n + 1]])
    clusters <- subset(clusters, select = clustering.var)

    results <- BaseBuildKNN(clusters=clusters, table.breaks=table.breaks, offset=offset,
                            output.graph=output.graph,n=n, k=k, min=min, max=max,
                            distance.metric=distance.metric)
    output.graph <- results$output.graph
    edgelist.save[[n]] <- results$edgelist.save
    knn.indexes[offset:table.breaks[n+1],11:20] <- results$indexes[offset:table.breaks[n+1],] #data.frame(cbind(knn.indexes, ))
    knn.indexes <- data.frame(dplyr::bind_rows(knn.indexes, results$indexes[(table.breaks[n+1]+1):table.breaks[n+2],]))
    knn.distances[offset:table.breaks[n+1],11:20] <- results$distances[offset:table.breaks[n+1],] #data.frame(cbind(knn.distances, ))
    knn.distances <- data.frame(dplyr::bind_rows(knn.distances, results$distances[(table.breaks[n+1]+1):table.breaks[n+2],]))
    global.knn.ls[[n]] <<- list(knn.indexes=knn.indexes,
                                knn.distances=knn.distances)
  }
  output.graph.final <- graph.empty()
  vertices.edges <- get.edgelist(output.graph)
  vertices.edges <- as.matrix(data.frame(lapply(data.frame(vertices.edges), function(x) as.numeric(as.character(x)))))
  output.graph.final <- graph_from_edgelist(vertices.edges, directed = FALSE)
  E(output.graph.final)$weight <- unlist(lapply(E(output.graph)$weight, function(x) 1/as.numeric(x)))
  output.graph.final <- AnnotateGraph(output.graph = output.graph.final,
                                      FLOWMAP.clusters = FLOWMAP.clusters)
  ##Remove duplicates and self-edges
  output.graph.final <- simplify(output.graph.final)
  global.output.graph.final <<- output.graph.final
  global.knn.indexes <<- knn.indexes
  global.knn.distances <<- knn.distances
  global.edgelist.save <<- edgelist.save
  return(list("output.graph" = output.graph.final,
              "edgelist.save" = edgelist.save,
              "knn.indexes" = knn.indexes,
              "knn.distances" = knn.distances))
}
