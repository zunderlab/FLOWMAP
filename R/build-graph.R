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
      #TODO WHY IS THIS DIFFERENT FROM VERTEX NAMES AFTER WHILE LOOP????????? *************
      output.graph.update <- set.vertex.attribute(output.graph.update,'name',
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
  #TODO <-make this and following three sections into a function, repeated above <OR JUST REMOVE THIS SECTION?>
  output.graph.update <- graph.empty()
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
  knn.indexes <- list()
  knn.distances <- list()
  for (i in names(normalized.densities)) {
    ## New not relying on cluster dist matrix
    global.normalized.densities <<- normalized.densities
    global.nn.ids.df <<- nn.ids.df
    global.nn.dists.df <<- nn.dists.df
    # Set density-based edge number and format as igraph edgelist
    tmp.edgelist <- data.frame(cbind(as.numeric(i), as.numeric(nn.ids.df[i,]),as.numeric(nn.dists.df[i,])))#[1:normalized.densities[i],]
    global.tmp.edgelist <<- tmp.edgelist
    colnames(tmp.edgelist) <- c("id","index", "dist")
    tmp.el.ord <- tmp.edgelist[with(tmp.edgelist, order(dist)), ] #need this because keeping all edges, just expanding distance for ones that originally would have been dropped ##############################################################################################
    if (nrow(tmp.el.ord) > normalized.densities[i]) {
      tmp.el.ord$dist[(normalized.densities[i]+1):nrow(tmp.el.ord)] <- (tmp.el.ord$dist[(normalized.densities[i]+1):nrow(tmp.el.ord)])^4
    }
    tmp.el.ord <- tmp.el.ord[complete.cases(tmp.el.ord), ]##############################################################################################
    final.edgelist.with.distances <- rbind(final.edgelist.with.distances, tmp.el.ord)##############################################################################################
    #final.edgelist.with.distances <- rbind(final.edgelist.with.distances, tmp.edgelist)
    
    knn.indexes <- rbind(knn.indexes,tmp.el.ord$index)
    knn.distances <- rbind(knn.distances,tmp.el.ord$dist)
  }
  #colnames(final.edgelist.with.distances) <- c("row.inds", "col.inds", "values")
  print("final")
  ## Make df with vertex names, make graph from edge.list df
  if (offset == 0) {
    output.graph.update <- graph.empty()
    vertices.edges <- final.edgelist.with.distances[,1:2]
    vertices.edges <- as.matrix(vertices.edges)
    output.graph.update <- graph_from_edgelist(vertices.edges, directed = FALSE)
    E(output.graph.update)$weight <- final.edgelist.with.distances[,3]
  } else {
    output.graph.update <- graph.empty()
    #vertices.edges <- get.edgelist(output.graph)
    #vertices.edges <- rbind(vertices.edges, final.edgelist.with.distances[,1:2])
    vertices.edges <- data.frame(get.edgelist(output.graph))##############################################################################################
    vertices.edges <- data.table::rbindlist(list(vertices.edges, final.edgelist.with.distances[,1:2]), use.names = FALSE)##############################################################################################
    print("vert edge rbind")
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
                  'edgelist.with.distances' = final.edgelist.with.distances,
                  'indexes' = knn.indexes,
                  'distances' = knn.distances)
  return(results)
}

BaseBuildKNN <- function(clusters, table.breaks, offset, n, k, min, max, 
                         output.graph, edgelist.save, distance.metric) {

  ## Find k+1 nearest neighbors
  cat(paste0("Clusters length:", as.character(dim(clusters)), "\n"))
  if (k < max) {
    cat(paste0("Input k (",k,") less than max edge number (",max,") -- changing k to max edge number"))
    k <- max
  }
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
  edgelist.save <- results$edgelist.with.distances
  knn.indexes <- results$indexes
  knn.distances <- results$distances
  
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
                                  edge.list = edgelist.save,
                                  offset = offset,
                                  table.breaks = table.breaks,
                                  n = n,
                                  distance.metric = distance.metric,
                                  clusters = clusters)
    }
    output.graph <- results$output.graph
    edgelist.save <- results$edgelist.with.distances
  }#end if
  
  return(list("output.graph" = output.graph,
              'indexes' = knn.indexes,
              'distances' = knn.distances))
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
  global.output.graph.first <<- output.graph
  knn.indexes <- list()
  knn.distances <- list()
  knn.indexes[[1]] <- first.results$indexes #
  knn.distances[[1]] <- first.results$distances #

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
    #output.graph <- simplify(output.graph) ### TEMP FOR TESTING!!!!!
    
    #First fill in new connections from cells in previous timepoint
    knn.indexes[[n]] <- data.frame(cbind(knn.indexes[[n]],results$indexes[1:(nrow(results$indexes)/2),])) #to make cols for adding connections back from second timepoint
    knn.distances[[n]] <- data.frame(cbind(knn.distances[[n]],results$distances[1:(nrow(results$distances)/2),]))
    #Then add for cells in new timepoint
    #knn.indexes[[n+1]] <- results$indexes[1:(nrow(results$indexes)/2),] 
    knn.indexes[[n+1]] <- results$indexes[(nrow(results$indexes)/2+1):nrow(results$indexes),]
    knn.distances[[n+1]] <- results$distances[(nrow(results$distances)/2+1):nrow(results$distances),]
  }
  # PREP KNN FOR OUTPUT ####################################### (should probably be simplified in the future, but works so ok for now)
  #make final timepoint data frame as well
  knn.indexes[[n+1]] <- data.frame(knn.indexes[[n+1]]) #make final timepoint a dataframe so it doesnt mess up below processing
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
  print("processing knn output")
  for(i in 1:ncol(knn.out$indexes)) knn.out$indexes[sapply(knn.out$indexes[,i], is.null),i] <- 0 #replace null index with 0 (so 'unlist' maintains full dim)
  for(i in 1:ncol(knn.out$distances)) knn.out$distances[sapply(knn.out$distances[,i], is.null),i] <- 10000 #replace null dist with high dist
  knn.out$indexes  <-  as.data.frame(matrix(unlist(knn.out$indexes), nrow=length(unlist(knn.out$indexes[1])))) #unnest 
  knn.out$distances  <-  as.data.frame(matrix(unlist(knn.out$distances), nrow=length(unlist(knn.out$distances[1])))) #unnest 
  for(i in 1:nrow(knn.out$indexes)) knn.out$indexes[i,] <- knn.out$indexes[i,order(knn.out$distances[i,])] #order indexes by increasing distance (necessary to incorporate second timepoint connections into overall order)
  for(i in 1:nrow(knn.out$distances)) knn.out$distances[i,] <- knn.out$distances[i,][order(knn.out$distances[i,])] #increasing order distances
  knn.out$indexes  <-  as.data.frame(matrix(unlist(knn.out$indexes), nrow=length(unlist(knn.out$indexes[1])))) 
  knn.out$distances  <-  as.data.frame(matrix(unlist(knn.out$distances), nrow=length(unlist(knn.out$distances[1]))))
  print("finished processing knn output")
  global.knn.out <<- knn.out
  
  # PREP GRAPH FOR OUTPUT #######################################
  output.graph.final <- graph.empty()
  vertices.edges <- get.edgelist(output.graph)
  vertices.edges <- as.matrix(data.frame(lapply(data.frame(vertices.edges), function(x) as.numeric(as.character(x)))))
  output.graph.final <- graph_from_edgelist(vertices.edges, directed = FALSE)
  E(output.graph.final)$weight <- unlist(lapply(E(output.graph)$weight, function(x) 1/as.numeric(x)))
  output.graph.final <- AnnotateGraph(output.graph = output.graph.final,
                                      FLOWMAP.clusters = FLOWMAP.clusters)
  print("3")
  ## Remove duplicates and self-edges
  output.graph.final.simplified <- simplify(output.graph.final)
  global.output.graph.final <<- output.graph.final
  global.output.graph.simplified <<- output.graph.final.simplified
  
  
  return(list("output.graph" = output.graph.final.simplified,
              "knn.out" = knn.out))
              
}
