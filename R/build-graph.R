# #############################################################################
# ##FUNCTIONS FOR KNN-BASED DENSITY GRAPH BUILDING====
# #############################################################################

#' @importFrom BBmisc normalize
# kNN density function ====
KnnDensity <- function(k, min, max, n, nn.ids.df, nn.dists.df,
                       numcluster = numcluster, #table.lengths
                       table.breaks, offset) {
  all.densities.list <- list()
  for(i in 1:nrow(nn.ids.df)) {
    ## Alternatives for calculating density
    ## Dist to kth nearest neighbor:
    #density.by.knn <- abs(nn.dists.df[i,k])
    ## Transformed metric from X.shift paper:
    #density.by.Xshift <- (1/(n*(longest.dist^d))) * (sum(c(1:k)^d)/sum(compile.dists))^d
    ## Sum of distances from 1:kth nearest neighbor:
    #density.by.knn <- sum(abs(nn.dists.df[i,1:k]))
    ## Mean of distances from 1:kth nearest neighbor:
    density.by.knn <- sum(abs(nn.dists.df[i,1:k]))/k
    all.densities.list[[paste(i,".nn", sep='')]] <- density.by.knn
  }
  densities.df <- rbind(matrix(unlist(all.densities.list), byrow=T))
  normalized.densities <- BBmisc::normalize(densities.df, method = "range", range = c(min, max), margin = 2, on.constant = "quiet")
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
  #rownames(clusters) <- c(offset:table.breaks[n + 2])
  ## Make edge.list a df and change colnames
  edge.list <- as.data.frame(edge.list)
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
  #Remember starting nrow of edgelist to use for adding new edges to full graph later
  orig.nrow.edge.list <- nrow(edge.list)
  ## Make graph from edge.list (this is to force extra connection to be within self or neighboring time point)
  #because edgelist only contains current round, while graph has accumulation of all prior rounds
  time.prox.graph <- graph.empty()
  edge.list$id_num <- edge.list$id-offset+1 #have to number from 1 so igraph doesnt count all numbers up to vertex "name" number
  edge.list$index_num <- edge.list$index-offset+1 #have to number from 1 so igraph doesnt count all numbers up to vertex "name" number
  time.prox.graph <- graph_from_edgelist(as.matrix(cbind(edge.list$id_num, edge.list$index_num)), directed = FALSE)
  E(time.prox.graph)$weight <- edge.list$dist
  time.prox.graph <- set.vertex.attribute(time.prox.graph,'name',index=V(time.prox.graph),as.character(offset:table.breaks[n + 2])) #or, edge.list$id

  ## Identify subgraphs of density-based nearest neighbor graph
  subgraphs.ls <- decompose.graph(time.prox.graph)
  ## Set index for while loop
  y <- 0
  to.add.df <- NA #WORKAROUND FOR WHEN WHILE LOOP IS ENTERED BUT GRAPH IS ACTUALLY CONNECTED
  while (length(subgraphs.ls) >= 2) {
    y <- y + 1
    print(y)
    ## Compute inter-subgraph edges
    if (length(get.edgelist(subgraphs.ls[[2]])) != 0) { #decompose.graph() returns list of two, second one empty, when graph is fully connected
      print("in loop")
      ## Change row names of cluster medians df to have correct indexing
      rownames(clusters) <- c(offset:table.breaks[n + 2])

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
          temp.inter.sub.nns.df <- data.frame("id" = V(si.graph)$name)
          temp.inter.sub.nns.df <- cbind(temp.inter.sub.nns.df,
                                         as.data.frame(V(sj.graph)$name[as.numeric(inter.sub.nns$nn.idx)]),
                                         as.data.frame(inter.sub.nns$nn.dists))
          colnames(temp.inter.sub.nns.df) <- c("id", "index", "dist")
          inter.sub.nns.df <-  rbind(inter.sub.nns.df, temp.inter.sub.nns.df)
        }#for loop
        ## Put list of nearest neighbors in order with shortest dist at top
        inter.sub.nns.df <- inter.sub.nns.df[order(inter.sub.nns.df$dist), ]
        to.add.df <- rbind(to.add.df, inter.sub.nns.df[1,])

      }#for loop
      ## Add new edges to edgelist and rebuild full graph
      #TODO make this and following three sections into a function, repeated below
      to.add.df <- data.frame(lapply(data.frame(to.add.df),function(x) as.numeric(as.character(x))))
      to.add.df$id_num <- to.add.df$id-offset+1 #have to number from 1 so igraph doesnt count all numbers up to vertex "name" number
      to.add.df$index_num <- to.add.df$index-offset+1 #have to number from 1 so igraph doesnt count all numbers up to vertex "name" number
      edge.list <- rbind(edge.list, to.add.df)
      ## rebuild graph from new updated edgelist
      output.graph.update <- graph.empty()
      output.graph.update <- graph_from_edgelist(as.matrix(edge.list[,c("id_num", "index_num")]), directed = FALSE)
      ## Add edge weights to graph
      E(output.graph.update)$weight <- edge.list$dist 
      ## Add vertex names to edgelist (still only has the two timepoints)
      output.graph.update <- set.vertex.attribute(output.graph.update,'name',
                                                  index=V(output.graph.update),
                                                  as.character(offset:table.breaks[n + 2]))

      ## Test for connectedness of new graph
      subgraphs.ls <- decompose.graph(output.graph.update)
      subgraphs.el.ls <- list()
      for (ind in 1:length(subgraphs.ls)) {
        subgraphs.el.ls[[ind]] <- get.edgelist(subgraphs.ls[[ind]])
      }
    }#if
    else { 
      break }
  }#while
  if (!is.na(to.add.df)) { #WORKAROUND FOR WHEN WHILE LOOP IS ENTERED BUT GRAPH IS ACTUALLY CONNECTED
    ## This section adds new edges to full graph
    output.graph.update <- graph.empty()
    og.el <- get.edgelist(output.graph)
    colnames(og.el) <- c("id", "index")
    #combine original edges in graph with new edges added here
    og.el <- rbind(og.el[,1:2], edge.list[(orig.nrow.edge.list+1):nrow(edge.list),c("id", "index")]) #can use actual numbers bc this has from 1:highest anyway
    output.graph.update <- graph_from_edgelist(as.matrix(og.el), directed = FALSE)
    
    ## Add corresponding new edge weights to full graph
    ogu.weights <- c(E(output.graph)$weight, edge.list[(orig.nrow.edge.list+1):nrow(edge.list),"dist"])#combine original weights in graph with new weights added here
    E(output.graph.update)$weight <- ogu.weights
    
    ##Add vertex names to graph containing all clusters from all timepoints so far
    V(output.graph.update)$name <- 1:length(V(output.graph.update))
    
    output.graph <- output.graph.update
  }
  ##return edge.list and output.graph as results
  edgelist.with.distances <- edge.list[,c("id", "index","dist")]
  colnames(edgelist.with.distances) <- c("row.inds", "col.inds", "values")
  results <- list('output.graph' = output.graph, 'edgelist.with.distances' = edge.list)
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
    # Set density-based edge number and format as igraph edgelist
    tmp.edgelist <- data.frame(cbind(as.numeric(i), as.numeric(nn.ids.df[i,]),as.numeric(nn.dists.df[i,])))#[1:normalized.densities[i],]
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
  cat(paste0("Clusters length:", as.character(dim(clusters)[1]), "\n"))
  if (k < max) {
    cat(paste0("Input k (",k,") less than max edge number (",max,") -- changing k to max edge number"))
    k <- max
  }
  if (distance.metric == 'manhattan') {
    cat("Manhattan distance no longer supported. Using euclinean distance.\n")
    nns <- RANN::nn2(data=clusters, k=k+1, searchtype="priority", eps=0.1) #k=k+1
  } else if (distance.metric == 'euclidean') {
    nns <- RANN::nn2(data=clusters, k=k+1, searchtype="priority", eps=0.1) #k=k+1
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
    knn.indexes[[n]] <- data.frame(cbind(knn.indexes[[n]],results$indexes[1:table.lengths[n],])) #to make cols for adding connections back from second timepoint
    knn.distances[[n]] <- data.frame(cbind(knn.distances[[n]],results$distances[1:table.lengths[n],])) #(nrow(results$indexes)/2) dplyr::bind_cols

    #Then add for cells in new timepoint
    #knn.indexes[[n+1]] <- results$indexes[1:(nrow(results$indexes)/2),] 
    knn.indexes[[n+1]] <- results$indexes[(table.lengths[n]+1):nrow(results$indexes),]
    knn.distances[[n+1]] <- results$distances[(table.lengths[n]+1):nrow(results$distances),]
  }
  # PREP KNN FOR OUTPUT ####################################### (should probably be simplified in the future, but works so ok for now)
  #make final timepoint data frame as well
  knn.indexes[[n+1]] <- data.frame(knn.indexes[[n+1]]) #make final timepoint a dataframe so it doesnt mess up below processing
  knn.distances[[n+1]] <- data.frame(knn.distances[[n+1]])
  ## Remove duplicate connections====
  keep_ind = lapply(knn.indexes, function(y) apply(y,1,function(x){
    to.keep <- c()
    for (i in unique(x)) {
      to.keep <- rbind(to.keep, head(which(x==i),1))
    }
    to.keep
  }))

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
  print("processing knn output")
  for(i in 1:ncol(knn.out$indexes)) knn.out$indexes[sapply(knn.out$indexes[,i], is.null),i] <- 0 #replace null index with 0 (so 'unlist' maintains full dim)
  for(i in 1:ncol(knn.out$distances)) knn.out$distances[sapply(knn.out$distances[,i], is.null),i] <- 10000 #replace null dist with high dist
  knn.out$indexes  <-  as.data.frame(matrix(unlist(knn.out$indexes), nrow=length(unlist(knn.out$indexes[1])))) #unnest 
  knn.out$distances  <-  as.data.frame(matrix(unlist(knn.out$distances), nrow=length(unlist(knn.out$distances[1])))) #unnest 
  dist_order = t(apply(knn.out$distances, 1, order))
  for(i in 1:nrow(knn.out$indexes)) knn.out$indexes[i,] <- knn.out$indexes[i,dist_order[i,]] #order indexes by increasing distance (necessary to incorporate second timepoint connections into overall order)
  for(i in 1:nrow(knn.out$distances)) knn.out$distances[i,] <- knn.out$distances[i,][dist_order[i,]] #increasing order distances
  knn.out$indexes  <-  as.data.frame(matrix(unlist(knn.out$indexes), nrow=length(unlist(knn.out$indexes[1])))) 
  knn.out$distances  <-  as.data.frame(matrix(unlist(knn.out$distances), nrow=length(unlist(knn.out$distances[1]))))
  print("finished processing knn output")
  
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
  
  
  return(list("output.graph" = output.graph.final.simplified,
              "knn.out" = knn.out))
              
}
