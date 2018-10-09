
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

RemoveDuplicateValues <- function(arr.inds.with.values) {
  x <- arr.inds.with.values[!duplicated(arr.inds.with.values[, "values"]),]
  return(x)
}

# THIS IS REMOVING EDGES OF CELLS THAT ARE ALL THE SAME
# BECAUSE THESE CELLS ARE IDENTICAL, THEY HAVE VALUES (DISTANCE)
# ALL EQUAL TO ZERO. THIS MEANS IT ONLY KEEPS ONE OF THOSE EDGES

RemoveWithinNEdges <- function(arr.inds.with.values, inds.in.n) {
  d <- arr.inds.with.values
  d <- d[!(!is.na(match(d[, "row.inds"], inds.in.n)) &
             !is.na(match(d[, "col.inds"], inds.in.n))), ]
  inds.no.n_n <- d
  return(inds.no.n_n)
}

# FindNormalized <- function(cluster.distances.matrix, per, min,
#                            max, numcluster, table.lengths,
#                            table.breaks, offset) {
#   # make fully connected graph from adjacency list
#   # note: "weight" here is really distance, calling it weight for the mst function later needs this
#   arr.inds <- ConvertIndex(inds = order(cluster.distances.matrix), orig.data.frame = cluster.distances.matrix, keep.names = TRUE)
#   edgelist.with.distances <- MergeValues(arr.inds = arr.inds, orig.data.frame = cluster.distances.matrix)
#   edgelist.with.distances <- RemoveDuplicateValues(edgelist.with.distances)
#   # remove last edge with Inf weight
#   edgelist.with.distances <- edgelist.with.distances[1:(nrow(edgelist.with.distances) - 1), ]
#   if (table.lengths[1] != FALSE) {
#     inds.in.n <- (offset + 1):(offset + table.lengths[1])
#     edgelist.with.distances <- RemoveWithinNEdges(edgelist.with.distances,
#                                                   inds.in.n = inds.in.n)
#   }
#   num.edges <- length(edgelist.with.distances[, 1])
#   val <- max(floor(num.edges * per / 100), 1)
#   trim.edgelist.with.distances <- edgelist.with.distances[(1:val), ]
#   # calculate "density" for each cluster and normalize
#   if (val == 1) {
#     trim.edgelist.with.distances <- t(trim.edgelist.with.distances)
#   }
#   densities.no.zeros <- table(trim.edgelist.with.distances[, 1:2])
#   # add in zeros for clusters with no edges (table function leaves these out)
#   densities <- rep(0, numcluster)
#   if (table.lengths[1] == FALSE) {
#     names(densities) <- (1:numcluster)
#   } else {
#     names(densities) <- (offset + 1):(offset + table.lengths[1] + table.lengths[2])
#   }
#   densities[names(densities.no.zeros)] <- densities.no.zeros
#   normalized.densities <- round(densities / max(densities) * (max - min) + min)
#   return(list(normalized.densities = normalized.densities,
#               edgelist.with.distances = edgelist.with.distances))
# }

## kNN density function ====
KnnDensity <- function(k, min, max, n, nn.ids.df, nn.dists.df,
                       numcluster = numcluster,
                       table.lengths, offset, table.breaks) {
  all.densities.list <- list()
  for(i in 1:nrow(nn.ids.df)) {
    compile.dists <- c()
    j <- 1
    # while (j <= k) {
    #   #print(paste("j = ", j))
    #   dist.j <- abs(nn.dists.df[i,j])
    #   #print(paste("dist.j = ", dist.j))
    #   compile.dists[[j]] <- dist.j
    #   #print(paste("compile.dists[[j]] = ", compile.dists[[j]]))
    #   j <- j + 1
    # }
    longest.dist <- compile.dists[[k]]
    #density.by.Xshift <- (1/(n*(longest.dist^d))) * (sum(c(1:k)^d)/sum(compile.dists))^d
    #density.by.knn <- sum(compile.dists)
    density.by.knn <- sum(nn.dists.df[i,1:k])#would distance ever be negative? maybe for one of the other dist measures? ====
    all.densities.list[[paste(i,".nn", sep='')]] <- density.by.knn
  }
  densities.df <- rbind(matrix(unlist(all.densities.list), byrow=T))
  global.densities.df <<- densities.df
  normalized.densities <- round(densities.df / max(densities.df) * (max - min) + min)
  global.norm.dens.noNames <<- normalized.densities
  if (table.lengths[1] == FALSE) {
    names(normalized.densities) <- (1:numcluster)
  } else {
    #names(normalized.densities) <- (offset + 1):table.breaks[n + 2]
    names(normalized.densities) <- (offset + 1):(offset + sum(table.lengths))
  }
  global.norm.dens <<- normalized.densities
  return(normalized.densities)
}

DrawNormalizedEdges <- function(output.graph, nn.ids.df, nn.dists.df,
                                normalized.densities, n, offset = FALSE) {
  final.edgelist.with.distances <- c()
  for (i in names(normalized.densities)) {
    ##old
    # matches.in.order <- order(cluster.distances.matrix[, i])
    # tmp.edgelist <- cbind(as.numeric(i), as.numeric(rownames(cluster.distances.matrix)[matches.in.order]),
    #                       sort(cluster.distances.matrix[, i]))[1:normalized.densities[i], ]
    # final.edgelist.with.distances <- rbind(final.edgelist.with.distances, tmp.edgelist)

    ##new not relying on cluster dist matrix
    tmp.edgelist <- cbind(as.numeric(i), as.numeric(nn.ids.df[i,]),as.numeric(nn.dists.df[i,]))[1:normalized.densities[i], ]
    final.edgelist.with.distances <- rbind(final.edgelist.with.distances, tmp.edgelist)
  }
  colnames(final.edgelist.with.distances) <- c("row.inds", "col.inds", "values")
  final.edgelist.with.distances <- RemoveDuplicateValues(final.edgelist.with.distances)
  # note: weight here is really distance, will be converted to weight after the graph is completed
  # try removing edges with 0 distance (Inf weight)
  # if cells are identical in all measurements, distance is 0
  # remove these edges
  vertices.edges <- as.vector(t(as.matrix(final.edgelist.with.distances[, 1:2])))
  output.graph <- add.edges(output.graph, edges = vertices.edges,
                            weight = final.edgelist.with.distances[, 3],
                            label = "EXTRA", sequence_assignment = n)
  return(list(output.graph = output.graph,
              final.edgelist.with.distances = final.edgelist.with.distances))
}
##TO DO: graph_from_data_frame to keep vertex names
ConnectSubgraphs <- function(output.graph, edge.list, offset,
                             table.breaks, n, distance.metric, clusters) {

  global.edge.list[[n]] <<- edge.list
  global.offset[[n]] <<- offset
  global.table.breaks[[n]] <<- table.breaks
  global.n[[n]] <<- n
  global.distance.metric[[n]] <<- distance.metric
  global.clusters[[n]] <<- clusters

  ##make df with vertex names, make graph from edge.list df
  v.df <- data.frame(names = (offset + 1):table.breaks[n + 2])
  time.prox.graph <- graph.data.frame(edge.list[,1:2], vertices = v.df)
  E(time.prox.graph)$weight <- edge.list[,3]
  global.output.graph.full[[n]] <<- time.prox.graph
  global.is.conn.pre[[n]] <<- is.connected(time.prox.graph)

  ##Identify subgraphs of density-based nearest neighbor graph
  subgraphs.ls <- decompose.graph(time.prox.graph)
  global.subgraphs.ls <<- subgraphs.ls

  ##change row names of cluster medians df to have correct indexing
  rownames(clusters) <- c((offset + 1):table.breaks[n + 2])

  ##Compute inter-subgraph edges
  i <- 1
  j <- 1
  to.add.df <- data.frame("v1"='', "v2"='', "dist"='')
  for (p in 1:length(subgraphs.ls)) {
    #print(paste('\n', 'in loop1', '\n'))
    si.graph <- subgraphs.ls[[p]]
    ##Prepare df to hold 1 nn from this subraph to other subgraphs
    inter.sub.nns.df <- data.frame("v1" = '',
                                    "v2" = '',
                                    "dist" = '')
    #print(length((p+1):length(subgraphs.ls)))
    #print(length(subgraphs.ls))
    for (m in 1:length(subgraphs.ls)){
      #print(paste('\n', 'in loop2', '\n'))
      #print(m)
      if (m > length(subgraphs.ls)){
        break
      }
      if (m == p){
        next
      }
      sj.graph <- subgraphs.ls[[m]]
      if (distance.metric == 'manhattan') {
        inter.sub.nns <- RANN.L1::nn2(data = clusters[V(sj.graph)$name,],
                                      query=clusters[V(si.graph)$name,],
                                      k=1, searchtype="priority", eps=0.1)
        #cat(paste('\n',"ran rann",'\n'))
      } else if (distance.metric == 'euclidean') {
        inter.sub.nns <- RANN::nn2(data = clusters[V(sj.graph)$name,],
                                    query=clusters[V(si.graph)$name,],
                                   k=1, searchtype="priority", eps=0.1)
        #cat(paste('\n',"ran rann",'\n'))
      }
      ##add nn output to temp df, then bind to full list
      temp.inter.sub.nns.df <- data.frame("v1" = 1:length(V(si.graph)))
      temp.inter.sub.nns.df <- cbind(temp.inter.sub.nns.df, as.data.frame(inter.sub.nns$nn.idx), as.data.frame(inter.sub.nns$nn.dists))
      colnames(temp.inter.sub.nns.df) <- c("v1", "v2", "dist")
      inter.sub.nns.df <- rbind(inter.sub.nns.df, temp.inter.sub.nns.df)
    }
    ##put list of nearest neighbors in order with shortest dist at top
    inter.sub.nns.df <- inter.sub.nns.df[ do.call(order, inter.sub.nns.df), ]
    to.add.df <- rbind(to.add.df, inter.sub.nns.df[1,])

  }
  ##try adding to edgelist then just making a new graph???? ----------------------------------------------------------------
  time.prox.graph <- add.edges(time.prox.graph, c(to.add.df[,1], to.add.df[,2]),
                            weight = (to.add.df[,3]))
  global.is.conn.pre[[n]] <<- is.connected(time.prox.graph)
  ##Add edge in first row (smallest dist edge between subgraph and any other subgraph) to graph
  ###Thinking there is prob something tricky here about the naming of the vertices...
  output.graph <- add.edges(output.graph, c(to.add.df[,1], to.add.df[,2]),
                            weight = (to.add.df[,3]), label = "proxyMST",
                            sequence_assignment = n)
  global.output.graph.iteration[[n]] <<- output.graph
  return(output.graph)
}


CheckMSTEdges <- function(output.graph, cluster.distances.matrix,
                          table.lengths, offset, n) {
  cluster.distances.matrix <- as.matrix(cluster.distances.matrix)
  # check that any distances between clusters are 0
  # if they are, set them to infinity
  # if two cells are identical in all measurements, they should
  # either be clustered together beforehand or excluded (keep one)
  # to avoid Inf edge weights, set their distance to Inf

  # For each timepoint:
  #   1. make a complete weighted edgelist (this is standard for building
  #      a regular MST)
  #   - This is an edgelist (col 1: vertex 1, col 2: vertex 2,
  #     col 3: distance/weight) for timepoint n+1 vertices to each other
  adjacency.n1n1 <- cluster.distances.matrix[(table.lengths[1] + 1):sum(table.lengths),
                                             (table.lengths[1] + 1):sum(table.lengths)]
  # adjacency.n1n1 is the cluster distances matrix of cells from
  # time n+1 to each other (if second 100 are from n+1, it's a 100 x 100
  # matrix of their distances to each other)
  el.n1n1.with.dist <- cbind((as.vector(row(adjacency.n1n1)) + table.lengths[1] + offset),
                             (as.vector(col(adjacency.n1n1)) + table.lengths[1] + offset),
                             as.vector(adjacency.n1n1))
  colnames(el.n1n1.with.dist) <- c("Vertex 1", "Vertex 2", "Distance")
  # remove permutations (order does not matter)
  # AKA take out edges that are node 1 --- node 3
  # if you already have node 3 --- node 1 edge
  el.n1n1.with.dist.sort <- t(apply(el.n1n1.with.dist, 1, sort))
  el.n1n1.with.dist <- el.n1n1.with.dist[!duplicated(el.n1n1.with.dist.sort), ]
  # remove edges with zero distance/Inf weight
  remove.ind <- which(el.n1n1.with.dist[, "Distance"] == 0)
  if (length(remove.ind) != 0) {
    el.n1n1.with.dist <- el.n1n1.with.dist[-remove.ind, ]
  }
  # el.n1n1.with.dist is the cluster distances matrix of cells from
  # time n+1 to each other in an edgelist format with three columns
  # first being the rows of the distance matrix, second being the
  # columns of the distance matrix, and third being the values
  # of the distance between them (node 1 --- node 2 --- distance 10
  # would be an example of a row of this matrix)

  #   2. for each point in timepoint n+1, create a 3-column matrix with a
  #      row for each cluster of timepoint n+1
  #   column 1: timepoint n+1 vertex index
  #   column 2: vertex in timepoint n that is closest to it
  #   column 3: distance/weight for this edge
  #   - Add to the edgelist from Step 1 i edges (where i is
  #     the number of vertices in timepoint n+1), one for each vertex
  #     connecting them to the closest vertex in timepoint n
  adjacency.n1n <- cluster.distances.matrix[(table.lengths[1] + 1):sum(table.lengths),
                                            1:table.lengths[1]]
  ##### MODIFIED #####
  el.n1n.with.dist <- cbind((as.vector(row(adjacency.n1n)) + table.lengths[1] + offset),
                            (as.vector(col(adjacency.n1n)) + offset),
                            as.vector(adjacency.n1n))
  ##### MODIFIED #####

  ##### ORIGINAL #####
  # el.n1n.with.dist <- cbind((as.vector(row(adjacency.n1n1)) + table.lengths[1] + offset),
  #                           (as.vector(col(adjacency.n1n1)) + offset),
  #                           as.vector(adjacency.n1n))
  ##### ORIGINAL #####

  # rows are n+1 so add offset
  colnames(el.n1n.with.dist) <- c("Vertex 1", "Vertex 2", "Distance")
  # convert row and col to match correct indexing
  el.n1n.with.dist.closest <- matrix(NA, ncol = 3, nrow = table.lengths[2])
  colnames(el.n1n.with.dist.closest) <- c("Vertex 1", "Vertex 2", "Distance")
  x <- 1
  # select only one edge for each vertex in timepoint n+1
  # iterate through edges that have vertex i in timepoint n+1
  # find the smallest distance and keep only that edge

  for (i in ((table.lengths[1] + offset + 1):(sum(table.lengths) + offset))) {
    i.ind <- which(el.n1n.with.dist[, 1] == i)
    i.edges <- el.n1n.with.dist[i.ind, ]
    closest.ind <- which.min(i.edges[, 3])
    closest.edge.i <- i.edges[closest.ind, ]
    el.n1n.with.dist.closest[x, ] <- closest.edge.i
    x <- x + 1
  }
  #   3. add these new n edges to the edgelist, with all n points
  #      given the same unique name/index "tmp.same.name" to distinguish them
  #   - "new n edges" means edges between n+1 and n whereas
  #     step 1 is only within n, give the column 2 "vertex" as
  #     some unique name instead of the vertex it actually is,
  #     so that it treats all those points as one entity/vertex
  original.n1n.closest <- el.n1n.with.dist.closest

  save.name <- el.n1n.with.dist.closest[, 2]
  tmp.same.name <- (sum(table.lengths) + offset + 1)
  el.n1n.with.dist.closest[, 2] <- rep(tmp.same.name,
                                       times = nrow(el.n1n.with.dist.closest))
  full.edgelist <- rbind(el.n1n1.with.dist, el.n1n.with.dist.closest)
  # add one Inf distance edge of timepoint n entity "tmp.same.name" to itself
  full.edgelist <- rbind(full.edgelist, c(tmp.same.name,
                                          tmp.same.name,
                                          Inf))

  # 4. use MST function to find MST from this complete + n edgelist
  dummy.graph <- graph.empty(n = length(unique(full.edgelist[, 1])), directed = FALSE)
  # dummy.graph is an empty graph with enough nodes for
  # all cells from time n + 1 and treat time n cells as one node
  dummy.graph <- add.edges(dummy.graph,
                           edges = (as.vector(t(full.edgelist[, 1:2])) - table.lengths[1] - offset), #why??????????
                           weight = full.edgelist[, 3])
  # has all timepoint n + 1 distances to each other, plus an added
  # row or column that is the one timepoint n entity
  # with its shortest distances to the timepoint n + 1 cells
  mst.dummy.graph <- minimum.spanning.tree(dummy.graph)

  # 5. go through MST edgelist and replace all indices "X" with the
  #    actual n index from column 2 of the 3-column matrix
  #   - take the MST we just calculated, and going back to draw
  #     the real edges to the t-1 timepoint. For MST building,
  #     all of n was just considered (and recorded) as a single
  #     vertex "X", so now we're putting back in the real n vertices.

  # find edge in MST graph (node 1 = node i + 1 or all cells
  # from timepoint n) and based on its partner and its distance
  # match it to an edge from full.edgelist and find out what the
  # one cell is, to generate
  mst.edgelist.fixed <- get.edgelist(mst.dummy.graph)
  mst.edge.weights <- get.edge.attribute(mst.dummy.graph,
                                         "weight", index = E(mst.dummy.graph))
  mst.edgelist.fixed <- cbind(mst.edgelist.fixed, mst.edge.weights)
  colnames(mst.edgelist.fixed) <- c("Vertex 1", "Vertex 2", "Distance")
  # fix vertex IDs back based on offset/table lengths
  mst.edgelist.fixed[, 1:2] <- mst.edgelist.fixed[, 1:2] + table.lengths[1] + offset

  # find what is the vertex from timepoint n, based on distance
  fix.needed.ind <- which(mst.edgelist.fixed[, 1:2] == tmp.same.name, arr.ind = TRUE)
  correct.names <- c("row", "col")
  # if ((colnames(fix.needed.ind) == c("row", "col")) & (nrow(fix.needed.ind) > 1)) {
  if (all.equal(colnames(fix.needed.ind), correct.names) & (nrow(fix.needed.ind) > 1)) {
    # print("more than one edge needs to be renamed")
    temp.el <- mst.edgelist.fixed[fix.needed.ind[, 1], ]
  } else if (nrow(fix.needed.ind) == 1) {
    warning("fix.needed.ind is only 1 row!")
    temp.el <- mst.edgelist.fixed[fix.needed.ind[, 1], ]
    temp.el <- t(as.matrix(temp.el))
  } else {
    stop("fix.needed.ind is only <1 row or different colnames!")
  }
  for (i in 1:nrow(temp.el)) {
    each.edge <- temp.el[i, ]
    each.edge <- as.matrix(each.edge)
    each.edge <- t(each.edge)
    dist.match <- which((original.n1n.closest[, "Distance"] == each.edge[, "Distance"]), arr.ind = TRUE)
    if ((length(dist.match) > 1) & (is.vector(dist.match)) & is.null(names(dist.match))) {
      save.match <- dist.match
      dist.match <- dist.match[1]
      fixed.vertices <- save.name[save.match]
    } else {
      fixed.vertices <- FALSE
    }
    fixed.vertex <- save.name[dist.match]
    # Vertex 2 column will always be the cells from timepoint n
    # that need to be correctd
    mst.edgelist.fixed[fix.needed.ind[i, 1], "Vertex 2"] <- fixed.vertex
  }
  # 6. add edges to actual FLOW-MAP graph
  # if they exist, label them as MST
  # if not, add them and label them as MST

  for (i in 1:nrow(mst.edgelist.fixed)) {
    v1 <- mst.edgelist.fixed[i, "Vertex 1"]
    v2 <- mst.edgelist.fixed[i, "Vertex 2"]
    if (are.connected(output.graph, v1, v2)) {
      E(output.graph, P = c(v1, v2))$label <- "MST"
      if ((v2 == fixed.vertex) & (fixed.vertices != FALSE)) {
        for (i in 2:length(fixed.vertices)) {
          v2 <- fixed.vertices[i]
          if (are.connected(output.graph, v1, v2)) {
            E(output.graph, P = c(v1, v2))$label <- "MST"
          } else {
            output.graph <- add.edges(output.graph, c(v1, v2),
                                      weight = mst.edgelist.fixed[i, "Distance"], label = "MST",
                                      sequence_assignment = n)
          }
        }
      }
    } else {
      output.graph <- add.edges(output.graph, c(v1, v2),
                                weight = mst.edgelist.fixed[i, "Distance"], label = "MST",
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

BuildFirstFLOWMAP <- function(FLOWMAP.clusters, k, min, max, distance.metric,
                              clustering.var) {
  n <- 0
  output.graph <- InitializeGraph(FLOWMAP.clusters = FLOWMAP.clusters)
  # This section creates a flowmap for the first time point
  cat("Building first FLOWMAP:\n")
  table.lengths <- FLOWMAP.clusters$table.lengths
  ## Old full dist matrix calculation
  # # get distance matrix from clusters
  # clusters.old <<- FLOWMAP.clusters$full.clusters[1:table.lengths[1], ]
  # clusters.old <<- subset(clusters.old, select = clustering.var)
  # cluster.distances <- dist(clusters.old, method = distance.metric, diag = TRUE, upper = TRUE)
  # cluster.distances.matrix <- as.matrix(cluster.distances)
  # cluster.distances.matrix <- as.data.frame(cluster.distances.matrix)
  # # set i-i distances to Inf instead of 0 so they aren't the closest neighbors
  # for (i in 1:ncol(cluster.distances.matrix)) {
  #   cluster.distances.matrix[i, i] <- Inf
  # }
  # numcluster <- nrow(clusters)
  # normalized.results <- FindNormalized(cluster.distances.matrix = cluster.distances.matrix,
  #                                      per = per, min = min, max = max, numcluster = numcluster,
  #                                      table.lengths = FALSE)

  ##Replace with approximate nearest neighbors
  clusters <- FLOWMAP.clusters$full.clusters[1:table.lengths[1], ]
  clusters.new <<- subset(clusters, select = clustering.var)
  if (distance.metric == 'manhattan') {
    nns <<- RANN.L1::nn2(data=clusters.new, k=(max + 1), searchtype="priority", eps=0.1)
  } else if (distance.metric == 'euclidean') {
    nns <<- RANN::nn2(data=clusters.new, k=(max + 1), searchtype="priority", eps=0.1)
  }
  temp_nnids.df <- as.data.frame(nns$nn.idx)
  temp_nndists.df <- as.data.frame(nns$nn.dists)
  nn.ids.df <<- temp_nnids.df[,2:length(temp_nnids.df)]
  nn.dists.df <<- temp_nndists.df[,2:length(temp_nndists.df)]
  numcluster <- nrow(clusters.new)
  normalized.densities <- KnnDensity(k=k, min, max, n=n,
                                   nn.ids.df = nn.ids.df,
                                   nn.dists.df = nn.dists.df,
                                   numcluster = numcluster,
                                   table.lengths = FALSE)
  #normalized.densities <- normalized.results$normalized.densities
  #edgelist.with.distances <- normalized.results$edgelist.with.distances <--- don't need this for anything...
  # build new edgelist with N edges for each cluster based on normalized density

  results <- DrawNormalizedEdges(output.graph = output.graph,
                                 nn.ids.df = nn.ids.df,
                                 nn.dists.df = nn.dists.df,
                                 normalized.densities = normalized.densities,
                                 n = n, offset = FALSE)
  output.graph <- results$output.graph

  edgelist.save <- results$final.edgelist.with.distances
  ##Temp version with no MST check 10.01.18
  # # now add all MST edges that are not yet included in the graph, and annotate all as "MST"
  # adjacency.graph <- graph.adjacency(as.matrix(cluster.distances.matrix),
  #                                    mode = "undirected", weighted = "weight", diag = FALSE)
  # mst.graph <- minimum.spanning.tree(adjacency.graph, algorithm = "prim")
  # # for each edge of the mst, if it exists in the graph then label MST, if it doesn't exist then add it
  # mst.graph.edgelist <- cbind(get.edgelist(mst.graph), E(mst.graph)$weight)
  # class(mst.graph.edgelist) <- "numeric"
  #
  # # for each edge of the mst, if it exists in the graph then label MST, if it doesn't exist then add it
  # for (i in 1:nrow(mst.graph.edgelist)) {
  #   if (are.connected(output.graph, mst.graph.edgelist[i, 1], mst.graph.edgelist[i, 2])) {
  #     E(output.graph, P = c(mst.graph.edgelist[i, 1],
  #                           mst.graph.edgelist[i, 2]))$label <- "MST"
  #   } else {
  #     output.graph <- add.edges(output.graph,
  #                               as.numeric(mst.graph.edgelist[i, 1:2]),
  #                               weight = mst.graph.edgelist[i, 3],
  #                               label = "MST", sequence_assignment = n)
  #   }
  # }
  return(list(output.graph = output.graph,
              edgelist.save = edgelist.save))
}

BuildFLOWMAP <- function(FLOWMAP.clusters, k, min, max,
                         distance.metric, clustering.var) {
  edgelist.save <- list()
  first.results <- BuildFirstFLOWMAP(FLOWMAP.clusters = FLOWMAP.clusters,
                                     k = k, min = min, max = max,
                                     distance.metric = distance.metric,
                                     clustering.var = clustering.var)
  output.graph <- first.results$output.graph
  edgelist.save[["first"]] <- first.results$edgelist.save
  table.breaks <<- c(0, FLOWMAP.clusters$table.breaks)
  # This section builds the flowmap one timepoint at a time
  for (n in 1:(length(FLOWMAP.clusters$cluster.medians) - 1)) {
    # offset value is used to correctly index the edges at each sequential step
    offset <- table.breaks[n]
    # go through sequential cluster sets, add edges for each a-a and a-a+1 set, and also label sequential MST
    cat("Building FLOWMAP from", n, "to", n + 1, "\n")
    # # # get clusters for time a and a+1
    #  clusters <- rbind(FLOWMAP.clusters$cluster.medians[[n]], FLOWMAP.clusters$cluster.medians[[n + 1]])
    #  clusters <- subset(clusters, select = clustering.var)
    # # numcluster <- nrow(clusters)
    # # make adjacency matrix from clusters
    # cluster.distances <- dist(clusters, method = distance.metric, diag = TRUE, upper = TRUE)
    # cluster.distances.matrix <- as.matrix(cluster.distances)
    # cluster.distances.matrix <- as.data.frame(cluster.distances.matrix)
    # rownames(cluster.distances.matrix) <- (offset + 1):table.breaks[n + 2]
    # colnames(cluster.distances.matrix) <- (offset + 1):table.breaks[n + 2]
    # # set i-i distances to Inf instead of 0 so they aren't the closest neighbors
    # for (i in 1:ncol(cluster.distances.matrix)) {
    #   cluster.distances.matrix[i, i] <- Inf
    # }
    # # This section adds the lowest distance n_n+1 and n+1_n+1 edges to the output graph
    n_n1.table.lengths <- FLOWMAP.clusters$table.lengths[n:(n + 1)]
    n_n1.table.breaks <- table.breaks[n:(n + 1)]
    # normalized.results <- FindNormalized(cluster.distances.matrix = cluster.distances.matrix,
    #                                      per = per, min = min, max = max, numcluster = numcluster,
    #                                      table.lengths = n_n1.table.lengths,
    #                                      table.breaks = n_n1.table.breaks,
    #                                      offset = offset)
    # normalized.densities <- normalized.results$normalized.densities
    #
    # edgelist.with.distances <- normalized.results$edgelist.with.distances

    ##Replace with approximate nearest neighbors
    clusters <- rbind(FLOWMAP.clusters$cluster.medians[[n]], FLOWMAP.clusters$cluster.medians[[n + 1]])
    clusters <- subset(clusters, select = clustering.var)
    if (distance.metric == 'manhattan') {
      nns <- RANN.L1::nn2(data=clusters, k=(max + 1), searchtype="priority", eps=0.1)
    } else if (distance.metric == 'euclidean') {
      nns <- RANN::nn2(data=clusters, k=(max + 1), searchtype="priority", eps=0.1)
    }
    temp_nnids.df <- as.data.frame(nns$nn.idx)
    temp_nndists.df <- as.data.frame(nns$nn.dists)
    nn.ids.df <- temp_nnids.df[,2:length(temp_nnids.df)]
    global.nn.ids.df <<- nn.ids.df
    nn.dists.df <- temp_nndists.df[,2:length(temp_nndists.df)]
    row.names(nn.ids.df) <- (offset + 1):table.breaks[n + 2]
    nn.ids.df[] <- lapply(nn.ids.df, function(x) x+offset)
    row.names(nn.dists.df) <- (offset + 1):table.breaks[n + 2]
    numcluster <- nrow(clusters)
    normalized.densities <- KnnDensity(k=k, min, max, n=n,
                                       nn.ids.df = nn.ids.df,
                                       nn.dists.df = nn.dists.df,
                                       numcluster = numcluster,
                                       table.lengths = n_n1.table.lengths,
                                       offset = offset,
                                       table.breaks = table.breaks)


    # build new edgelist with N edges for each cluster based on normalized density
    results <- DrawNormalizedEdges(output.graph = output.graph,
                                   nn.ids.df = nn.ids.df,
                                   nn.dists.df = nn.dists.df,
                                   normalized.densities = normalized.densities,
                                   n = n, offset = offset)
    output.graph <- results$output.graph
    edgelist.save[[n]] <- results$final.edgelist.with.distances
    ##Temp version with no MST check 10.01.18
    # # This section adds the "MST" for n_n+1 and n+1_n+1 nodes
    # output.graph <- CheckMSTEdges(output.graph = output.graph,
    #                               cluster.distances.matrix = cluster.distances.matrix,
    #                               table.lengths = n_n1.table.lengths, n = n, offset = offset)

    ##Temp version with new ConnectSubgraphs function, only used if graph is not connected
    if (!is_connected(output.graph)) {
      cat("Graph has disconnected components!")
      output.graph <- ConnectSubgraphs(output.graph = output.graph,
                                       edge.list = edgelist.save[[n]],
                                       offset = offset,
                                       table.breaks = table.breaks,
                                       n = n,
                                       distance.metric = distance.metric,
                                       clusters = clusters)
    }
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
  weights <- (1 / distances)
  E(output.graph)$weight <- weights
  output.graph <- AnnotateGraph(output.graph = output.graph,
                                FLOWMAP.clusters = FLOWMAP.clusters)
  global.output.graph <<- output.graph
  return(list(output.graph = output.graph,
              edgelist.save = edgelist.save))
}

AnnotateGraph <- function(output.graph, FLOWMAP.clusters) {
  # This section annotates the graph
  anno.cat <- c()
  anno <- list()
  # iterate through all times and annotate
  for (f in 1:length(FLOWMAP.clusters$cluster.medians)) {
    cat("Annotating graph for file", f, "\n")
    # get medians for all parameters and counts for all clusters
    counts <- FLOWMAP.clusters$cluster.counts[[f]]$Counts
    anno$count <- counts
    total.cell <- sum(counts)
    anno$percent.total <- data.frame(percent.total = c(counts / total.cell))
    anno$medians  <- FLOWMAP.clusters$cluster.medians[[f]]
    # add time information column
    time.matrix <- matrix(f, nrow = length(anno$count))
    colnames(time.matrix) <- c("Timepoint")
    anno$medians <- cbind(anno$medians, time.matrix)
    # add median and percent values
    for (col in c("medians", "percent.total")) {
      anno.cat[[col]] <- rbind(anno.cat[[col]], anno[[col]])
    }
  }
  # combine anno.cat matrices
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

AnnotateSpecialGraph <- function(output.graph, FLOWMAP.clusters,
                                 label.key.special) {
  # This section annotates the graph
  anno.cat <- c()
  anno <- list()
  # iterate through all times and annotate
  for (f in 1:length(FLOWMAP.clusters$cluster.medians)) {
    cat("Annotating graph for file", f, "\n")
    # get medians for all parameters and counts for all clusters
    counts <- FLOWMAP.clusters$cluster.counts[[f]]$Counts
    anno$count <- counts
    total.cell <- sum(counts)
    anno$percent.total <- data.frame(percent.total = c(counts / total.cell))
    anno$medians  <- FLOWMAP.clusters$cluster.medians[[f]]
    # add time information column
    time.matrix <- matrix(f, nrow = length(anno$count))
    colnames(time.matrix) <- c("Timepoint")
    anno$medians <- cbind(anno$medians, time.matrix)
    # add median and percent values
    for (col in c("medians", "percent.total")) {
      anno.cat[[col]] <- rbind(anno.cat[[col]], anno[[col]])
    }
  }
  # combine anno.cat matrices
  output.anno <- cbind(anno.cat[[1]], anno.cat[[2]])
  print("label.key.special")
  print(label.key.special)
  output.anno <- ConvertCharacterLabelSpecial(output.anno, label.key.special)
  for (c in colnames(output.anno)) {
    if (c == "Condition") {
      output.graph <- set.vertex.attribute(output.graph, c,
                                           index = as.numeric(1:dim(output.anno)[1]),
                                           value = as.character(output.anno[, c]))
    } else {
      output.graph <- set.vertex.attribute(output.graph, c,
                                           index = as.numeric(1:dim(output.anno)[1]),
                                           value = output.anno[, c])
    }
  }
  # add name attribute
  V(output.graph)$name <- 1:length(V(output.graph))
  return(output.graph)
}
