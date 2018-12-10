#############################################################################
##FUNCTIONS FOR RADIUS-BASED DENSITY ====
#############################################################################

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

FindNormalized <- function(cluster.distances.matrix, per, min,
                           max, numcluster, table.lengths,
                           table.breaks, offset) {
  # make fully connected graph from adjacency list
  # note: "weight" here is really distance, calling it weight for the mst function later needs this
  arr.inds <- ConvertIndex(inds = order(cluster.distances.matrix), orig.data.frame = cluster.distances.matrix, keep.names = TRUE)
  edgelist.with.distances <- MergeValues(arr.inds = arr.inds, orig.data.frame = cluster.distances.matrix)
  edgelist.with.distances <- RemoveDuplicateValues(edgelist.with.distances)
  # remove last edge with Inf weight
  edgelist.with.distances <- edgelist.with.distances[1:(nrow(edgelist.with.distances) - 1), ]
  if (table.lengths[1] != FALSE) {
    inds.in.n <- (offset + 1):(offset + table.lengths[1])
    edgelist.with.distances <- RemoveWithinNEdges(edgelist.with.distances,
                                                  inds.in.n = inds.in.n)
  }
  num.edges <- length(edgelist.with.distances[, 1])
  val <- max(floor(num.edges * per / 100), 1)
  trim.edgelist.with.distances <- edgelist.with.distances[(1:val), ]
  # calculate "density" for each cluster and normalize
  if (val == 1) {
    trim.edgelist.with.distances <- t(trim.edgelist.with.distances)
  }
  densities.no.zeros <- table(trim.edgelist.with.distances[, 1:2])
  # add in zeros for clusters with no edges (table function leaves these out)
  densities <- rep(0, numcluster)
  if (table.lengths[1] == FALSE) {
    names(densities) <- (1:numcluster)
  } else {
    names(densities) <- (offset + 1):(offset + table.lengths[1] + table.lengths[2])
  }
  densities[names(densities.no.zeros)] <- densities.no.zeros
  normalized.densities <- round(densities / max(densities) * (max - min) + min)
  return(list(normalized.densities = normalized.densities,
              edgelist.with.distances = edgelist.with.distances))
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
                           edges = (as.vector(t(full.edgelist[, 1:2])) - table.lengths[1] - offset),
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

BuildFirstFLOWMAP <- function(FLOWMAP.clusters, per, min, max, distance.metric,
                              clustering.var) {
  n <- 0
  output.graph <- InitializeGraph(FLOWMAP.clusters = FLOWMAP.clusters)
  # This section creates a flowmap for the first time point
  cat("Building first FLOWMAP:\n")
  table.lengths <- FLOWMAP.clusters$table.lengths
  # get distance matrix from clusters
  clusters <- FLOWMAP.clusters$full.clusters[1:table.lengths[1], ]
  clusters <- subset(clusters, select = clustering.var)
  cluster.distances <- dist(clusters, method = distance.metric, diag = TRUE, upper = TRUE)
  cluster.distances.matrix <- as.matrix(cluster.distances)
  cluster.distances.matrix <- as.data.frame(cluster.distances.matrix)
  # set i-i distances to Inf instead of 0 so they aren't the closest neighbors
  for (i in 1:ncol(cluster.distances.matrix)) {
    cluster.distances.matrix[i, i] <- Inf
  }
  numcluster <- nrow(clusters)
  normalized.results <- FindNormalized(cluster.distances.matrix = cluster.distances.matrix,
                                       per = per, min = min, max = max, numcluster = numcluster,
                                       table.lengths = FALSE)
  normalized.densities <- normalized.results$normalized.densities
  edgelist.with.distances <- normalized.results$edgelist.with.distances
  # build new edgelist with N edges for each cluster based on normalized density

  results <- DrawNormalizedEdges(output.graph = output.graph,
                                 cluster.distances.matrix = cluster.distances.matrix,
                                 normalized.densities = normalized.densities,
                                 n = n, offset = FALSE)
  output.graph <- results$output.graph

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
  return(list(output.graph = output.graph,
              edgelist.save = edgelist.save))
}

BuildFLOWMAP <- function(FLOWMAP.clusters, per, min, max,
                         distance.metric, clustering.var) {
  edgelist.save <- list()
  first.results <- BuildFirstFLOWMAP(FLOWMAP.clusters = FLOWMAP.clusters,
                                     per = per, min = min, max = max,
                                     distance.metric = distance.metric,
                                     clustering.var = clustering.var)
  output.graph <- first.results$output.graph
  edgelist.save[["first"]] <- first.results$edgelist.save
  table.breaks <- c(0, FLOWMAP.clusters$table.breaks)
  # This section builds the flowmap one timepoint at a time
  for (n in 1:(length(FLOWMAP.clusters$cluster.medians) - 1)) {
    # offset value is used to correctly index the edges at each sequential step
    offset <- table.breaks[n]
    # go through sequential cluster sets, add edges for each a-a and a-a+1 set, and also label sequential MST
    cat("Building FLOWMAP from", n, "to", n + 1, "\n")
    # get clusters for time a and a+1
    clusters <- rbind(FLOWMAP.clusters$cluster.medians[[n]], FLOWMAP.clusters$cluster.medians[[n + 1]])
    clusters <- subset(clusters, select = clustering.var)
    numcluster <- nrow(clusters)
    # make adjacency matrix from clusters
    cluster.distances <- dist(clusters, method = distance.metric, diag = TRUE, upper = TRUE)
    cluster.distances.matrix <- as.matrix(cluster.distances)
    cluster.distances.matrix <- as.data.frame(cluster.distances.matrix)
    rownames(cluster.distances.matrix) <- (offset + 1):table.breaks[n + 2]
    colnames(cluster.distances.matrix) <- (offset + 1):table.breaks[n + 2]
    # set i-i distances to Inf instead of 0 so they aren't the closest neighbors
    for (i in 1:ncol(cluster.distances.matrix)) {
      cluster.distances.matrix[i, i] <- Inf
    }
    # This section adds the lowest distance n_n+1 and n+1_n+1 edges to the output graph
    n_n1.table.lengths <- FLOWMAP.clusters$table.lengths[n:(n + 1)]
    n_n1.table.breaks <- table.breaks[n:(n + 1)]
    normalized.results <- FindNormalized(cluster.distances.matrix = cluster.distances.matrix,
                                         per = per, min = min, max = max, numcluster = numcluster,
                                         table.lengths = n_n1.table.lengths,
                                         table.breaks = n_n1.table.breaks,
                                         offset = offset)
    normalized.densities <- normalized.results$normalized.densities

    edgelist.with.distances <- normalized.results$edgelist.with.distances
    # build new edgelist with N edges for each cluster based on normalized density
    results <- DrawNormalizedEdges(output.graph = output.graph,
                                   cluster.distances.matrix = cluster.distances.matrix,
                                   normalized.densities = normalized.densities,
                                   n = n, offset = offset)
    output.graph <- results$output.graph
    edgelist.save[[n]] <- results$final.edgelist.with.distances
    # This section adds the "MST" for n_n+1 and n+1_n+1 nodes
    output.graph <- CheckMSTEdges(output.graph = output.graph,
                                  cluster.distances.matrix = cluster.distances.matrix,
                                  table.lengths = n_n1.table.lengths, n = n, offset = offset)
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
  return(list(output.graph = output.graph,
              edgelist.save = edgelist.save))
}

# #############################################################################
# ##FUNCTIONS FOR KNN-BASED DENSITY ====
# #############################################################################

## kNN density function ====
KnnDensity <- function(k, min, max, n, nn.ids.df, nn.dists.df,
                       numcluster = numcluster,
                       table.lengths, offset) {
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
  if (offset == FALSE) {
    names(normalized.densities) <- (1:numcluster)
  } else {
    names(normalized.densities) <- (offset + 1):(offset + sum(table.lengths))
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
  rownames(clusters) <- c((offset + 1):table.breaks[n + 2])
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
            inter.sub.nns <- RANN.L1::nn2(data = clusters[V(sj.graph)$name,],
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
                                    as.character((offset + 1):table.breaks[n + 2]))
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
  time.prox.graph <- set.vertex.attribute(time.prox.graph,'name',index=V(time.prox.graph),as.character((offset + 1):table.breaks[n + 2]))
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
      rownames(clusters) <- c((offset + 1):table.breaks[n + 2])

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
            inter.sub.nns <- RANN.L1::nn2(data = clusters[V(sj.graph)$name,],
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
                                                  as.character((offset + 1):table.breaks[n + 2]))
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
                                offset = FALSE) {
  final.edgelist.with.distances <- c()
  for (i in names(normalized.densities)) {
    ##new not relying on cluster dist matrix
    tmp.edgelist <- cbind(as.numeric(i), as.numeric(nn.ids.df[i,]),as.numeric(nn.dists.df[i,]))[1:normalized.densities[i], ]
    final.edgelist.with.distances <- rbind(final.edgelist.with.distances, tmp.edgelist)
  }
  colnames(final.edgelist.with.distances) <- c("row.inds", "col.inds", "values")
  #global.draw.norm.edgelist <<- final.edgelist.with.distances
  ##make df with vertex names, make graph from edge.list df
  if (offset == FALSE) {
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

BuildFirstFLOWMAPkNN <- function(FLOWMAP.clusters, k, min, max, distance.metric,
                                 clustering.var) {
  #global.flowmap.clusters <<- FLOWMAP.clusters
  n <- 0
  # This section creates a flowmap for the first time point
  cat("Building first FLOWMAP:\n")
  ##Read in info from FLOWMAP.clusters
  table.lengths <- FLOWMAP.clusters$table.lengths
  table.breaks <- c(0, FLOWMAP.clusters$table.breaks)
  clusters <- FLOWMAP.clusters$full.clusters[1:table.lengths[1], ]
  clusters <- subset(clusters, select = clustering.var)
  ##Calculate k nearest neighboors for clusters/cells
  if (distance.metric == 'manhattan') {
    nns <- RANN.L1::nn2(data=clusters, k=k+1, searchtype="priority", eps=0.1)
  } else if (distance.metric == 'euclidean') {
    nns <- RANN::nn2(data=clusters, k=k+1, searchtype="priority", eps=0.1)
  }
  temp_nnids.df <- as.data.frame(nns$nn.idx)
  temp_nndists.df <- as.data.frame(nns$nn.dists)
  nn.ids.df <- temp_nnids.df[,2:length(temp_nnids.df)]
  nn.dists.df <- temp_nndists.df[,2:length(temp_nndists.df)]
  numcluster <- nrow(clusters)
  ##Calculate density based on kNN
  normalized.densities <- KnnDensity(k=k, min, max, n=n,
                                     nn.ids.df = nn.ids.df,
                                     nn.dists.df = nn.dists.df,
                                     numcluster = numcluster,
                                     table.lengths = table.lengths,
                                     offset = FALSE)
  ##Build new edgelist with N edges for each cluster based on normalized density
  results <- DrawNormalizedEdgesKnn(output.graph = graph.empty() ,
                                    nn.ids.df = nn.ids.df,
                                    nn.dists.df = nn.dists.df,
                                    normalized.densities = normalized.densities,
                                    n = n, table.breaks = table.breaks, offset = FALSE)
  output.graph <- results$output.graph
  edgelist.save <- results$final.edgelist.with.distances
  offset <- 0
  ##Test whether graph is connected, connect by kNN if not
  if (!is_connected(output.graph)) {
    print("Graph has disconnected components!")
    results <- FirstConnectSubgraphs(output.graph = output.graph,
                                     edge.list = edgelist.save,
                                     offset = offset,
                                     table.breaks = table.breaks,
                                     n = n,
                                     distance.metric = distance.metric,
                                     clusters = clusters)
  }
  output.graph <- results$output.graph
  return(list(output.graph = output.graph,
              edgelist.save = edgelist.save))
}

BuildFLOWMAPkNN <- function(FLOWMAP.clusters, k, min, max,
                            distance.metric, clustering.var) {
  edgelist.save <- list()
  first.results <- BuildFirstFLOWMAPkNN(FLOWMAP.clusters = FLOWMAP.clusters,
                                        k = k, min = min, max = max,
                                        distance.metric = distance.metric,
                                        clustering.var = clustering.var)
  output.graph <- first.results$output.graph
  #global.first.results <<- first.results
  edgelist.save[["first"]] <- first.results$edgelist.save
  table.breaks <- c(0, FLOWMAP.clusters$table.breaks)
  # This section builds the flowmap one timepoint at a time
  for (n in 1:(length(FLOWMAP.clusters$cluster.medians) - 1)) {
    # offset value is used to correctly index the edges at each sequential step
    offset <- table.breaks[n]
    # go through sequential cluster sets, add edges for each a-a and a-a+1 set
    print(paste("Building FLOWMAP from", n, "to", n + 1, "\n", sep = ''))

    n_n1.table.lengths <- FLOWMAP.clusters$table.lengths[n:(n + 1)]
    n_n1.table.breaks <- table.breaks[n:(n + 1)]

    ##get clusters for time a and a+1
    clusters <- rbind(FLOWMAP.clusters$cluster.medians[[n]], FLOWMAP.clusters$cluster.medians[[n + 1]])
    clusters <- subset(clusters, select = clustering.var)
    if (distance.metric == 'manhattan') {
      nns <- RANN.L1::nn2(data=clusters, k=k+1, searchtype="priority", eps=0.1)
    } else if (distance.metric == 'euclidean') {
      nns <- RANN::nn2(data=clusters, k=k+1, searchtype="priority", eps=0.1)
    }
    temp_nnids.df <- as.data.frame(nns$nn.idx)
    temp_nndists.df <- as.data.frame(nns$nn.dists)
    nn.ids.df <- temp_nnids.df[,2:length(temp_nnids.df)]
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
                                       offset = offset)


    # build new edgelist with N edges for each cluster based on normalized density
    results <- DrawNormalizedEdgesKnn(output.graph = output.graph,
                                   nn.ids.df = nn.ids.df,
                                   nn.dists.df = nn.dists.df,
                                   normalized.densities = normalized.densities,
                                   n = n, table.breaks = table.breaks, offset = offset)
    output.graph <- results$output.graph
    #global.post.DrawNorm <<- output.graph
    edgelist.save[[n]] <- results$final.edgelist.with.distances

    ##Temp version with new ConnectSubgraphs function, only used if graph is not connected
    if (!is_connected(output.graph)) {
      print("Graph has disconnected components!")
      results <- ConnectSubgraphs(output.graph = output.graph,
                                  edge.list = edgelist.save[[n]],
                                  offset = offset,
                                  table.breaks = table.breaks,
                                  n = n,
                                  distance.metric = distance.metric,
                                  clusters = clusters)
      output.graph <- results$output.graph
      #global.post.connect <<- output.graph
      edgelist.save[[n]] <- results$edgelist.with.distances
    }
  }
  output.graph.final <- graph.empty()
  vertices.edges <- get.edgelist(output.graph)
  vertices.edges <- as.matrix(data.frame(lapply(data.frame(vertices.edges), function(x) as.numeric(as.character(x)))))
  output.graph.final <- graph_from_edgelist(vertices.edges, directed = FALSE)
  E(output.graph.final)$weight <- unlist(lapply(E(output.graph)$weight, function(x) 1/as.numeric(x)))
  #global.output.graph.final.b <<- output.graph.final
  output.graph.final <- AnnotateGraph(output.graph = output.graph.final,
                                      FLOWMAP.clusters = FLOWMAP.clusters)
  #global.pre.simplify.el <<- get.edgelist(output.graph.final)
  ##Remove duplicates and self-edges
  output.graph.final <- simplify(output.graph.final)
  #global.post.simplify.el <<- get.edgelist(output.graph.final)
  #global.output.graph.final.c <<- output.graph.final
  return(list(output.graph = output.graph.final,
              edgelist.save = edgelist.save))
}

#############################################################################
##FUNCTIONS SHARED FOR RADIUS AND KNN DENSITY ====
#############################################################################

InitializeGraph <- function(FLOWMAP.clusters) {
  # This section initializes the graph
  total.nodes <- length(FLOWMAP.clusters$full.clusters[, 1])
  # create empty graph with the right number of nodes - will fill in the edges later
  initial.graph <- graph.empty(n = total.nodes, directed = FALSE)
  return(initial.graph)
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
