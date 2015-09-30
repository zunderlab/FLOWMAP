rm(list=ls(all=T))
library(igraph)
library(spade)

CLUSTERS_TABLE <- "clusters.table"
DISTANCE_METHOD <- "manhattan"
FCS_FILE_PATTERN     <- "*.fcs.density.fcs.cluster.fcs$"  # File pattern for FCS files with cluster assignment

# set up percentage of edges to take
PERCENT_EDGES <- 1
MAX_EDGES <- 10
MIN_EDGES <- 2

################################################################################
# This section collects the cluster data and initializes the graph
################################################################################

# get clusters output directories
cluster_table_files <- list.files(pattern=CLUSTERS_TABLE,recursive=TRUE,full.names=TRUE)
cluster_table_files <- sort(cluster_table_files)

# make full table of clusters
clusters_full <- c()
table_breaks   <- c() # this will be used later in the annotation section
table_lengths <- c()
for (i in 1:length(cluster_table_files)) {
  cat("Reading cluster data from:",cluster_table_files[i],"\n")
  tmp_table <- read.table(cluster_table_files[i], header=TRUE, check.names=FALSE)
  table_lengths <- append(table_lengths,nrow(tmp_table)) # used for node numbering/index
  clusters_full <- rbind(clusters_full,tmp_table)
  table_breaks <- append(table_breaks, nrow(clusters_full)) # used for node numbering/index
}

# create empty graph with the right number of nodes - will fill in the edges later
total_nodes <- length(clusters_full[,1])
output_graph <- graph.empty(n=total_nodes,directed=FALSE)

sequence_assignments <- c()
for (i in 1:length(cluster_table_files)) {
  sequence_assignments <- c(sequence_assignments,rep(i,table_lengths[i]))
}

E(output_graph)$sequence_assignemt <- sequence_assignments

################################################################################
# This section creates a flowmap for the first time point
################################################################################

cat("building first flowmap\n")
# get distance matrix from clusters
clusters <- clusters_full[1:table_lengths[1],]
cluster_distances <- dist(clusters, method=DISTANCE_METHOD)
cluster_distances_matrix <- as.matrix(cluster_distances)
# set i-i distances to Inf instead of 0 so they aren't the closest neighbors
for (i in 1:ncol(cluster_distances_matrix)) {
  cluster_distances_matrix[i,i] <- Inf
}
# make fully connected graph from adjacency list
# note: "weight" here is really distance, calling it weight for the mst function later needs this
edgelist_with_distances <- cbind(as.vector(row(cluster_distances_matrix)),
        as.vector(col(cluster_distances_matrix)),as.vector(cluster_distances_matrix))
# take only the upper triangle of the matrix
edgelist_with_distances <- edgelist_with_distances[upper.tri(matrix(data=1:length(cluster_distances_matrix),
        nrow=nrow(cluster_distances_matrix),ncol=ncol(cluster_distances_matrix))),]
# convert strings to numeric
class(edgelist_with_distances) <- "numeric"
# sort edges by distance (column 3)
edgelist_with_distances <- edgelist_with_distances[order(edgelist_with_distances[,3]),]
# take only the top X-% edges by distance based on PERCENT_TOTAL
trim_edgelist_with_distances <- edgelist_with_distances[1:floor(length(edgelist_with_distances[,1])*PERCENT_EDGES/100),]
# calculate "density" for each cluster and normalize
densities_no_zeros <- table(trim_edgelist_with_distances[,1:2])
# add in zeros for clusters with no edges (table function leaves these out)
densities <- rep(0,nrow(clusters))
densities[as.numeric(names(densities_no_zeros))] <- densities_no_zeros
normalized_densities <- round(densities/max(densities)*(MAX_EDGES - MIN_EDGES) + MIN_EDGES)

# build new edgelist with N edges for each cluster based on normalized density
cat("building first edgelist\n")
final_edgelist_with_distances <- c()
for (n in 1:length(normalized_densities)) {
  tmp_final_edgelist_with_distances <- cbind(n,order(cluster_distances_matrix[,n]),
                                             sort(cluster_distances_matrix[,n]))[1:normalized_densities[n],]
  final_edgelist_with_distances <- rbind(final_edgelist_with_distances,tmp_final_edgelist_with_distances)
}

# remove duplicate edges from edgelist
final_edgelist_with_distances <- unique(cbind(t(apply(final_edgelist_with_distances[,1:2],1,sort)),final_edgelist_with_distances[,3]))


# add edges to graph
# note: weight here is really distance, will be converted to weight after the graph is completed
# note: use -1 because igraph graphs are zero indexed
output_graph <- add.edges(output_graph,edges=as.vector(t(final_edgelist_with_distances[,1:2])-1),
                          weight=final_edgelist_with_distances[,3],label="EXTRA",sequence_assignment=1)

# now add all MST edges that are not yet included in the graph, and annotate all as "MST"
adjacency_graph <- graph.adjacency(cluster_distances_matrix, mode="undirected", weighted="weight", diag=FALSE)
mst_graph <- minimum.spanning.tree(adjacency_graph, algorithm="prim")
# for each edge of the mst, if it exists in the graph then label MST, if it doesn't exist then add it
mst_graph_edgelist <- cbind(get.edgelist(mst_graph),E(mst_graph)$weight)
class(mst_graph_edgelist) <- "numeric"
# for each edge of the mst, if it exists in the graph then label MST, if it doesn't exist then add it
for (i in 1:nrow(mst_graph_edgelist)) {
  if (are.connected(output_graph,mst_graph_edgelist[i,1]-1,mst_graph_edgelist[i,2]-1)) {
    E(output_graph,P=c(mst_graph_edgelist[i,1]-1,mst_graph_edgelist[i,2]-1))$label <- "MST"
  } else {
    output_graph <- add.edges(output_graph,as.numeric(mst_graph_edgelist[i,1:2]-1),weight=mst_graph_edgelist[i,3],label="MST",sequence_assignment=1)
  }
}

################################################################################
# This section builds the flowmap one timepoint at a time
################################################################################

# offset value is used to correctly index the edges at each sequential step
offset=0

# go through sequential cluster sets, add edges for each a-a and a-a+1 set, and also label sequential MST
for (n in 1:(length(cluster_table_files)-1)) {
  cat(paste("Build FlowMap from",n,"to",n+1,"\n"))
  # get clusters for time a and a+1
  clusters <- rbind(read.table(cluster_table_files[n], header=TRUE, check.names=FALSE),read.table(cluster_table_files[n+1], header=TRUE, check.names=FALSE))
  # make adjacency matrix from clusters
  cluster_distances <- dist(clusters, method=DISTANCE_METHOD)
  cluster_distances_matrix <- as.matrix(cluster_distances)
  # set i-i distances to Inf instead of 0 so they aren't the closest neighbors
  for (i in 1:ncol(cluster_distances_matrix)) {
    cluster_distances_matrix[i,i] <- Inf
  }
  
  ###################################################################################
  # This section adds the lowest distance n_n+1 and n+1_n+1 edges to the output graph
  ###################################################################################
  
  # take the bottom half of the matrix, giving (n to n+1) and (n+1 to n+1) distances
  # as separate matrices
  n_n1 <- cluster_distances_matrix[(table_lengths[n]+1):(table_lengths[n]+table_lengths[n+1]),
          1:table_lengths[n]]
  n1_n1 <- cluster_distances_matrix[(table_lengths[n]+1):(table_lengths[n]+table_lengths[n+1]),
          (table_lengths[n]+1):(table_lengths[n]+table_lengths[n+1])]
  # convert the distance matrices to lists of coordinates with distances as the 3rd row
  edgelist_with_distances_n_n1 <- cbind(as.vector(row(n_n1)+table_breaks[n]),
          as.vector(col(n_n1)+table_breaks[n]-table_lengths[n]),as.vector(n_n1))
  edgelist_with_distances_n1_n1 <- cbind(as.vector(row(n1_n1)+table_breaks[n]),
          as.vector(col(n1_n1)+table_breaks[n]),as.vector(n1_n1))
  # only take upper triangle from n+1 to n+1 matrix
  edgelist_with_distances_n1_n1 <- edgelist_with_distances_n1_n1[upper.tri(matrix(data=1:length(n1_n1),
          nrow=nrow(n1_n1),ncol=ncol(n1_n1))),]
  #combine n_n1 and n1_n1 distances into a single "edgelist" with distances
  edgelist_with_distances <- rbind(edgelist_with_distances_n_n1,edgelist_with_distances_n1_n1)
  # convert strings to numeric
  class(edgelist_with_distances) <- "numeric"
  # sort edges by distance (column 3)
  edgelist_with_distances <- edgelist_with_distances[order(edgelist_with_distances[,3]),]
  # take only the top X-% edges by distance based on PERCENT_TOTAL
  trim_edgelist_with_distances <- edgelist_with_distances[1:floor(length(edgelist_with_distances[,1])*PERCENT_EDGES/100),]
  # add edges to graph
  # note: weight here is really distance, will be converted to weight after the graph is completed
  
  # calculate "density" for each cluster and normalize
  densities_no_zeros <- table(trim_edgelist_with_distances[,1:2])
  # add in zeros for clusters with no edges (table function leaves these out)
  densities <- rep(0,nrow(clusters))
  names(densities) <- (offset+1):(offset+table_lengths[n]+table_lengths[n+1])
  densities[names(densities_no_zeros)] <- densities_no_zeros
  normalized_densities <- round(densities/max(densities)*(MAX_EDGES - MIN_EDGES) + MIN_EDGES)
  
  
  # build new edgelist with N edges for each cluster based on normalized density
  cat("building edgelist\n")
  final_edgelist_with_distances <- c()
  for (i in 1:length(normalized_densities)) {
    tmp_final_edgelist_with_distances <- cbind(offset+i,order(cluster_distances_matrix[,i])+offset,
                                               sort(cluster_distances_matrix[,i]))[1:normalized_densities[i],]
    final_edgelist_with_distances <- rbind(final_edgelist_with_distances,tmp_final_edgelist_with_distances)
  }
  
  # remove duplicate edges from edgelist
  final_edgelist_with_distances <- unique(cbind(t(apply(final_edgelist_with_distances[,1:2],1,sort)),final_edgelist_with_distances[,3]))
  
  output_graph <- add.edges(output_graph,edges=as.vector(t(final_edgelist_with_distances[,1:2])-1),
          weight=final_edgelist_with_distances[,3],label="EXTRA",sequence_assignment=n+1)
  
  ###################################################################################
  # This section adds the "MST" for n_n+1 and n+1_n+1 nodes, a little complicated. . .
  ###################################################################################
  
  # create dummy graph contains n and n+1 vertices, but n_n section is fully connected so
  # it will act as a single connected component while n_n+1 and n+1_n+1 sections are not
  # connected, so each n+1 vertex will be its own disconnected component
  adjacency_nn <- cluster_distances_matrix[1:table_lengths[n],1:table_lengths[n]]
  el_nn_with_dist <- cbind(as.vector(row(adjacency_nn)),as.vector(col(adjacency_nn)),as.vector(adjacency_nn))
  dummy_graph <- graph.empty(n=nrow(cluster_distances_matrix), directed=FALSE)
  dummy_graph <- add.edges(dummy_graph,edges=as.vector(t(el_nn_with_dist[,1:2])-1))
  clu <- clusters(dummy_graph)
  # get the cluter that each node belongs to
  members <- clu$membership
  # initialize matrix to hold distances between graph components
  linkweights <- matrix(nrow=clu$no,ncol=clu$no)
  # loop through all components vs. all components to fill in distances
  for (i in 1:clu$no) {
    members_i <- which(members==i-1)
    for (j in 1:clu$no) {
      members_j <- which(members==j-1)
      # get the minimum distance between components i and j
      # (note members has zero index)
      tmpmatrix <- cluster_distances_matrix[members_i,members_j]
      linkweights[i,j] <- min(tmpmatrix)
    }
  }
  # set i-i distances to Inf instead of 0 so they aren't picked
  # as the closest neighbors
  for (i in 1:ncol(linkweights)) {
    linkweights[i,i] <- Inf
  }
  # make the minimum spanning tree for graph components
  components_adjacency <- graph.adjacency(linkweights,mode="undirected",weighted=TRUE,diag=TRUE)
  components_mst <- minimum.spanning.tree(components_adjacency)
  components_mst_el <- get.edgelist(components_mst,names=FALSE)
  
  # go from MST edges to the actual edges between individual nodes
  for (i in 1:nrow(components_mst_el)) {
    # get index of shortest connection between components
    # tmp <- which.min(adjacency[which(members==components_mst_el[i,1]),which(members==components_mst_el[i,2])])
    members_x <- which(members==components_mst_el[i,1])
    members_y <- which(members==components_mst_el[i,2])
    tmpmatrix <- as.matrix(cluster_distances_matrix[members_x,members_y])
    if(length(members_x)==1){
      tmpmatrix <- t(tmpmatrix)
    }
    tmp_min <- min(tmpmatrix) # get the largest weight = shortest distance
    tmp_index <- which(tmpmatrix==tmp_min,arr.ind=TRUE)
    
    # add new edge to graph with the shortest connection vertices and the weight from linkweights
    if (are.connected(output_graph,members_x[tmp_index[1]]-1+offset,members_y[tmp_index[2]]-1+offset)) {
      # cat(paste("\ntrue = ",i))
      E(output_graph,P=c(members_x[tmp_index[1]]-1+offset,members_y[tmp_index[2]]-1+offset))$label <- "MST"
    } else {
      output_graph <- add.edges(output_graph,c(members_x[tmp_index[1]]-1+offset,
              members_y[tmp_index[2]]-1+offset),weight=tmp_min,label="MST",sequence_assignment=n+1)
    }
  }
  
  # now that all edges are added and labeled for a-a and a-a+1,
  # update offset index so it will work for a+1-a+1 and a+1-a+2
  offset = offset+table_lengths[n]
}

#remove duplicate edges and loops just in case
#simplify(output_graph)

# convert graph distances to weights (low distance = high weight and vice versa)
distances <- E(output_graph)$weight
#weights <- -distances+max(distances)+min(distances)
weights <- 1/distances
E(output_graph)$weight <- weights

################################################################################
# This section annotates the graph
################################################################################

# get fcs files to use for annotation values
fcs_files = list.files(pattern=FCS_FILE_PATTERN,recursive=TRUE,full.names=TRUE)
fcs_files <- sort(fcs_files)

anno_cat <- c()
# iterate through all times and annotate
for (i in 1:length(fcs_files)) {
  
  # get fcs file corresponding to time[i]
  fcs_file = fcs_files[i]
  
  cat("Annotating graph for",fcs_file,"\n")
  
  # get medians (and count) for all parameters
  anno  <- SPADE.markerMedians(fcs_file, table_lengths[i])
  
  # add time[i] information column
  time_matrix <- matrix(i,nrow=length(anno[[1]]))
  colnames(time_matrix) <- c("Timepoint")
  anno$medians <- cbind(anno$medians,time_matrix)
    
  # add median and percent values
  for (a in c("medians", "percenttotal")) {

    # set rownames = the clusters for this time[i]
    # note: I am making this zero-indexed so it works better with zero-indexed graph later
    rownames(anno[[a]]) <- (table_breaks[i]-table_lengths[i]):(table_breaks[i]-1)
    
    # anno_cat is all the annos concatenated, will be used to make "anno.Rsave" file
    anno_cat[[a]] <- rbind(anno_cat[[a]],anno[[a]])
  }
}

output_anno <- c()
# combine anno_cat matrices
for (i in 1:length(anno_cat)) {
  output_anno <- cbind(output_anno,anno_cat[[i]])
}

# don't include barcode assignment and time (per cell) values because these are boring
output_anno <- output_anno[,colnames(output_anno)!="barcode"&colnames(output_anno)!="Time"]

for (c in colnames(output_anno)) {
  output_graph <- set.vertex.attribute(output_graph,c,index=as.numeric(rownames(output_anno)), value=output_anno[,c])
}

# add name attribute
V(output_graph)$name <- 1:length(V(output_graph))

write.graph(output_graph, paste("output_",MIN_EDGES,"_",MAX_EDGES,"_",PERCENT_EDGES,".graphml",sep=""), format="graphml")
