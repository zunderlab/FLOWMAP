rm(list=ls(all=T))
library(spade)


DISTANCE_METHOD <- "manhattan"
CLUSTER_TABLE_PATTERN <- "clusters.table$"
CLUSTER_FCS_PATTERN <- ".density.fcs.cluster.fcs$"
PERCENT_EDGES <- c(1,2,5)
MAX_EDGES <- c(5,10,15)
MIN_EDGES <- c(2,3)


# get clusters from file
clusters_table_files <- list.files(pattern=CLUSTER_TABLE_PATTERN,recursive=TRUE,full.names=TRUE)

for (c in clusters_table_files) {
  tmp_dir <- sub(paste("/",CLUSTER_TABLE_PATTERN,sep=""),"",c)
  
  cat("calculating densities\n")
  # get distance matrix from clusters
  clusters <- read.table(c, header=TRUE, check.names=FALSE)
  cluster_distances <- dist(clusters, method=DISTANCE_METHOD)
  cluster_distances_matrix <- as.matrix(cluster_distances)
  # set i-i distances to Inf instead of 0 so they aren't the closest neighbors
  for (i in 1:ncol(cluster_distances_matrix)) {
    cluster_distances_matrix[i,i] <- Inf
  }
  
  # build initial graph
  adjacency_graph <- graph.adjacency(cluster_distances_matrix, mode="undirected", weighted="weight", diag=FALSE)
  mst_graph <- minimum.spanning.tree(adjacency_graph, algorithm="prim")
  # for each edge of the mst, if it exists in the graph then label MST, if it doesn't exist then add it
  mst_graph_edgelist <- cbind(t(apply(get.edgelist(mst_graph),1,as.numeric)),E(mst_graph)$weight)
  # for each edge of the mst, if it exists in the graph then label MST, if it doesn't exist then add it
  
  
  # make edgelist with distances as the third column
  edgelist_with_distances <- cbind(as.vector(row(cluster_distances_matrix)),
              as.vector(col(cluster_distances_matrix)),as.vector(cluster_distances_matrix))
  edgelist_with_distances <- edgelist_with_distances[upper.tri(matrix(data=1:length(cluster_distances_matrix),
              nrow=nrow(cluster_distances_matrix),ncol=ncol(cluster_distances_matrix))),]
  # sort edges by distance (column 3)
  edgelist_with_distances <- edgelist_with_distances[order(edgelist_with_distances[,3]),]
  
  for (per in PERCENT_EDGES) {
    for (min in MIN_EDGES) {
      for (max in MAX_EDGES) {
        # take only the top X-% edges by distance based on PERCENT_TOTAL
        trim_edgelist_with_distances <- edgelist_with_distances[1:floor(length(edgelist_with_distances[,1])*per/100),]
        
        # calculate "density" for each cluster and normalize
        densities_no_zeros <- table(trim_edgelist_with_distances[,1:2])
        # add in zeros for clusters with no edges (table function leaves these out)
        densities <- rep(0,nrow(clusters))
        densities[as.numeric(names(densities_no_zeros))] <- densities_no_zeros
        normalized_densities <- round(densities/max(densities)*(max - min) + min)
        
        
        # build new edgelist with N edges for each cluster based on normalized density
        cat("building edgelist\n")
        final_edgelist_with_distances <- c()
        for (n in 1:length(normalized_densities)) {
          tmp_final_edgelist_with_distances <- cbind(n,order(cluster_distances_matrix[,n]),
                                                     sort(cluster_distances_matrix[,n]))[1:normalized_densities[n],]
          final_edgelist_with_distances <- rbind(final_edgelist_with_distances,tmp_final_edgelist_with_distances)
        }
        
        # remove duplicate edges from edgelist
        final_edgelist_with_distances <- unique(cbind(t(apply(final_edgelist_with_distances[,1:2],1,sort)),final_edgelist_with_distances[,3]))
        
        # make a new graph from the edgelist and add distances (here they are called "weights"
        # but they are really distances - will convert them after MST edges are added
        map_graph <- graph.edgelist(final_edgelist_with_distances[,1:2], directed=FALSE)
        distances <- final_edgelist_with_distances[,3]
        # weights <- -distances+max(distances)+min(distances)
        E(map_graph)$weight <- distances
        E(map_graph)$label <- "EXTRA"
        
        for (i in 1:nrow(mst_graph_edgelist)) {
          if (are.connected(map_graph,mst_graph_edgelist[i,1],mst_graph_edgelist[i,2])) {
            E(map_graph,P=c(mst_graph_edgelist[i,1],mst_graph_edgelist[i,2]))$label <- "MST"
          } else {
            map_graph <- add.edges(graph=map_graph,edges=mst_graph_edgelist[i,1:2],weight=mst_graph_edgelist[i,3],label="MST")
          }
        }
        # remove the "0" node
        map_graph <- delete.vertices(map_graph,c(0))
        
        #convert distances to weight
        E(map_graph)$weight <- 1/as.numeric(E(map_graph)$weight)
        
        # add name attribute
        V(map_graph)$name <- 1:length(V(map_graph))
        
        for (f in list.files(path=tmp_dir,pattern=CLUSTER_FCS_PATTERN,full.names=TRUE)) {
          cat("Annotating medians for",f,"\n")
          
          # get medians and percenttotal
          anno <- SPADE.markerMedians(f,nrow(clusters))
          anno_cat <- cbind(anno[["medians"]],anno[["percenttotal"]])
          
          tmp_graph <- map_graph  
          rownames(anno_cat) <- 1:(nrow(clusters))  
          for (n in colnames(anno_cat)) {
            tmp_graph <- set.vertex.attribute(tmp_graph,n,index=as.numeric(rownames(anno_cat))-1, value=anno_cat[,n])
          }
          
          # write out graph in graphml format and use GEPHI command line API to make new layout
          output_name <- paste(f,".medians_",per,"_",min,"_",max,".graphml",sep="")
          write.graph(tmp_graph, output_name ,format="graphml")
        }
      }
    }
  }
}