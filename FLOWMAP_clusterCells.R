
FLOWMAPcluster <- function(fullclusters, table_breaks, table_lengths,
                           cluster_medians, cluster_counts, cellassgn,
                           ...)  {
  object <- list(fullclusters = fullclusters,
                 table_breaks = table_breaks,
                 table_lengths = table_lengths,
                 cluster_medians = cluster_medians,
                 cluster_counts = cluster_counts,
                 cellassgn = cellassgn,
                 ...)
  # use only cluster_medians object and cellassgn ???
  # make method to concatenate cluster_medians into fullclusters
  # make method to get table breaks and table lengths from cluster_medians
  class(object) <- "FLOWMAPcluster"
  return (object)
}

clusterFCS <- function(fcs_files, channel_cluster, numcluster,
                       distance_metric, ...) {
  fullclusters <- data.frame()
  table_breaks <- c()
  table_lengths <- c()
  cluster_medians <- list()
  cluster_counts <- list()
  cellassgn <- list()
  for (i in 1:length(fcs_files)) {
    currentfile <- fcs_files[[i]]
    # print FCS file that is being read
    cat("Clustering data from file", i, "\n")
    cluster_counts[[i]] <- data.frame()
    cluster_medians[[i]] <- data.frame()
    tmp_FCSforCluster <- subset(x = currentfile, select = channel_cluster)
    cat("Subsetting for clustering channels only", "\n")
    cluster_results <- hclustClustering(currentfile = currentfile, tmp_FCSforCluster = tmp_FCSforCluster,
                                        distance_metric = distance_metric, numcluster = numcluster)
    cellassgn[[i]] <- cluster_results$tmp_cellassgn
    cluster_medians[[i]] <- cluster_results$new_medians
    # cat("cluster medians are", dim(cluster_medians[[i]]), "\n")
    cluster_counts[[i]] <- cluster_results$new_counts
    colnames(cluster_counts[[i]]) <- c("Counts")
    # cat("cluster counts are", dim(cluster_counts[[i]]), "\n")
    # store the number of clusters in the file being read currently
    # used for node numbering/index
    table_lengths <- append(table_lengths, nrow(cluster_medians[[i]]))
    # cat("table lengths are", table_lengths, "\n")
    # store the currently read clusters and their median values to a master array
    fullclusters <- rbind(fullclusters, cluster_medians[[i]])
    # cat("dim of full clusters are", dim(fullclusters), "\n")
    # store where the array breaks between clusters from one file and those from the next file
    # used for node numbering/index
    table_breaks <- append(table_breaks, nrow(fullclusters))
    # cat("table breaks are", table_breaks, "\n")
  }
  for (i in 1:length(cluster_medians)) {
    rownames(cluster_medians[[i]]) <- seq(1, table_lengths[i])
  }    
  rownames(fullclusters) <- seq(1, dim(fullclusters)[1])
  FLOWMAPclusters <- FLOWMAPcluster(fullclusters = fullclusters,
                                    table_breaks = table_breaks,
                                    table_lengths = table_lengths,
                                    cluster_medians = cluster_medians,
                                    cluster_counts = cluster_counts,
                                    cellassgn = cellassgn)
  return(FLOWMAPclusters)  
}


multiClusterFCS <- function(listOfTreatFiles, channel_cluster, numcluster,
                            clustalgorithm = "clara") {
  listofFLOWMAPclusters <- list()
  for (treat in names(listOfTreatFiles)) {
    cat("Clustering all files from", treat, "\n")
    fcs_files <- listOfTreatFiles[[treat]]
    file_clusters <- clusterFCS(fcs_files, channel_cluster,
                                numcluster, clustalgorithm = "clara")
    listofFLOWMAPclusters[[treat]] <- file_clusters
  }
  return(listofFLOWMAPclusters)
}


hclustClustering <- function(currentfile, tmp_FCSforCluster, distance_metric, numcluster) {
  # print("currentfile")
  # print(head(currentfile))
  # print(colnames(currentfile))
  if (distance_metric == "euclidean") {
    method <- "ward"
  } else {
    method <- "single"
  }
  FCSclusters <- hclust.vector(tmp_FCSforCluster, method = method,
                               metric = distance_metric)
  # print("tmp_FCSforCluster")
  # print(head(tmp_FCSforCluster))
  # print(colnames(tmp_FCSforCluster))
  clust <- list(assgn = cutree(FCSclusters, k = numcluster))
  new_counts <- data.frame()
  new_medians <- data.frame()
  # Invalid clusters have assgn == 0
  is.na(clust$assgn) <- which(clust$assgn == 0)
  for (p in c(1:max(clust$assgn, na.rm = TRUE))) {  
    obs <- which(clust$assgn == p)
    # cat("length(obs) is", length(obs), "\n")
    if (length(obs) > 1) {
      # Finding clusters containing more than one cell
      # Saving counts and medians of data
      new_counts <- rbind(new_counts, data.frame(length(obs)))
      new_median <- colMedians(as.matrix(currentfile[obs, ]))
      new_median <- as.data.frame(new_median)
      # print("more than one cell")
      # print(new_median)
      if(colnames(new_median) != colnames(currentfile)) {
        new_median <- t(new_median)
        colnames(new_median) <- colnames(currentfile)
      }
      # Saving cluster data
    } else {
      new_counts <- rbind(new_counts, data.frame(length(obs)))
      new_median <- currentfile[obs, ]
      # print("only one cell")
      # print(new_median)
      new_median <- as.data.frame(new_median)
      if(colnames(new_median) != colnames(currentfile)) {
        new_median <- t(new_median)
        colnames(new_median) <- colnames(currentfile)
      }
    }
    # print("woof")
    # cat("colnames(new_medians) is", colnames(new_medians), "\n")
    # cat("colnames(new_median) is", colnames(new_median), "\n")
    new_medians <- rbind(new_medians, new_median)
    # print("meow")
  }
  tmp_cellassgn <- data.frame(clust$assgn)
  tmp_cellassgn <- as.data.frame(tmp_cellassgn[complete.cases(tmp_cellassgn), ])
  colnames(tmp_cellassgn) <- c("Cluster")
  return(list(tmp_cellassgn = tmp_cellassgn,
              new_medians = new_medians,
              new_counts = new_counts))
}



#### UNFINISHED FUNCTION ####
flowmapFCSFiles <- function(fcs_files, FLOWMAPclusters, output_folder) {
  # write CSV of FCS files with flowmap cluster assignments
  for (i in 1:length(fcs_files)) {
    flowmap_fcs_files <- list()
    FLOWMAPclusters$cellassgn
  }
}


