library(Rclusterpp)
library(cluster)

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
                       clustalgorithm = "clara",
                       ...) {
  fullclusters <- data.frame()
  table_breaks <- c()
  table_lengths <- c()
  cluster_medians <- list()
  cluster_counts <- list()
  cellassgn <- list()
  SAMPLE_NUM <- max(length(fcs_files) * 2, 5)
  SAMPLE_SIZE <- numcluster * 2    
  if (is.data.frame(fcs_files) == TRUE) {
    print("Only 1 file given in fcs_files")
    currentfile <- fcs_files
    i <- 1
    cat("Clustering data from file.\n")
    cluster_counts[[i]] <- data.frame()
    cluster_medians[[i]] <- data.frame()
    if (CLUSTERING_VAR == "all") {
      tmp_FCSforCluster <- currentfile
    } else {
      tmp_FCSforCluster <- subset(currentfile, select = colnames(currentfile)[colnames(currentfile) %in% channel_cluster])
    }
    if (clustalgorithm == "clara") {
      DISTANCE_METHOD <- "manhattan"
      FCSclusters <- clara(tmp_FCSforCluster, numcluster,
                           metric = DISTANCE_METHOD, 
                           samples = SAMPLE_NUM, sampsize = SAMPLE_SIZE)
      for (n in 1:numcluster) {
        nclustdata <- currentfile[FCSclusters$clustering == n, ]
        newmedian <- as.data.frame(colMedians(as.matrix(nclustdata)))
        if(colnames(newmedian) != colnames(currentfile)) {
          newmedian <- t(newmedian)
        }
        cluster_medians[[i]] <- rbind(cluster_medians[[i]], newmedian)
        cluster_counts[[i]] <- rbind(cluster_counts[[i]], sum(FCSclusters$clustering == n))
      }
      tmp_cellassgn <- data.frame(FCSclusters$clustering)
      tmp_cellassgn <- as.data.frame(tmp_cellassgn[complete.cases(tmp_cellassgn),])
      colnames(tmp_cellassgn) <- c("Cluster")
    }
    if (clustalgorithm == "hclust") {
      FCSclusters <- Rclusterpp.hclust(tmp_FCSforCluster)
      clust <- list(assgn = cutree(FCSclusters, k = numcluster))
      # Invalid clusters have assgn == 0
      centers = c()
      is.na(clust$assgn) <- which(clust$assgn == 0)
      for (p in c(1:max(clust$assgn, na.rm = TRUE))) {  
        obs <- which(clust$assgn == p)
        if (length(obs) > 1) {
          # Finding clusters containing more than one cell
          # Saving counts and medians of data
          cluster_counts[[i]] <- rbind(cluster_counts[[i]], data.frame(length(obs)))
          cluster_medians[[i]] <- rbind(cluster_medians[[i]], colMedians(currentfile[obs, keep.names = TRUE]))
          # Saving cluster data
        } else {
          cluster_counts[[i]] <- rbind(cluster_counts[[i]], data.frame(length(obs)))
          cluster_medians[[i]] <- rbind(cluster_medians[[i]], currentfile[obs, ])
        }
      }
      tmp_cellassgn <- data.frame(clust$assgn)
      tmp_cellassgn <- as.data.frame(tmp_cellassgn[complete.cases(tmp_cellassgn),])
      colnames(tmp_cellassgn) <- c("Cluster")
    }
    cellassgn[[i]] <- tmp_cellassgn
    # store the number of clusters in the file being read currently
    # used for node numbering/index
    table_lengths <- append(table_lengths, nrow(cluster_medians[[i]]))
    # store the currently read clusters and their median values to a master array
    fullclusters <- cluster_medians[[i]]
    # store where the array breaks between clusters from one file and those from the next file
    # used for node numbering/index
    table_breaks <- append(table_breaks, nrow(fullclusters)) 
  }
  
  else {
    for (i in 1:length(fcs_files)) {
      currentfile <- fcs_files[[i]]
      # print FCS file that is being read
      cat("Clustering data from file", i, "\n")
      cluster_counts[[i]] <- data.frame()
      cluster_medians[[i]] <- data.frame()
      if (CLUSTERING_VAR == "all") {
        tmp_FCSforCluster <- currentfile
      } else {
        tmp_FCSforCluster <- subset(currentfile, select = colnames(currentfile)[colnames(currentfile) %in% channel_cluster])
      }
      if (clustalgorithm == "clara") {
        DISTANCE_METHOD <- "manhattan"
        FCSclusters <- clara(tmp_FCSforCluster, numcluster,
                             metric = DISTANCE_METHOD, 
                             samples = SAMPLE_NUM, sampsize = SAMPLE_SIZE)
        for (n in 1:numcluster) {
          nclustdata <- currentfile[FCSclusters$clustering == n, ]
          if (!is.na(match("Treat", colnames(nclustdata)))) {
            Treat <- nclustdata[, "Treat"]
#             print(head(Treat))
            nclustdata <- subset(nclustdata, select = colnames(nclustdata)[colnames(nclustdata) != "Treat"])
          }
          newmedian <- as.data.frame(colMedians(as.matrix(nclustdata)))
          if(colnames(newmedian) != colnames(currentfile)) {
            newmedian <- t(newmedian)
          }
          newmedian <- cbind(newmedian, Treat)
#           print(Treat)
#           print(colnames(newmedian))
          cluster_medians[[i]] <- rbind(cluster_medians[[i]], newmedian)
          cluster_counts[[i]] <- rbind(cluster_counts[[i]], sum(FCSclusters$clustering == n))
        }
        tmp_cellassgn <- data.frame(FCSclusters$clustering)
        tmp_cellassgn <- as.data.frame(tmp_cellassgn[complete.cases(tmp_cellassgn),])
        colnames(tmp_cellassgn) <- c("Cluster")
      }
      if (clustalgorithm == "hclust") {
        FCSclusters <- Rclusterpp.hclust(tmp_FCSforCluster)
        clust <- list(assgn = cutree(FCSclusters, k = numcluster))
        # Invalid clusters have assgn == 0
        centers = c()
        is.na(clust$assgn) <- which(clust$assgn == 0)
        for (p in c(1:max(clust$assgn, na.rm = TRUE))) {  
          obs <- which(clust$assgn == p)
          if (length(obs) > 1) {
            # Finding clusters containing more than one cell
            # Saving counts and medians of data
            cluster_counts[[i]] <- rbind(cluster_counts[[i]], data.frame(length(obs)))
            cluster_medians[[i]] <- rbind(cluster_medians[[i]], colMedians(currentfile[obs, keep.names = TRUE]))
            # Saving cluster data
          } else {
            cluster_counts[[i]] <- rbind(cluster_counts[[i]], data.frame(length(obs)))
            cluster_medians[[i]] <- rbind(cluster_medians[[i]], currentfile[obs, ])
          }
        }
        tmp_cellassgn <- data.frame(clust$assgn)
        tmp_cellassgn <- as.data.frame(tmp_cellassgn[complete.cases(tmp_cellassgn),])
        colnames(tmp_cellassgn) <- c("Cluster")
      }
      cellassgn[[i]] <- tmp_cellassgn
      colnames(cluster_counts[[i]]) <- c("Counts")
      # store the number of clusters in the file being read currently
      # used for node numbering/index
      table_lengths <- append(table_lengths, nrow(cluster_medians[[i]]))
      # store the currently read clusters and their median values to a master array
      fullclusters <- rbind(fullclusters, cluster_medians[[i]])
      # store where the array breaks between clusters from one file and those from the next file
      # used for node numbering/index
      table_breaks <- append(table_breaks, nrow(fullclusters)) 
    }
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


#### UNFINISHED FUNCTION ####
flowmapFCSFiles <- function(fcs_files, FLOWMAPclusters, output_folder) {
  # write CSV of FCS files with flowmap cluster assignments
  for (i in 1:length(fcs_files)) {
    flowmap_fcs_files <- list()
    FLOWMAPclusters$cellassgn
  }
}


