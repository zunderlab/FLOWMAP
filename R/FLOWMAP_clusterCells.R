
#' @export
FLOWMAPcluster <- function(full.clusters, table.breaks, table.lengths,
                           cluster.medians, cluster.counts, cell.assgn)  {
  object <- list(full.clusters = full.clusters,
                 table.breaks = table.breaks,
                 table.lengths = table.lengths,
                 cluster.medians = cluster.medians,
                 cluster.counts = cluster.counts,
                 cell.assgn = cell.assgn)
  # use only cluster_medians object and cellassgn ???
  # make method to concatenate cluster_medians into fullclusters
  # make method to get table breaks and table lengths from cluster_medians
  class(object) <- "FLOWMAPcluster"
  return (object)
}


ClusterFCS <- function(fcs.files, clustering.var, numcluster,
                       distance.metric) {
  full.clusters <- data.frame()
  table.breaks <- c()
  table.lengths <- c()
  cluster.medians <- list()
  cluster.counts <- list()
  cell.assgn <- list()
  if (length(numcluster) == 1) {
    cat("Clustering all files to:", numcluster, "\n")
    numcluster.new <- rep(numcluster, times = length(fcs.files))
    numcluster <- numcluster.new
  }
  if (length(numcluster) != length(fcs.files)) {
    stop("Cluster number not specified for all FCS files!")
  }
  for (i in 1:length(fcs.files)) {
    current.file <- fcs.files[[i]]
    cat("Clustering data from file", i, "\n")
    cluster.counts[[i]] <- data.frame()
    cluster.medians[[i]] <- data.frame()
    tmp.FCS.for.cluster <- subset(x = current.file, select = clustering.var)
    cat("Subsetting for clustering channels only", "\n")
    cluster.results <- HclustClustering(current.file = current.file, tmp.FCS.for.cluster = tmp.FCS.for.cluster,
                                        distance.metric = distance.metric, numcluster = numcluster[i])
    cell.assgn[[i]] <- cluster.results$tmp.cell.assgn
    cluster.medians[[i]] <- cluster.results$new.medians
    cluster.counts[[i]] <- cluster.results$new.counts
    colnames(cluster.counts[[i]]) <- c("Counts")
    # store the number of clusters in the file being read currently
    # used for node numbering/index
    table.lengths <- append(table.lengths, nrow(cluster.medians[[i]]))
    # store the currently read clusters and their median values to a master array
    full.clusters <- rbind(full.clusters, cluster.medians[[i]])
    # store where the array breaks between clusters from one file and those from the next file
    # used for node numbering/index
    table.breaks <- append(table.breaks, nrow(full.clusters))
  }
  for (i in 1:length(cluster.medians)) {
    rownames(cluster.medians[[i]]) <- seq(1, table.lengths[i])
  }    
  rownames(full.clusters) <- seq(1, dim(full.clusters)[1])
  FLOWMAP.clusters <- FLOWMAPcluster(full.clusters = full.clusters,
                                     table.breaks = table.breaks,
                                     table.lengths = table.lengths,
                                     cluster.medians = cluster.medians,
                                     cluster.counts = cluster.counts,
                                     cell.assgn = cell.assgn)
  return(FLOWMAP.clusters)  
}


MultiClusterFCS <- function(list.of.time.files, clustering.var, numcluster, distance.metric) {
  list.of.FLOWMAP.clusters <- list()
  numcluster.orig <- numcluster
  if (length(numcluster.orig) == 1) {
    cat("Clustering all files to:", numcluster, "\n")
  } 
  for (time in names(list.of.time.files)) {
    cat("Clustering all files from time", time, "\n")
    fcs.files <- list.of.time.files[[time]]
    if (length(numcluster.orig) == 1) {
      numcluster.new <- rep(numcluster.orig, times = length(fcs.files))
      numcluster <- numcluster.new
    } else {
      numcluster <- numcluster.orig[[time]]
    }
    if (length(numcluster) != length(fcs.files)) {
      stop("Cluster number not specified for all FCS files!")
    }
    file.clusters <- ClusterFCS(fcs.files, clustering.var, numcluster, distance.metric)
    list.of.FLOWMAP.clusters[[time]] <- file.clusters
  }
  return(list.of.FLOWMAP.clusters)
}


HclustClustering <- function(current.file, tmp.FCS.for.cluster, distance.metric, numcluster) {
  if (distance.metric == "euclidean") {
    method <- "ward"
  } else {
    method <- "single"
  }
  FCS.clusters <- Rclusterpp.hclust(tmp.FCS.for.cluster, method = method,
                                    distance = distance.metric)
  clust <- list(assgn = cutree(FCS.clusters, k = numcluster))
  new.counts <- data.frame()
  new.medians <- data.frame()
  # Invalid clusters have assgn == 0
  is.na(clust$assgn) <- which(clust$assgn == 0)
  for (p in c(1:max(clust$assgn, na.rm = TRUE))) {  
    obs <- which(clust$assgn == p)
    if (length(obs) > 1) {
      # Finding clusters containing more than one cell
      # Saving counts and medians of data
      new.counts <- rbind(new.counts, data.frame(length(obs)))
      new.median <- colMedians(as.matrix(current.file[obs, ]))
      new.median <- as.data.frame(new.median)
      if (!identical(colnames(new.median), colnames(current.file))) {
        new.median <- t(new.median)
        colnames(new.median) <- colnames(current.file)
      }
      # Saving cluster data
    } else {
      new.counts <- rbind(new.counts, data.frame(length(obs)))
      new.median <- current.file[obs, ]
      new.median <- as.data.frame(new.median)
      if (!identical(colnames(new.median), colnames(current.file))) {
        new.median <- t(new.median)
        colnames(new.median) <- colnames(current.file)
      }
    }
    new.medians <- rbind(new.medians, new.median)
  }
  tmp.cell.assgn <- data.frame(clust$assgn)
  tmp.cell.assgn <- as.data.frame(tmp.cell.assgn[complete.cases(tmp.cell.assgn), ])
  colnames(tmp.cell.assgn) <- c("Cluster")
  return(list(tmp.cell.assgn = tmp.cell.assgn,
              new.medians = new.medians,
              new.counts = new.counts))
}


