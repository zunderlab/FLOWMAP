
FLOWMAPcluster <- function(full.clusters, table.breaks, table.lengths,
                           cluster.medians, cluster.counts, cell.assgn)  {
  object <- list(full.clusters = full.clusters,
                 table.breaks = table.breaks,
                 table.lengths = table.lengths,
                 cluster.medians = cluster.medians,
                 cluster.counts = cluster.counts,
                 cell.assgn = cell.assgn)
  class(object) <- "FLOWMAPcluster"
  return (object)
}

RemodelFLOWMAPClusterList <- function(list.of.FLOWMAP.clusters) {
  # take FLOWMAP of conditions with timeseries and make into one timeseries
  # combine FLOWMAP conditions
  full.clusters <- data.frame()
  table.breaks <- c()
  table.lengths <- c()
  cluster.medians <- list()
  cluster.counts <- list()
  cell.assgn <- list()
  for (t in 1:length(list.of.FLOWMAP.clusters)) {
    temp.medians <- data.frame()
    temp.cell.assgn <- data.frame()
    temp.counts <- data.frame()
    for (c in 1:length(list.of.FLOWMAP.clusters[[t]]$cluster.medians)) {
      temp.medians <- rbind(temp.medians, list.of.FLOWMAP.clusters[[t]]$cluster.medians[[c]])
      temp.cell.assgn <- rbind(temp.cell.assgn, list.of.FLOWMAP.clusters[[t]]$cell.assgn[[c]])
      temp.counts <- rbind(temp.counts, list.of.FLOWMAP.clusters[[t]]$cluster.counts[[c]])
    }
    cluster.medians[[t]] <- temp.medians
    cluster.counts[[t]] <- temp.counts
    cell.assgn[[t]] <- temp.cell.assgn
    table.lengths <- c(table.lengths, dim(temp.medians)[1])
    table.breaks <- c(table.breaks, sum(table.lengths))
    full.clusters <- rbind(full.clusters, temp.medians)
  }
  remodeled.FLOWMAP.clusters <- FLOWMAPcluster(full.clusters, table.breaks, table.lengths,
                                               cluster.medians, cluster.counts, cell.assgn)
  return(remodeled.FLOWMAP.clusters)
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
    # cluster.results <- SPADEClustering(current.file = current.file, tmp.FCS.for.cluster = tmp.FCS.for.cluster,
    #                                     distance.metric = distance.metric, numcluster = numcluster[i])
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

MultiClusterFCS <- function(list.of.files, clustering.var, numcluster, distance.metric) {
  list.of.FLOWMAP.clusters <- list()
  numcluster.orig <- numcluster
  if (length(numcluster.orig) == 1) {
    cat("Clustering all files to:", numcluster, "\n")
  } 
  for (t in 1:length(list.of.files)) {
    cat("Clustering all files from time", t, "\n")
    fcs.files <- list.of.files[[t]]
    if (length(numcluster.orig) == 1) {
      numcluster.new <- rep(numcluster.orig, times = length(fcs.files))
      numcluster <- numcluster.new
    } else {
      numcluster <- numcluster.orig[[t]]
    }
    if (length(numcluster) != length(fcs.files)) {
      stop("Cluster number not specified for all FCS files!")
    }
    file.clusters <- ClusterFCS(fcs.files, clustering.var, numcluster, distance.metric)
    list.of.FLOWMAP.clusters[[t]] <- file.clusters
  }
  return(list.of.FLOWMAP.clusters)
}

HclustClustering <- function(current.file, tmp.FCS.for.cluster, distance.metric = "manhattan", numcluster) {
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

SPADEClustering <- function(current.file, tmp.FCS.for.cluster, distance.metric = "manhattan", numcluster) {
  if (nrow(tmp.FCS.for.cluster) > 60000) {
    warning("Potentially too many observations for the clustering step",immediate=TRUE);
  }
  if (nrow(tmp.FCS.for.cluster) < numcluster) {
    stop("Number of requested clusters exceeds number of events")
  }
  # Transpose table before call into row major order
  cluster <- Rclusterpp.hclust(tmp.FCS.for.cluster, distance = distance.metric)
  clust <- list(assgn = cutree(cluster, k = numcluster))
  # Invalid clusters have assgn == 0
  centers <- c()
  is.na(clust$assgn) <- which(clust$assgn == 0)
  for (i in c(1:max(clust$assgn, na.rm = TRUE))) {  
    obs <- which(clust$assgn == i)
    if (length(obs) > 1) {
      centers <- rbind(centers, colMeans(current.file[obs, , drop = FALSE]))
      clust$assgn[obs] <- nrow(centers)
    } else {
      is.na(clust$assgn) <- obs
    }
  }
  tmp.cell.assgn <- as.data.frame(clust$assgn)
  colnames(tmp.cell.assgn) <- c("Cluster")
  new.counts <- as.matrix(table(clust$assgn))
  new.counts <- as.data.frame(new.counts)
  colnames(new.counts) <- c("Counts")
  return(list(tmp.cell.assgn = tmp.cell.assgn,
              new.medians = centers,
              new.counts = new.counts))
}


Upsample_OLD <- function(FLOWMAP.clusters) {
  fixed.FLOWMAP.clusters <- FLOWMAP.clusters
  for (f in 1:length(FLOWMAP.clusters$cluster.counts)) {
    counts <- FLOWMAP.clusters$cluster.counts[[f]]$Counts
    densities <- FLOWMAP.clusters$cluster.medians[[f]]$density
    fixed.counts <- round((counts * densities))
    fixed.FLOWMAP.clusters$cluster.counts[[f]]$Counts <- fixed.counts
  }
  return(fixed.FLOWMAP.clusters)
}

Upsample <- function(file.names, FLOWMAP.clusters, var.remove, var.annotate) {
  fixed.FLOWMAP.clusters <- FLOWMAP.clusters
  for (f in 1:length(file.names)) {
    print("file.names[[f]]")
    print(file.names[[f]])
    fcs.file <- LoadCleanFCS(file.names[[f]], channel.remove = var.remove,
                             channel.annotate = var.annotate, subsamples = FALSE)
    print("class(fcs.file[[1]])")
    print(class(fcs.file[[1]])
    print("head(FLOWMAP.clusters$cell.assgn[[f]])")
    print(head(FLOWMAP.clusters$cell.assgn[[f]]))
    print("head(as.integer(FLOWMAP.clusters$cell.assgn[[f]]))")
    print(head(as.integer(FLOWMAP.clusters$cell.assgn[[f]])))
    all.cells.assign <- spade:::SPADE.assignToCluster(fcs.file[[1]], 
                                                      FLOWMAP.clusters$cluster.medians[[f]],
                                                      FLOWMAP.clusters$cell.assgn[[f]])
    print("all.cells.assign")
    print(all.cells.assign)
    fixed.counts <- table(all.cells.assign)
    print("fixed.counts")
    print(fixed.counts)
    fixed.FLOWMAP.clusters$cluster.counts[[f]]$Counts <- fixed.counts
  }
  stop("TESTING SPADE UPSAMPLING")
  return(fixed.FLOWMAP.clusters)
}

#### DF SPECIFIC METHODS ####

ConstructOneFLOWMAPCluster <- function(df) {
  full.clusters <- df
  cluster.medians <- list()
  cluster.medians[[1]] <- df
  table.breaks <- nrow(df)
  table.lengths <- nrow(df)
  cluster.counts <- list()
  cluster.counts[[1]] <- as.data.frame(matrix(rep(1, times = nrow(df))))
  colnames(cluster.counts[[1]]) <- c("Counts")
  cell.assgn <- list()
  cell.assgn[[1]] <- as.data.frame(matrix(seq(1, nrow(df))))
  colnames(cell.assgn[[1]]) <- c("Cluster")
  file.clusters <- FLOWMAPcluster(full.clusters, table.breaks, table.lengths,
                                  cluster.medians, cluster.counts, cell.assgn)
  return(file.clusters)
}

ConstructSingleFLOWMAPCluster <- function(list.of.df) {
  full.clusters <- data.frame()
  cluster.medians <- list()
  table.breaks <- c()
  table.lengths <- c()
  cluster.counts <- list()
  cell.assgn <- list()
  for (i in 1:length(list.of.df)) {
    full.clusters <- rbind(full.clusters, list.of.df[[i]])
    cluster.medians[[i]] <- list.of.df[[i]]
    table.lengths <- c(table.lengths, nrow(list.of.df[[i]]))
    table.breaks <- c(table.breaks, sum(table.lengths))
    cluster.counts[[i]] <- as.data.frame(matrix(rep(1, times = nrow(list.of.df[[i]]))))
    colnames(cluster.counts[[i]]) <- c("Counts")
    cell.assgn[[i]] <- as.data.frame(matrix(seq(1, nrow(list.of.df[[i]]))))
    colnames(cell.assgn[[i]]) <- c("Cluster")
  }
  file.clusters <- FLOWMAPcluster(full.clusters, table.breaks, table.lengths,
                                  cluster.medians, cluster.counts, cell.assgn)
  return(file.clusters)
}

ConstructMultiFLOWMAPCluster <- function(multi.list.df) {
  list.of.FLOWMAP.clusters <- list()
  for (i in 1:length(multi.list.df)) {
    list.of.FLOWMAP.clusters[[i]] <- ConstructSingleFLOWMAPCluster(multi.list.df[[i]])
  }
  return(list.of.FLOWMAP.clusters)
}
