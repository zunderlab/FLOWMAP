library(Rclusterpp)
library(cluster)
library(kernlab)


getFilesWithInfo <- function(loaded_fcs_files, FLOWMAPClusters_stages, cellassgn) {
  fcs_files_with_clusters <- getFilesWithClusters(loaded_fcs_files, cellassgn)
  fcs_files_with_stage <- getFilesWithStages(fcs_files_with_clusters, FLOWMAPClusters_stages)
  return(fcs_files_with_stage)
}


getFilesWithStages <- function(fcs_files_with_clusters, FLOWMAPClusters_stages) {
  #   num_stages <- 0
  #   for (i in 1:length(FLOWMAPClusters_stages)) {
  #     file_max <- max(FLOWMAPClusters_stages[[i]]$Stage)
  #     num_stages <- max(num_stages, file_max)
  #   }
  fcs_files_with_stage <- fcs_files_with_clusters
  for (i in 1:length(fcs_files_with_clusters)) {
    num_clus <- max(fcs_files_with_clusters[[i]]$Cluster)
    ncell <- dim(fcs_files_with_clusters[[i]])[1]
    Stage <- rep(NA, times = ncell)
    for (n in 1:num_clus) {
      ind <- which(fcs_files_with_clusters[[i]]$Cluster == n)
      stage_assgn <- FLOWMAPClusters_stages[[i]]$Stage[n]
      Stage[ind] <- stage_assgn
    }
    fcs_files_with_stage[[i]] <- cbind(fcs_files_with_clusters[[i]], Stage)
  }
  return(fcs_files_with_stage)
}


getFilesWithClusters <- function(loaded_fcs_files, cellassgn) {
  for (i in 1:length(loaded_fcs_files)) {
    numcells <- dim(loaded_fcs_files[[i]])[1]
    Cluster <- rep(NA, times = numcells)
    numclus <- max(cellassgn[[i]]$Cluster)
    for (clus in 1:numclus) {
      ind <- which(cellassgn[[i]]$Cluster == clus)
      Cluster[ind] <- clus
    }
    loaded_fcs_files[[i]] <- cbind(loaded_fcs_files[[i]], Cluster)
  }
  fcs_files_with_clusters <- loaded_fcs_files
  return(fcs_files_with_clusters)
}


getClusFiles <- function(fcs_file_with_clusters, clust_num) {
  ind <- which(fcs_file_with_clusters$Cluster == clust_num)
  cluster_fcs_file <- fcs_file_with_clusters[ind, ]
  return(cluster_fcs_file)
}


getAllClusFiles <- function(fcs_file_with_clusters) {
  cluster_fcs_files <- list()
  clust_num <- max(fcs_file_with_clusters$Cluster)
  for (clus in 1:clust_num) {
    cluster_fcs_files[[clus]] <- getClusFiles(fcs_file_with_clusters, clus)
  }
  return(cluster_fcs_files)
}


writeFLOWMAPFCSFiles <- function(loaded_fcs_files, FLOWMAPclusters, output_folder) {
  # write CSV of FCS files with flowmap cluster assignments
  cellassgn <- FLOWMAPclusters$cellassgn
  fcs_files_with_clusters <- getFilesWithClusters(loaded_fcs_files, cellassgn)
  all_fcs_files <- list()
  out_folder <- paste(output_folder, "/clusterfiles/", sep = "")
  cat("Making output folder:", out_folder)
  dir.create(out_folder)
  for (i in 1:length(fcs_files_with_clusters)) {
    base_filename <- paste(output_folder, "/", "file", i, sep = "")
    cluster_fcs_files <- getAllClusFiles(fcs_files_with_clusters[[i]])
    for (f in 1:length(cluster_fcs_files)) {
      filename <- paste(base_filename, "_cluster", i, ".csv", sep = "")
      write.csv(cluster_fcs_files[[f]], filename)
    }
    all_fcs_files[[i]] <- cluster_fcs_files
  }
  return(all_fcs_files)
}


writeStageFCSFiles <- function(fcs_files_with_info, output_folder) {
  # fcs_files_with_info has stage and cluster info per cell
  # write CSV of FCS files with flowmap stage assignments
  num_stages <- 0
  for (i in 1:length(fcs_files_with_info)) {
    file_max <- max(fcs_files_with_info[[i]]$Stage)
    num_stages <- max(num_stages, file_max)
  }
  all_fcs_files <- list()
  out_folder <- paste(output_folder, "/stagefiles/", sep = "")
  cat("Making output folder:", out_folder)
  dir.create(out_folder)
  for (i in 1:num_stages) {
    filename <- paste(out_folder, "/", "stage", i, ".csv", sep = "")
    stage_fcs_file <- data.frame()
    for (f in 1:length(fcs_files_with_info)) {
      ind <- which(fcs_files_with_info[[f]]$Stage == i)
      newcells <- fcs_files_with_info[[f]][ind, ]
      File <- rep(f, times = length(ind))
      newcells <- cbind(newcells, File)
      stage_fcs_file <- rbind(stage_fcs_file, newcells)
    }
    all_fcs_files[[i]] <- stage_fcs_file
    write.csv(stage_fcs_file, filename)
    rm(stage_fcs_file)
  }
#   return(all_fcs_files)
}




#   for (i in 1:length(fcs_files)) {
#     x <- max(FLOWMAPclusters[[i]]$cellassgn)
#     for (num in 1:x) {
#       ind <- FLOWMAPclusters[[i]]$cellassgn == num
#       
#     }
#     flowmap_fcs_files <- list()
#     FLOWMAPclusters[[i]]$cellassgn
#   }
# }
# 


#   numclusters <- 0
#   for (l in 1:length(fcs_files)) {
#     numclusters <- sum(numclusters, length(file_clusters$cellassgn))
#   }

#   for (i in 1:length(file_clusters$cellassgn)) {
#     # numrow <- nrow(file_clusters$cellassgn[[i]])
#     numclus <- max(file_clusters$cellassgn[[i]]$Cluster)
#     for (clus in 1:numclus) {
#       ind <- which(file_clusters$cellassgn[[i]]$Cluster == clus)
#       cluscells <- loaded_fcs_files[[i]][ind, ]
#     }
#   }


