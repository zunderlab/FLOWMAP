
RunForceDirectedLayout <- function(mode, file.name, graph, orig.times=NULL, which.palette=NULL) {
  file.name <- paste(unlist(strsplit(basename(file.name), "\\."))[1], "FLOW-MAP", sep = "_")
  ConvertToGraphML(output.graph = graph, file.name = file.name)
  graph.xy <- ForceDirectedXY(graph = graph)
  file.name.xy <- paste(file.name, "xy", sep = "_")
  final.file.name <- ConvertToGraphML(output.graph = graph.xy, file.name = file.name.xy)
  fixed.file.name <- paste(file.name.xy, "orig_time", sep = "_")
  if (mode != "one" && mode != "one-special") {
    fixed.graph <- ConvertOrigTime(graph.xy, orig.times)
  } else {
    fixed.graph <- graph.xy
  }
  fixed.file <- ConvertToGraphML(output.graph = fixed.graph, file.name = fixed.file.name)
  ExportClusterTables(output.graph = fixed.graph, file.name = fixed.file.name)
  MakeFLOWMAPRFile(env = parent.frame())
  if (savePDFs) {
    cat("Printing pdfs.", "\n")
    ConvertToPDF(graphml.file = fixed.file, which.palette = which.palette)
  }
  return(graph.xy)
}


#' Generate UMAP layout for FLOWMAP graph
#' @import ggplot2
#' @import ggfortify
RunUMAPlayout <- function(graph, file.clusters, file.name=file.name, umap.settings) { #file.clusters, knn.indexes, knn.distances

  ###NEED OUTPUT KNN IN FORMAT:
  #  knn = list(indexes=rbind(umap$knn$indexes, spectator.knn$indexes),
  #             distances=rbind(umap$knn$distances, spectator.knn$distances))


  # global.edgelist.save.df <- data.frame(global.edgelist.save)
  # index.list <- list()
  # dist.list <- list()
  # for (i in unique(global.edgelist.save.df$X0.row.inds)) {
  #   index.list[[i]] <- global.edgelist.save.df[ which(global.edgelist.save.df$X0.row.inds==i), "X0.col.inds"]
  #   dist.list[[i]] <- global.edgelist.save.df[ which(global.edgelist.save.df$X0.row.inds==i), "X0.values"]
  # }

  # #Read in knn info into correct format for umap
  # umap.settings$knn <- list(indexes=knn.indexes,
  #                          distances=knn.distances)
  # #class(umap.settings$knn) = "umap.knn"

  # #Read in df with median values for clusters from all files/timepoints
  # clusters.data <- file.clusters$full.clusters
  # #Add column to identify file origin/timepoint
  # diff.file.var <- c()
  # for (i in 1:length(file.clusters$table.lengths)) {
  #   diff.file.var <- append(diff.file.var, rep(i, file.clusters$table.lengths[i]))
  # }
  # clusters.data$file_var <- diff.file.var
  # umap.in <- clusters.data

  edgelist.with.weights <- cbind( get.edgelist(graph) , round( E(graph)$weight, 3 ))

  full.dims <- length(unique(append(edgelist.with.weights[,1], edgelist.with.weights[,2])))

  matSparse <- Matrix::sparseMatrix(i = edgelist.with.weights[,1],
                            j = edgelist.with.weights[,2],
                            x = edgelist.with.weights[,3],
                            dims = c(full.dims, full.dims),
                            dimnames = list(1:full.dims, 1:full.dims))

  input_dists <- Matrix::as.matrix(matSparse)
  #input_dists <- output.graph.dist.mat

  ## Use distance matrix to generate umap layout of clusters=======================================
  # currently "dist" are actually weights (this is what is reqd for igraph)
  # set the self-self distances back to 1, then do 1.0001-weights to get distances (so
  # nothing ends up being zero)
  filler.index <- which(input_dists==0,arr.ind = T)
  input_dists[filler.index] <-  0.0001 #1/0.0001=10000, much higher than any actual distances, so these will not be neighbors
  input_dists <- 1/input_dists #to get back to dists, was weight for igraph/force directed
  global.input.dists <<- input_dists
  #Run UMAP
  print("running UMAP")
  # run umap with default setting except for specifying that input is distance matrix
  umap.settings$input <- 'dist'
  umap.out <- umap(input_dists,config = umap.settings)

  print("ran UMAP, outputting files")

  #Save UMAP layout and knn graph
  UMAP.LAYOUT.FILENAME <- "UMAP_layout.csv"
  UMAP.INDEXES.FILENAME <- "UMAP_knn_indexes.csv"
  UMAP.DISTANCES.FILENAME <- "UMAP_knn_distances.csv"
  data.table::fwrite(umap.out$knn$indexes,file = UMAP.INDEXES.FILENAME,row.names = FALSE,col.names = FALSE,
         sep = ",")
  data.table::fwrite(umap.out$knn$distances,file = UMAP.DISTANCES.FILENAME,row.names = FALSE,col.names = FALSE,
         sep = ",")
  colnames(umap.out$layout) <- c("umap_x","umap_y")
  data.table::fwrite(umap.out$layout,file = UMAP.LAYOUT.FILENAME,row.names = FALSE,col.names = TRUE, sep = ",")

  #Add file var to umap layout
  diff.file.var <- c()
  for (i in 1:length(file.clusters$table.lengths)) {
    diff.file.var <- append(diff.file.var, rep(i, file.clusters$table.lengths[i]))
  }
  umap.layout <- data.frame(cbind(umap.out$layout, diff.file.var))
  global.umap.layout <<- umap.layout
  colnames(umap.layout) <- c("umap_x","umap_y","file_var")
  umap.layout$file_var <- as.factor(umap.layout$file_var)

  #Plot UMAP
  POINT.SIZE <- 1
  PLOT.HEIGHT <- 7
  PLOT.WIDTH <- 7

  ggplot2::ggsave(paste0("basic_umap",".png"),
                  plot = ggplot2::ggplot(global.umap.layout,
                                         ggplot2::aes_string(x="umap_x",y="umap_y",color="file_var")) + #factor()
                    ggplot2::geom_point(size = POINT.SIZE) +
                    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                   panel.grid.minor = ggplot2::element_blank(),
                                   panel.background = ggplot2::element_blank(),
                                   axis.line = ggplot2::element_line(colour = "black")),
                  height = PLOT.HEIGHT, width = PLOT.WIDTH)
}









