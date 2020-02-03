
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
RunUMAPlayout <- function(graph, file.clusters, file.name=file.name, umap.settings) {

  ## Generate dist matrix from FLOWMAP graph
  # non-existent edges will be set to high distance
  # to reflect density-based n neighbors
  edgelist.with.weights <- cbind( get.edgelist(graph) , round( E(graph)$weight, 3 ))
  full.dims <- length(unique(append(edgelist.with.weights[,1], edgelist.with.weights[,2])))
  matSparse <- Matrix::sparseMatrix(i = edgelist.with.weights[,1],
                            j = edgelist.with.weights[,2],
                            x = edgelist.with.weights[,3],
                            dims = c(full.dims, full.dims),
                            dimnames = list(1:full.dims, 1:full.dims))
  input_dists <- Matrix::as.matrix(matSparse)

  ## Use distance matrix to generate umap layout of clusters
  # currently "dist" are actually weights (this is what is reqd for igraph)
  # set the self-self distances back to 1, then do 1.0001-weights to get distances (so
  # nothing ends up being zero)
  filler.index <- which(input_dists==0,arr.ind = T)
  input_dists[filler.index] <-  0.0001 #1/0.0001=10000, much higher than any actual distances, so these will not be neighbors
  input_dists <- 1/input_dists #to get back to dists, was weight for igraph/force directed
  global.input.dists <<- input_dists
  ## Run UMAP
  cat("running UMAP\n")
  umap.settings$input <- 'dist' #run umap with input settings except for specifying that input is distance matrix
  umap.out <- umap(input_dists,config = umap.settings)
  cat("ran UMAP, outputting files\n")

  ## Save UMAP layout and knn graph
  data.table::fwrite(umap.out$knn$indexes,file = "UMAP_knn_indexes.csv",row.names = FALSE,col.names = FALSE,sep = ",")
  data.table::fwrite(umap.out$knn$distances,file = "UMAP_knn_distances.csv",row.names = FALSE,col.names = FALSE,sep = ",")
  colnames(umap.out$layout) <- c("umap_x","umap_y")
  data.table::fwrite(umap.out$layout,file = "UMAP_layout.csv",row.names = FALSE,col.names = TRUE, sep = ",")

  ## Add file var to umap layout
  #TODO identify better way to track timepoint/other potentially needed info?
  diff.file.var <- c()
  for (i in 1:length(file.clusters$table.lengths)) {
    diff.file.var <- append(diff.file.var, rep(i, file.clusters$table.lengths[i]))
  }
  umap.layout <- data.frame(cbind(umap.out$layout, diff.file.var))
  colnames(umap.layout) <- c("umap_x","umap_y","file_var")
  umap.layout$file_var <- as.factor(umap.layout$file_var)

  #Plot UMAP
  POINT.SIZE <- 0.5
  PLOT.HEIGHT <- 7
  PLOT.WIDTH <- 7

  ggplot2::ggsave(paste0("basic_umap",".png"),
                  plot = ggplot2::ggplot(umap.layout,
                                         ggplot2::aes_string(x="umap_x",y="umap_y",color="file_var")) + #factor()
                    ggplot2::geom_point(size = POINT.SIZE) +
                    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                   panel.grid.minor = ggplot2::element_blank(),
                                   panel.background = ggplot2::element_blank(),
                                   axis.line = ggplot2::element_line(colour = "black")),
                  height = PLOT.HEIGHT, width = PLOT.WIDTH)
}









