
# ForceDirectedXY <- function(graph) {
#   V(graph)$size <- rep(20, vcount(graph))
#   force.graph1 <- scaffold::layout.forceatlas2(graph, iter = 10000,
#                                                 stopping_tolerance = 0.001,
#                                                 prevent.overlap = FALSE)
#   graph.with.xy <- graph
#   V(graph.with.xy)$x <- force.graph1$lay[, 1]
#   V(graph.with.xy)$y <- force.graph1$lay[, 2]
#   force.graph2 <- scaffold::layout.forceatlas2(graph.with.xy, iter = 1000,
#                                                 stopping_tolerance = 0.001,
#                                                 prevent.overlap = TRUE)
#   V(graph.with.xy)$x <- force.graph2$lay[, 1]
#   V(graph.with.xy)$y <- force.graph2$lay[, 2]
#   return(graph.with.xy)
# }

RunForceDirectedLayout <- function(mode, file.name, graph, orig.times=NULL, which.palette=NULL) {
  global.graph.pre.ml <<- graph
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
  
  if (savePDFs) {
    cat("Printing pdfs.", "\n")
    ConvertToPDF(graphml.file = fixed.file, which.palette = which.palette)
  }
  return(graph.xy)
}


#' Generate UMAP layout for FLOWMAP graph
#' @import ggplot2
#' @import ggfortify
RunUMAPlayout <- function(graph, knn.in, file.clusters, clustering.var, file.name=file.name, 
                          umap_n_neighbors, k, umap_n_components, mode) {

  global.file.clusters <<- file.clusters
  # #Set up UMAP settings
  # umap.settings <- umap::umap.defaults
  # umap.settings$verbose <- TRUE
  # umap.settings$n_neighbors <- umap_n_neighbors
  # umap.settings$n_components <- umap_n_components
  
  ## Run UMAP
  cat("running UMAP\n")
  
  self_col_ind <- c(1:nrow(knn.in$indexes))
  self_col_dist <- rep(0.0, nrow(knn.in$indexes))
  knn.in$indexes <- data.frame(cbind(self_col_ind,data.frame(knn.in$indexes)))
  knn.in$indexes <- as.matrix(knn.in$indexes)
  knn.in$distances <- data.frame(cbind(self_col_dist,knn.in$distances))
  knn.in$distances <- as.matrix(knn.in$distances)
  knn <- list()
  if (umap_n_neighbors > k) {
    umap_n_neighbors <- k
    print("umap_n_neighboors must be <= k, has been changed to k")
  }
  umap_n_neighbors <- umap_n_neighbors+1 #because uwot requires self connection to be in knn matrix
  knn[['idx']] <- as.matrix(knn.in$indexes[,1:umap_n_neighbors]) #make sure ordered
  knn[['dist']] <- as.matrix(knn.in$distances[,1:umap_n_neighbors]) #make sure ordered
  umap.out <- uwot::umap(file.clusters$full.clusters, ret_nn = TRUE, nn_method = knn, verbose=TRUE, n_components=umap_n_components)
  global.umap.out <<- umap.out
  cat("ran UMAP, outputting files\n")

  ## Add file var to umap layout
  #TODO identify better way to track timepoint/other potentially needed info?
  # Idea: timepoint attribute... need to figure out how exactly it works though...
  diff.file.var <- c()
  for (i in 1:length(file.clusters$table.lengths)) {
    diff.file.var <- append(diff.file.var, rep(i, file.clusters$table.lengths[i]))
  }
  
  ## Outputs for 2D UMAP ====
  if (umap_n_components == 2) {
    ## Save UMAP layout and knn graph
    #data.table::fwrite(umap.out$knn$indexes,file = "UMAP_knn_indexes.csv",row.names = FALSE,col.names = FALSE,sep = ",")
    data.table::fwrite(umap.out$nn$precomputed$idx,file = "UMAP_knn_indexes.csv",row.names = FALSE,col.names = FALSE,sep = ",")
    #data.table::fwrite(umap.out$knn$distances,file = "UMAP_knn_distances.csv",row.names = FALSE,col.names = FALSE,sep = ",")
    data.table::fwrite(umap.out$nn$precomputed$dist,file = "UMAP_knn_distances.csv",row.names = FALSE,col.names = FALSE,sep = ",")
    #colnames(umap.out$layout) <- c("umap_x","umap_y")
    #data.table::fwrite(umap.out$layout,file = "UMAP_layout.csv",row.names = FALSE,col.names = TRUE, sep = ",") #umap package
    colnames(umap.out$embedding) <- c("umap_x","umap_y")
    data.table::fwrite(umap.out$embedding,file = "UMAP_layout.csv",row.names = FALSE,col.names = TRUE, sep = ",") #uwot package
    
    print("Generating 2D layouts colored by cluster") 
    # Add file.var to layout to color by file number
    #umap.layout <- data.frame(cbind(umap.out$layout, diff.file.var)) #umap package
    umap.layout <- data.frame(cbind(umap.out$embedding, diff.file.var)) #uwot package
    colnames(umap.layout) <- c("umap_x","umap_y","file_var")
    umap.layout$file_var <- as.factor(umap.layout$file_var)
    # Extract cluster cell number percent total attribute to assign size of points in layout plot
    cluster.size <- igraph::get.vertex.attribute(graph, "percent.total", index = V(graph))
    cluster.size <- (cluster.size / max(cluster.size))*(4 - 1) + 1 ## normalize to smaller range to make reasonably-sized points
    umap.layout$cluster_size <- cluster.size
    
    global.umap.layout <<- umap.layout
    
    #Plot UMAP
    PLOT.HEIGHT <- 7
    PLOT.WIDTH <- 7

    # Plot colored by time point
    ggplot2::ggsave(paste0("umap_layout_TIME",".png"),
                    plot = ggplot2::ggplot(umap.layout,
                                           ggplot2::aes_string(x="umap_x",y="umap_y",color="file_var")) + #factor()
                      ggplot2::geom_point(size = umap.layout$cluster_size*0.2, alpha = 0.4) +
                      ggplot2::geom_jitter() +
                      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                     panel.grid.minor = ggplot2::element_blank(),
                                     panel.background = ggplot2::element_blank(),
                                     axis.line = ggplot2::element_line(colour = "black")),
                    height = PLOT.HEIGHT, width = PLOT.WIDTH)

    if (mode %in% c('multi', 'static-multi')) {
      # Extract "Condition"
      condition.chr.id <- igraph::get.vertex.attribute(graph, "Condition", index = V(graph))
      umap.layout$condition.chr.id <- as.factor(condition.chr.id)
      # Plot colored by condition
      ggplot2::ggsave(paste0("umap_layout_CONDITION",".png"),
                      plot = ggplot2::ggplot(umap.layout,
                                             ggplot2::aes_string(x="umap_x",y="umap_y",color="condition.chr.id")) + #factor()
                        ggplot2::geom_point(size = umap.layout$cluster_size*0.2, alpha = 0.4) +
                        ggplot2::geom_jitter() +
                        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                       panel.grid.minor = ggplot2::element_blank(),
                                       panel.background = ggplot2::element_blank(),
                                       axis.line = ggplot2::element_line(colour = "black")),
                      height = PLOT.HEIGHT, width = PLOT.WIDTH)
    }
    
    # Generate plots for clustering markers expression
    #First, make data frame out of cluster.medians
    print("Generating 2D layouts colored by marker expression")
    #cluster.medians <- dplyr::bind_rows(file.clusters$cluster.medians)
    clusters.exprs <- file.clusters$full.clusters[,which(colnames(file.clusters$full.clusters) != 'Condition')]
    #clusters.exprs <- igraph::get.vertex.attribute(graph, index = V(graph))
    umap.layout <- data.frame(cbind(umap.layout, clusters.exprs))
    umap.layout <- reshape2::melt(umap.layout, id.vars = c("umap_x","umap_y","file_var", "cluster_size"), measure.vars = colnames(data.frame(clusters.exprs)))
    
    basic_plot <- function(umap.layout) {
      ggplot2::ggplot(umap.layout,ggplot2::aes_string(x="umap_x",y="umap_y",color="value")) + #factor()
        ggplot2::geom_point(size = umap.layout$cluster_size/10, alpha = 0.5) +
        ggplot2::scale_color_viridis_c() +
        ggplot2::ggtitle(umap.layout$variable) +
        #ggplot2::labs(color = paste0(variable,"\n")) + #umap.layout$
        #ggplot2::facet_wrap(~variable, scales = "free") +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(colour = "black"))
    } 

    # Open file for writing
    pdf("umap_feature_heatmaps.pdf", width = 20, height = 2*round(length(clusters.exprs)/4)) # Open a new pdf file
    # Write the grid.arrange in the file
    do.call(gridExtra::grid.arrange, 
            args=list(grobs=by(umap.layout, umap.layout$variable, basic_plot), 
                      nrow = round((length(clusters.exprs)/5)+0.49)))
    dev.off() # Close the file
    # #make individual plots per marker
    # for (marker in colnames(clusters.exprs)) {
    #   print(marker)
    #   umap.layout.temp <- data.frame(cbind(umap.layout, clusters.exprs[marker]))
    #   colnames(umap.layout.temp) <- c("umap_x","umap_y","file_var","cluster_size","marker")
    #   print(umap.layout.temp[1:10,4])
    #   ggplot2::ggsave(paste0("basic_umap_",marker,".png"),
    #                   plot = ggplot2::ggplot(umap.layout.temp,
    #                                          ggplot2::aes_string(x="umap_x",y="umap_y",color="marker")) + #factor()
    #                     ggplot2::geom_point(size = umap.layout.temp$cluster_size, alpha = 0.8) +
    #                     ggplot2::geom_jitter() +
    #                     ggplot2::scale_color_viridis_c() +
    #                     ggplot2::labs(color = paste0(marker,"\n")) +
    #                     ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
    #                                    panel.grid.minor = ggplot2::element_blank(),
    #                                    panel.background = ggplot2::element_blank(),
    #                                    axis.line = ggplot2::element_line(colour = "black")),
    #                   height = PLOT.HEIGHT, width = PLOT.WIDTH)
    # }
  }
  
  ## Outputs for 3D UMAP ====
  else if (umap_n_components == 3) {
    
    ## Save UMAP layout and knn graph
    #data.table::fwrite(umap.out$knn$indexes,file = "UMAP_knn_indexes.csv",row.names = FALSE,col.names = FALSE,sep = ",")
    data.table::fwrite(umap.out$nn$precomputed$idx,file = "UMAP_knn_indexes_3D.csv",row.names = FALSE,col.names = FALSE,sep = ",")
    #data.table::fwrite(umap.out$knn$distances,file = "UMAP_knn_distances.csv",row.names = FALSE,col.names = FALSE,sep = ",")
    data.table::fwrite(umap.out$nn$precomputed$dist,file = "UMAP_knn_distances_3D.csv",row.names = FALSE,col.names = FALSE,sep = ",")
    #colnames(umap.out$layout) <- c("umap_x","umap_y")
    #data.table::fwrite(umap.out$layout,file = "UMAP_layout.csv",row.names = FALSE,col.names = TRUE, sep = ",") #umap package
    colnames(umap.out$embedding) <- c("umap_x","umap_y", "umap_z")
    data.table::fwrite(umap.out$embedding,file = "UMAP_layout_3D.csv",row.names = FALSE,col.names = TRUE, sep = ",") #uwot package
    
    print("Generating 3D layouts") 
    # Add file.var to layout to color by file number
    umap.layout <- data.frame(cbind(umap.out$embedding, diff.file.var)) #uwot package
    colnames(umap.layout) <- c("umap_x","umap_y","umap_z","file_var")
    umap.layout$file_var <- as.factor(umap.layout$file_var)
    
    # Extract cluster cell number percent total attribute to assign size of points in layout plot
    cluster.size <- igraph::get.vertex.attribute(graph, "percent.total", index = V(graph))
    cluster.size <- (cluster.size / max(cluster.size))*(3) + 3 ## normalize to smaller range to make reasonably-sized points
    umap.layout$cluster_size <- cluster.size
    
    global.umap.layout <<- umap.layout
    
    print("Generating 3D layouts colored by cluster")
    ## 3D PLOTTING
    # Plot your data, in this example my Seurat object had 21 clusters (0-20)
    p <- plotly::plot_ly(data = umap.layout, 
                         x = ~umap_x, y = ~umap_y, z = ~umap_z, 
                         color = ~file_var, 
                         opacity = 0.5,
                         colors = c("lightseagreen",
                                    "gray50",
                                    "darkgreen",
                                    "red4",
                                    "red",
                                    "turquoise4",
                                    "black",
                                    "yellow4",
                                    "royalblue1",
                                    "lightcyan3",
                                    "peachpuff3",
                                    "khaki3",
                                    "gray20",
                                    "orange2",
                                    "royalblue4",
                                    "yellow3",
                                    "gray80",
                                    "darkorchid1",
                                    "lawngreen",
                                    "plum2",
                                    "darkmagenta"),
                         type = "scatter3d", 
                         mode = "markers", 
                         marker = list(size = cluster.size, width=1))
    htmlwidgets::saveWidget(plotly::as.widget(p), paste0("umap_layout_3d_plotly",".html"))
  }
}










