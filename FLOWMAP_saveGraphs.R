library(igraph)

convertToGraphML <- function(output_graph, output_folder, file_name, ...) {
  cat("Converting graph to graphml file:", file_name, "\n")
  file_name <- paste(Sys.Date(), file_name, gsub(":", ".", format(Sys.time(), "%X")), sep = "_")
  write.graph(output_graph, paste(output_folder, "/", file_name,".graphml", sep = ""),
              format = "graphml")
}


convertToPDF <- function(in_folder, file_pattern = ".graphml", out_folder = FALSE,
                         scale = NULL, normalize = "none", node_size_scale = 2,
                         min_node_size = 12, max_node_size = 24, pdf_width = 100, pdf_height = 100,
                         text_color = "black", edge_color = "grey", PALETTE = "jet", listOfTreatments) {
  #                          treatInvisible = TRUE, timeInvisible = TRUE, visibleChannels
  # visibleChannels should be a vector of strings that are channels to be colored and displayed
  # with different treatments/times invisible
  PCTILE_COLOR = c(0.2, 0.98)
  #   PCTILE_COLOR = c(0.25, 0.75)
  BARE = FALSE
  TEST = FALSE
  if (file_pattern != "DEFAULT") {
    graphml_file <- list.files(path = in_folder, pattern = file_pattern, full.names = TRUE)
    graph <- read.graph(graphml_file, format = "graphml")
  }
  else {
    graphml_file <- list.files(path = in_folder, pattern = "xy", full.names = TRUE)
    notind <- grep(".graphmlpdf", graphml_file)
    ind <- setdiff((1:length(graphml_file)), notind)
    print(graphml_file[1])
    graph <- read.graph(graphml_file[1], format = "graphml")
  }
  if (out_folder == FALSE) {
    out_folder <- paste(in_folder, "/", basename(graphml_file), "pdf/", sep = "")
    cat("Making output folder:", out_folder, "\n")
    dir.create(out_folder)
  }
  else {
    cat("Finding output folder:", out_folder, "\n")
    setwd(out_folder)
  }
  # get matrix/table of all vertex attribute values
  attrs <- c()
  attrs_colnames <- c()
  if (TEST) {
    a <- "percenttotal"
    attrs <- matrix(data = as.numeric(get.vertex.attribute(graph, a, index = V(graph))),
                    ncol = 1)
    attrs_colnames <- c(attrs_colnames, a)
  } else {
    for (a in list.vertex.attributes(graph)) {
      if (is.numeric(get.vertex.attribute(graph, a, index = V(graph)))) {
        attrs <- cbind(attrs, as.numeric(get.vertex.attribute(graph, a, index = V(graph))))
        attrs_colnames <- c(attrs_colnames, a)
      }
    }
  }
  colnames(attrs) <- attrs_colnames
  rownames(attrs) <- V(graph)
  # get x-y graph layout
  graph_l <- matrix(data = c(V(graph)$x, V(graph)$y), nrow = length(V(graph)$x), ncol = 2)
  # for global normalization 
  boundaries <- NULL
  if (normalize == "global") {
    boundaries <- c()
    all_attrs <- c()
    for (f in graphml_files) {
      tmp_gr <- read.graph(f, format = "graphml")
      tmp_attrs <- c()
      tmp_attrs_colnames <- c()
      for (a in list.vertex.attributes(tmp_gr)) {
        if (is.numeric(get.vertex.attribute(graph, a, index = V(graph)))) {
          tmp_attrs <- cbind(tmp_attrs, as.numeric(get.vertex.attribute(tmp_gr,
                                                                        a, index = V(tmp_gr))))
          tmp_attrs_colnames <- c(tmp_attrs_colnames, a)
        }
      }
      colnames(tmp_attrs) <- tmp_attrs_colnames
      rownames(tmp_attrs) <- V(tmp_gr)
      for (c in colnames(tmp_attrs)) {
        all_attrs[[c]] <- c(all_attrs[[c]], tmp_attrs[, c])
      }
    }
    for (i in seq_along(all_attrs)) {
      boundaries[[names(all_attrs)[i]]] <- quantile(all_attrs[[i]], 
                                                    probs = PCTILE_COLOR, na.rm = TRUE)
    }
  }
  # set up color scale
  if (PALETTE == "jet")
    palette <- colorRampPalette(c("#00007F","blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red","#7F0000"))
  if (PALETTE == "bluered")
    palette <- colorRampPalette(c("blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red"))
  colorSCALE <- palette(100)
  #   cat("color scale is", colorSCALE, "\n")
  # set up node size
  vsize <- attrs[, "percenttotal"]
  vsize <- (vsize - min(vsize, na.rm = TRUE)) / (max(vsize, na.rm = TRUE) ^ (1/node_size_scale)) * 
    ((max_node_size) ^ 0.5/pi) + ((min_node_size) ^ 0.5/pi)
  vsize[is.na(vsize) | (attrs[, "percenttotal"] == 0)] <- (min_node_size) ^ 0.5/pi
  
  # print out one pdf for each attribute
  for (name in colnames(attrs)) {
    #get attribute name and data
    attr <- attrs[, name]
    # set up color boundaries
    ifelse (!is.null(scale), 
            boundary <- scale,
            ifelse (normalize == "global",
                    boundary <- boundaries[[name]],
                    boundary <- quantile(attr, probs = PCTILE_COLOR, na.rm = TRUE)
            )
    )
    ifelse (length(grep("^medians|percent|cvs|Dd|Timepoint|Treatment|Treat", name)), 
            boundary <- c(min(boundary), max(boundary)),
            #             boundary <- c(-max(abs(boundary)), max(abs(boundary)))
            boundary <- c(min(boundary), max(boundary))
    )
    boundary <- round(boundary, 2)
    if (boundary[1] == boundary[2]) {
      #       print("boop")
      boundary <- c(boundary[1] - 1, boundary[2] + 1)
    }
    #     cat("boundaries for", name, "are:", boundary, "\n")
    grad <- seq(boundary[1], boundary[2], length.out = length(colorSCALE))
    #     cat("gradation for", name, "is:", grad, "\n")
    color <- colorSCALE[findInterval(attr, grad, all.inside = TRUE)]
    #     cat("attribute values for", name, "is:", head(sort(attr)), "\n")
    #     cat("color for", name, "is:\n")
    #     print(table(color))
    color[is.na(attr) | (attrs[,"percenttotal"] == 0)] <- "grey"
    if (grepl("^percenttotalratiolog$", name)) {
      color[is.na(attr) & attrs[,"percenttotal"] > 0] <- tail(colorSCALE, 1)
    }
    fill_color <- color
    is.na(fill_color) <- is.na(attr)
    frame_color <- color
    pdf(file = paste(out_folder, "/", name, ".pdf", sep = ""),
        width = pdf_width, height = pdf_height, pointsize = 12,
        bg = "transparent")
    graph_aspect <- ((max(graph_l[, 2]) - min(graph_l[, 2]))/(max(graph_l[, 1]) - min(graph_l[, 1])))
    par(mar = c(1.5, 0, 0, 0))
    plot(graph, layout = graph_l, vertex.shape = "circle", 
         vertex.color = fill_color, vertex.frame.color = frame_color, 
         edge.color = edge_color, vertex.size = vsize, edge.label = NA, 
         vertex.label = NA, edge.arrow.size = 0.25, edge.arrow.width = 1, 
         asp = graph_aspect)
    dev.off()
    
    #     if (treatInvisible | timeInvisible & name %in% visibleChannels) {
    #       #       out_folder <- paste(out_folder, "/invisible_pdf/", sep = "")
    #       #       cat("Making output folder for pdfs with treat or time invisible:", out_folder, "\n")
    #       #       dir.create(out_folder)
    #       color_copy <- color
    #       
    #       if (treatInvisible) {
    #         color <- color_copy
    #         treatments <- unique(attrs[, "Treatment"])
    #         for (treat in treatments) {
    #           treatNum <- which(listOfTreatments == treat)
    #           ind <- which(as.numeric(get.vertex.attribute(graph, "Treatment", index = V(graph))) == treatNum)
    #           allind <- seq(1:length(V(graph)))
    #           ind <- length(setdiff(allind, ind))
    #           color[ind] <- "#FF000000"
    #           fill_color <- color
    #           is.na(fill_color) <- is.na(attr)
    #           frame_color <- color
    #           pdf(file = paste(out_folder, "/", name, "with treat", treat, "invisble.pdf", sep = ""),
    #               width = pdf_width, height = pdf_height, pointsize = 12,
    #               bg = "transparent")
    #           graph_aspect <- ((max(graph_l[, 2]) - min(graph_l[, 2]))/(max(graph_l[, 1]) - min(graph_l[, 1])))
    #           par(mar = c(1.5, 0, 0, 0))
    #           plot(graph, layout = graph_l, vertex.shape = "circle", 
    #                vertex.color = fill_color, vertex.frame.color = frame_color, 
    #                edge.color = "#FF000000", vertex.size = vsize, edge.label = NA, 
    #                vertex.label = NA, edge.arrow.size = 0.25, edge.arrow.width = 1, 
    #                asp = graph_aspect)
    #           dev.off()
    #         }
    #       }
    #       
    #       if (timeInvisible) {
    #         color <- color_copy
    #         time <- unique(attrs[, "Timepoint"])
    #         for (t in time) {
    #           ind <- which(as.numeric(get.vertex.attribute(graph, "Timepoint", index = V(graph))) == t)
    #           allind <- seq(1:length(V(graph)))
    #           ind <- length(setdiff(allind, ind))
    #           color[ind] <- "#FF000000"
    #           fill_color <- color
    #           is.na(fill_color) <- is.na(attr)
    #           frame_color <- color
    #           pdf(file = paste(out_folder, "/", name, "with time", t, "invisble.pdf", sep = ""),
    #               width = pdf_width, height = pdf_height, pointsize = 12,
    #               bg = "transparent")
    #           graph_aspect <- ((max(graph_l[, 2]) - min(graph_l[, 2]))/(max(graph_l[, 1]) - min(graph_l[, 1])))
    #           par(mar = c(1.5, 0, 0, 0))
    #           plot(graph, layout = graph_l, vertex.shape = "circle", 
    #                vertex.color = fill_color, vertex.frame.color = frame_color, 
    #                edge.color = "#FF000000", vertex.size = vsize, edge.label = NA, 
    #                vertex.label = NA, edge.arrow.size = 0.25, edge.arrow.width = 1, 
    #                asp = graph_aspect)
    #           dev.off()
    #         }
    #       }
    #     }
  }
}



#   
#   if (treatInvisible) {
#     #   # make all nodes for invisible treatment size = 0
#     #   # set up node size
#     #   vsize <- attrs[, "percenttotal"]
#     #   vsize <- (vsize - min(vsize, na.rm = TRUE)) / (max(vsize, na.rm = TRUE) ^ (1/node_size_scale)) * 
#     #     ((max_node_size) ^ 0.5/pi) + ((min_node_size) ^ 0.5/pi)
#     #   vsize[is.na(vsize) | (attrs[, "percenttotal"] == 0)] <- (min_node_size) ^ 0.5/pi
#     
#     
#     # make edges invisible
#     # invisible color ? "#FF000000"
#     
#     # output visible channels with one of each treat
#     # visible or one of each time visible, permutations thereof
#     for (name in colnames(attrs)) {
#       #get attribute name and data
#       attr <- attrs[, name]
#       # set up color boundaries
#       ifelse (!is.null(scale), 
#               boundary <- scale,
#               ifelse (normalize == "global",
#                       boundary <- boundaries[[name]],
#                       boundary <- quantile(attr, probs = PCTILE_COLOR, na.rm = TRUE)
#               )
#       )
#       ifelse (length(grep("^medians|percent|cvs|Dd|Timepoint|Treatment|Treat", name)), 
#               boundary <- c(min(boundary), max(boundary)),
#               #             boundary <- c(-max(abs(boundary)), max(abs(boundary)))
#               boundary <- c(min(boundary), max(boundary))
#       )
#       boundary <- round(boundary, 2)
#       if (boundary[1] == boundary[2]) {
#         #       print("boop")
#         boundary <- c(boundary[1] - 1, boundary[2] + 1)
#       }
#       #     cat("boundaries for", name, "are:", boundary, "\n")
#       grad <- seq(boundary[1], boundary[2], length.out = length(colorSCALE))
#       #     cat("gradation for", name, "is:", grad, "\n")
#       color <- colorSCALE[findInterval(attr, grad, all.inside = TRUE)]
#       #     cat("attribute values for", name, "is:", head(sort(attr)), "\n")
#       #     cat("color for", name, "is:\n")
#       #     print(table(color))
#       color[is.na(attr) | (attrs[,"percenttotal"] == 0)] <- "grey"
#       if (grepl("^percenttotalratiolog$", name)) {
#         color[is.na(attr) & attrs[,"percenttotal"] > 0] <- tail(colorSCALE, 1)
#       }
#       
#       treatments <- attrs[, "Treatment"]
#       
#       for (treat in treatments) {
#       
#       # iterate through all vertex/nodes and set nodes that are in the
#       # excluded time/treatment to color invisible
#       color <- "#FF000000"
#       
#       fill_color <- color
#       is.na(fill_color) <- is.na(attr)
#       frame_color <- color
#       pdf(file = paste(out_folder, "/", name, ".pdf", sep = ""),
#           width = pdf_width, height = pdf_height, pointsize = 12,
#           bg = "transparent")
#       graph_aspect <- ((max(graph_l[, 2]) - min(graph_l[, 2]))/(max(graph_l[, 1]) - min(graph_l[, 1])))
#       par(mar = c(1.5, 0, 0, 0))
#       plot(graph, layout = graph_l, vertex.shape = "circle", 
#            vertex.color = fill_color, vertex.frame.color = frame_color, 
#            edge.color = "#FF000000", vertex.size = vsize, edge.label = NA, 
#            vertex.label = NA, edge.arrow.size = 0.25, edge.arrow.width = 1, 
#            asp = graph_aspect)
#       dev.off()
#       }
#     }
#   }
#   
#   if (timeInvisible) {
#     #   # make all nodes for invisible time size = 0
#     #   # set up node size
#     #   vsize <- attrs[, "percenttotal"]
#     #   vsize <- (vsize - min(vsize, na.rm = TRUE)) / (max(vsize, na.rm = TRUE) ^ (1/node_size_scale)) * 
#     #     ((max_node_size) ^ 0.5/pi) + ((min_node_size) ^ 0.5/pi)
#     #   vsize[is.na(vsize) | (attrs[, "percenttotal"] == 0)] <- (min_node_size) ^ 0.5/pi
#     
#     
#     # make edges invisible
#     # invisible color ? "#FF000000"
#     
#     # output visible channels with one of each treat
#     # visible or one of each time visible, permutations thereof
#     for (name in colnames(attrs)) {
#       #get attribute name and data
#       attr <- attrs[, name]
#       # set up color boundaries
#       ifelse (!is.null(scale), 
#               boundary <- scale,
#               ifelse (normalize == "global",
#                       boundary <- boundaries[[name]],
#                       boundary <- quantile(attr, probs = PCTILE_COLOR, na.rm = TRUE)
#               )
#       )
#       ifelse (length(grep("^medians|percent|cvs|Dd|Timepoint|Treatment|Treat", name)), 
#               boundary <- c(min(boundary), max(boundary)),
#               #             boundary <- c(-max(abs(boundary)), max(abs(boundary)))
#               boundary <- c(min(boundary), max(boundary))
#       )
#       boundary <- round(boundary, 2)
#       if (boundary[1] == boundary[2]) {
#         #       print("boop")
#         boundary <- c(boundary[1] - 1, boundary[2] + 1)
#       }
#       #     cat("boundaries for", name, "are:", boundary, "\n")
#       grad <- seq(boundary[1], boundary[2], length.out = length(colorSCALE))
#       #     cat("gradation for", name, "is:", grad, "\n")
#       color <- colorSCALE[findInterval(attr, grad, all.inside = TRUE)]
#       #     cat("attribute values for", name, "is:", head(sort(attr)), "\n")
#       #     cat("color for", name, "is:\n")
#       #     print(table(color))
#       color[is.na(attr) | (attrs[,"percenttotal"] == 0)] <- "grey"
#       if (grepl("^percenttotalratiolog$", name)) {
#         color[is.na(attr) & attrs[,"percenttotal"] > 0] <- tail(colorSCALE, 1)
#       }
#       
#       
#       # iterate through all vertex/nodes and set nodes that are in the
#       # excluded time/treatment to color invisible
#       color <- "#FF000000"
#       
#       fill_color <- color
#       is.na(fill_color) <- is.na(attr)
#       frame_color <- color
#       pdf(file = paste(out_folder, "/", name, ".pdf", sep = ""),
#           width = pdf_width, height = pdf_height, pointsize = 12,
#           bg = "transparent")
#       graph_aspect <- ((max(graph_l[, 2]) - min(graph_l[, 2]))/(max(graph_l[, 1]) - min(graph_l[, 1])))
#       par(mar = c(1.5, 0, 0, 0))
#       plot(graph, layout = graph_l, vertex.shape = "circle", 
#            vertex.color = fill_color, vertex.frame.color = frame_color, 
#            edge.color = "#FF000000", vertex.size = vsize, edge.label = NA, 
#            vertex.label = NA, edge.arrow.size = 0.25, edge.arrow.width = 1, 
#            asp = graph_aspect)
#       dev.off()
#     }
#   }


printSummary <- function(out_folder, ...) {
  summary <- matrix()
  print("Printing summary.")
  if (exists("MULTI_FOLDER")) {
    summary["FCS file source folder:"] <- toString(MULTI_FOLDER)
  }
  if (exists("listOfTreatments")) {
    summary["Multiple treatments include:"] <- toString(listOfTreatments)
  }
  if (exists("FOLDER")) {
    summary["FCS file source folder:"] <- toString(FOLDER)
  }
  if (exists("VAR_ANNOTATE")) {
    summary["annotated variables:"] <- toString(VAR_ANNOTATE)
  }
  if (exists("VAR_REMOVE")) {
    summary["removed variables:"] <- toString(VAR_REMOVE)
  }
  if (exists("VAR_UPSAMPLE")) {
    summary["upsampled variables:"] <- toString(VAR_UPSAMPLE)
  }  
  if (exists("CLUSTERING_VAR")) {
    summary["clustering variables:"] <- toString(CLUSTERING_VAR)
  } 
  if (exists("per")) {
    summary["distance for calculated density (n percent):"] <- toString(per)
  } 
  if (exists("minimum")) {
    summary["min number of edges:"] <- toString(minimum)
  } 
  if (exists("maximum")) {
    summary["max number of edges:"] <- toString(maximum)
  } 
  if (exists("distance_metric")) {
    summary["distance metric:"] <- toString(distance_metric)
  } 
  if (exists("SUBSAMPLE")) {
    summary["subsample per each FCS file:"] <- toString(SUBSAMPLE)
  } 
  if (exists("subsampleRand")) {
    summary["random subsample:"] <- toString(subsampleRand)
  } 
  if (exists("seedX")) {
    summary["set seed value:"] <- toString(seedX)
  }
  if (exists("CLUSTNUM")) {
    summary["number of clusters per each FCS file:"] <- toString(CLUSTNUM)
  } 
  summary <- as.data.frame(summary)
  file_name <- gsub(":", ".", gsub(" ", "_", Sys.time(), fixed = TRUE), fixed = TRUE)
  file_name <- paste(file_name, "summary", sep = "_")
  file_name <- paste(file_name, ".xls", sep = "")
  out_file <- paste(out_folder, file_name, sep = "/") 
  write.csv(summary, file = out_file, row.names = TRUE, na = "", ...)
}
