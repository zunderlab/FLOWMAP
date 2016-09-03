
MakeOutFolder <- function(runtype) {
  name <- gsub(" ", "_", Sys.time(), fixed = TRUE)
  name <- gsub(":", ".", name, fixed = TRUE)
  output.folder <- paste(name, "_", runtype, "_run", sep = "")
  dir.create(output.folder)
  cat("output.folder is", output.folder, "\n")
  return (output.folder)
}

ConvertToGraphML <- function(output.graph, file.name) {
  cat("Converting graph to graphml file:", file.name, "\n")
  file.name <- paste(Sys.Date(), file.name, gsub(":", ".", format(Sys.time(), "%X")), sep = "_")
  file.name <- paste(file.name, ".graphml", sep = "")
  write.graph(output.graph, file.name, format = "graphml")
  return(file.name)
}

ConvertToPDF <- function(graphml.file, scale = NULL, normalize = "none", node.size.scale = 2,
                         min.node.size = 12, max.node.size = 24, pdf.width = 100, pdf.height = 100,
                         text.color = "black", edge.color = "grey", which.palette = "jet") {
  pctile.color = c(0.2, 0.98) # PCTILE_COLOR = c(0.25, 0.75)
  graph <- read.graph(graphml.file, format = "graphml")
  out.folder <- paste(basename(graphml.file), "_pdf", sep = "")
  cat("Making output folder:", out.folder, "\n")
  dir.create(out.folder)
  setwd(out.folder)
  # get matrix/table of all vertex attribute values
  attrs <- c()
  attrs.colnames <- c()
  for (a in list.vertex.attributes(graph)) {
    if (is.numeric(get.vertex.attribute(graph, a, index = V(graph)))) {
      attrs <- cbind(attrs, as.numeric(get.vertex.attribute(graph, a, index = V(graph))))
      attrs.colnames <- c(attrs.colnames, a)
    }
  }
  colnames(attrs) <- attrs.colnames
  rownames(attrs) <- V(graph)
  # get x-y graph layout
  graph.l <- matrix(data = c(V(graph)$x, V(graph)$y), nrow = length(V(graph)$x), ncol = 2)
  # for global normalization 
  boundaries <- NULL
  if (normalize == "global") {
    boundaries <- c()
    all.attrs <- c()
    for (f in graphml.files) {
      tmp.gr <- read.graph(f, format = "graphml")
      tmp.attrs <- c()
      tmp.attrs.colnames <- c()
      for (a in list.vertex.attributes(tmp.gr)) {
        if (is.numeric(get.vertex.attribute(graph, a, index = V(graph)))) {
          tmp.attrs <- cbind(tmp.attrs, as.numeric(get.vertex.attribute(tmp.gr,
                                                                        a, index = V(tmp.gr))))
          tmp.attrs.colnames <- c(tmp.attrs.colnames, a)
        }
      }
      colnames(tmp.attrs) <- tmp.attrs.colnames
      rownames(tmp.attrs) <- V(tmp.gr)
      for (c in colnames(tmp.attrs)) {
        all.attrs[[c]] <- c(all.attrs[[c]], tmp.attrs[, c])
      }
    }
    for (i in seq_along(all.attrs)) {
      boundaries[[names(all.attrs)[i]]] <- quantile(all.attrs[[i]], 
                                                    probs = pctile.color, na.rm = TRUE)
    }
  }
  # set up color scale
  if (which.palette == "jet") {
    my.palette <- colorRampPalette(c("#00007F","blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red","#7F0000"))
  }
  if (which.palette == "bluered") {
    my.palette <- colorRampPalette(c("blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red"))
  }
  color.scale <- my.palette(100)
  # set up node size
  vsize <- attrs[, "percent.total"]
  vsize <- (vsize - min(vsize, na.rm = TRUE)) / (max(vsize, na.rm = TRUE) ^ (1 / node.size.scale)) * 
    ((max.node.size) ^ 0.5/pi) + ((min.node.size) ^ (0.5 / pi))
  vsize[is.na(vsize) | (attrs[, "percent.total"] == 0)] <- (min.node.size) ^ (0.5 / pi)
  
  # print out one pdf for each attribute
  for (name in colnames(attrs)) {
    # get attribute name and data
    attr <- attrs[, name]
    # set up color boundaries
    ifelse (!is.null(scale), 
            boundary <- scale,
            ifelse (normalize == "global",
                    boundary <- boundaries[[name]],
                    boundary <- quantile(attr, probs = pctile.color, na.rm = TRUE)
            )
    )
    ifelse (length(grep("^medians|percent|cvs|Dd|Timepoint|Treatment|Treat", name)), 
            boundary <- c(min(boundary), max(boundary)),
            # boundary <- c(-max(abs(boundary)), max(abs(boundary)))
            boundary <- c(min(boundary), max(boundary))
    )
    boundary <- round(boundary, 2)
    if (boundary[1] == boundary[2]) {
      boundary <- c(boundary[1] - 1, boundary[2] + 1)
    }
    # cat("boundaries for", name, "are:", boundary, "\n")
    grad <- seq(boundary[1], boundary[2], length.out = length(color.scale))
    # cat("gradation for", name, "is:", grad, "\n")
    color <- color.scale[findInterval(attr, grad, all.inside = TRUE)]
    # cat("attribute values for", name, "is:", head(sort(attr)), "\n")
    # cat("color for", name, "is:\n")
    # print(table(color))
    color[is.na(attr) | (attrs[,"percent.total"] == 0)] <- "grey"
    if (grepl("^percenttotalratiolog$", name)) {
      color[is.na(attr) & attrs[,"percent.total"] > 0] <- tail(color.scale, 1)
    }
    fill.color <- color
    is.na(fill.color) <- is.na(attr)
    frame.color <- color
    pdf(file = paste(name, ".pdf", sep = ""),
        width = pdf.width, height = pdf.height, pointsize = 12,
        bg = "transparent")
    graph.aspect <- ((max(graph.l[, 2]) - min(graph.l[, 2])) / (max(graph.l[, 1]) - min(graph.l[, 1])))
    par(mar = c(1.5, 0, 0, 0))
    plot(graph, layout = graph.l, vertex.shape = "circle", 
         vertex.color = fill.color, vertex.frame.color = frame.color, 
         edge.color = edge.color, vertex.size = vsize, edge.label = NA, 
         vertex.label = NA, edge.arrow.size = 0.25, edge.arrow.width = 1, 
         asp = graph.aspect)
    dev.off()
  }
}


printSummary <- function(...) {
  summary <- matrix()
  cat("Printing summary.", "\n")
  if (exists("multi.folder")) {
    summary["FCS file source folder:"] <- toString(multi.folder)
  }
  if (exists("list.of.treatments")) {
    summary["Multiple treatments include:"] <- toString(list.of.treatments)
  }
  if (exists("folder")) {
    summary["FCS file source folder:"] <- toString(folder)
  }
  if (exists("var.annotate")) {
    summary["annotated variables:"] <- toString(var.annotate)
  }
  if (exists("var.remove")) {
    summary["removed variables:"] <- toString(var.remove)
  }
  if (exists("clustering.var")) {
    summary["clustering variables:"] <- toString(clustering.var)
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
  if (exists("distance.metric")) {
    summary["distance metric:"] <- toString(distance.metric)
  } 
  if (exists("subsample")) {
    summary["subsample per each FCS file:"] <- toString(subsample)
  } 
  if (exists("subsample.rand")) {
    summary["random subsample:"] <- toString(subsample.rand)
  } 
  if (exists("seed.X")) {
    summary["set seed value:"] <- toString(seed.X)
  }
  if (exists("cluster.number")) {
    summary["number of clusters per each FCS file:"] <- toString(cluster.number)
  } 
  summary <- as.data.frame(summary)
  file.name <- gsub(":", ".", gsub(" ", "_", Sys.time(), fixed = TRUE), fixed = TRUE)
  file.name <- paste(file.name, "summary", sep = "_")
  file.name <- paste(file.name, ".xls", sep = "")
  write.csv(summary, file = file.name, row.names = TRUE, na = "")
}




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

