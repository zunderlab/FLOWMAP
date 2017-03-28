
MakeOutFolder <- function(runtype) {
  name <- gsub(" ", "_", Sys.time(), fixed = TRUE)
  name <- gsub(":", ".", name, fixed = TRUE)
  output.folder <- paste(name, "_", runtype, "_run", sep = "")
  dir.create(output.folder)
  cat("output.folder is", output.folder, "\n")
  return (output.folder)
}

#' @export
ConvertToGraphML <- function(output.graph, file.name) {
  cat("Converting graph to graphml file:", file.name, "\n")
  file.name <- paste(Sys.Date(), file.name, gsub(":", ".", format(Sys.time(), "%X")), sep = "_")
  file.name <- paste(file.name, ".graphml", sep = "")
  write.graph(output.graph, file.name, format = "graphml")
  return(file.name)
}

ConvertToPDF <- function(graphml.file, scale = NULL, node.size.scale = 2,
                         min.node.size = 12, max.node.size = 24, pdf.width = 100,
                         pdf.height = 100, text.color = "black",
                         edge.color = "grey", which.palette = "jet") {
  pctile.color = c(0.2, 0.98)
  graph <- read.graph(graphml.file, format = "graphml")
  out.folder <- paste(basename(graphml.file), "_pdf", sep = "")
  cat("Making output folder:", out.folder, "\n")
  dir.create(out.folder)
  setwd(out.folder)
  # get matrix/table of all vertex attribute values
  all.attributes <- c()
  attrs.colnames <- c()
  remember.attr <- c()
  for (attribute in list.vertex.attributes(graph)) {
    if (is.numeric(get.vertex.attribute(graph, attribute, index = V(graph)))) {
      all.attributes <- cbind(all.attributes, get.vertex.attribute(graph, attribute, index = V(graph)))
      attrs.colnames <- c(attrs.colnames, attribute) 
    } else if (is.character(get.vertex.attribute(graph, attribute, index = V(graph)))) {
      remember.attr <- c(remember.attr, attribute)
    }
  }
  colnames(all.attributes) <- attrs.colnames
  rownames(all.attributes) <- V(graph)
  # get x-y graph layout
  graph.l <- matrix(data = c(V(graph)$x, V(graph)$y), nrow = length(V(graph)$x), ncol = 2)
  # set up color scale
  if (which.palette == "jet") {
    my.palette <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  }
  if (which.palette == "bluered") {
    my.palette <- colorRampPalette(c("blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red"))
  }
  color.scale <- my.palette(100)
  # set up node size
  vsize <- all.attributes[, "percent.total"]
  vsize <- (vsize - min(vsize, na.rm = TRUE)) / (max(vsize, na.rm = TRUE) ^ (1 / node.size.scale)) * 
    ((max.node.size) ^ 0.5 / pi) + ((min.node.size) ^ (0.5 / pi))
  vsize[is.na(vsize) | (all.attributes[, "percent.total"] == 0)] <- (min.node.size) ^ (0.5 / pi)
  # print out one pdf for each attribute
  for (name in colnames(all.attributes)) {
    # get attribute name and data
    attribute <- all.attributes[, name]
    if (name == "Time") {
      num.unique <- length(unique(attribute))
      color.scale <- my.palette(num.unique)
      color <- color.scale[attribute]
    } else {
      # set up color boundaries
      ifelse (!is.null(scale), 
              boundary <- scale,
              boundary <- quantile(attribute, probs = pctile.color, na.rm = TRUE)
      )
      boundary <- c(min(boundary), max(boundary))
      boundary <- round(boundary, 2)
      if (boundary[1] == boundary[2]) {
        boundary <- c(boundary[1] - 1, boundary[2] + 1)
      }
      if (Inf %in% boundary) {
        next()
      }
      grad <- seq(boundary[1], boundary[2], length.out = length(color.scale))
      color <- color.scale[findInterval(attribute, grad, all.inside = TRUE)]
      color[is.na(attribute) | (all.attributes[, "percent.total"] == 0)] <- "grey"
    }
    fill.color <- color
    is.na(fill.color) <- is.na(attribute)
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
    pnts <- cbind(x = c(0.80, 0.875, 0.875, 0.80), y = c(1.1, 1.1, 0.8, 0.8))
    legend.gradient(pnts = pnts, cols = my.palette(20), title = name, round(c(min(attribute), max(attribute)), 4), cex = 10)
    dev.off()
  }
  remember.attr <- setdiff(remember.attr, "id")
  if (length(remember.attr) > 0) {
    my.palette <- colorRampPalette(c("#00007F","blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red","#7F0000"))
    for (name in remember.attr) {
      cat("categorical attribute is", name, "\n")
      # get attribute name and data
      attribute <- get.vertex.attribute(graph, name, index = V(graph))
      cat("options are", unique(attribute), "\n")
      num.unique <- length(unique(attribute))
      color.scale <- my.palette(num.unique)
      color <- rep(NA, times = length(attribute))
      for (i in 1:length(color.scale)) {
        fix.ind <- which(attribute == unique(attribute)[i])
        color[fix.ind] <- color.scale[i]
      }
      fill.color <- color
      is.na(fill.color) <- is.na(attribute)
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
      legend(0.75, 1, legend = unique(attribute),
             fill = color.scale, cex = 10)
      dev.off()
    }
  }
}

PrintPanel <- function(var.annotate) {
  level1 <- names(unlist(var.annotate))
  level2 <- unname(unlist(var.annotate))
  arrange.ind <- order(c(seq_along(level1), seq_along(level2)))
  old.panel <- c(level1, level2)[arrange.ind]
  panel <- c()
  for (i in seq.int(from = 2, to = length(old.panel), by = 2)) {
    first <- paste("'", old.panel[(i - 1)], "'", sep = "")
    second <- paste("'", old.panel[i], "'", sep = "")
    next.i <- paste(first, second, sep = " = ")
    panel <- c(panel, next.i)
  }
  return(panel)
}

PrintSummary <- function() {
  summary <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Variable", "Value", "Description"))
  cat("Printing summary.", "\n")
  # starting.files = c("FCS", "cluster_matrix")
  
  # if (exists("mode")) {
  #   summary[(dim(summary)[1] + 1), ] <- c("mode", mode,
  #                                         "selected FLOW-MAP mode")
  # }
  if (exists("files")) {
    summary[(dim(summary)[1] + 1), ] <- c("files", toString(files),
                                          "files")
  }  
  if (exists("output.folder")) {
    summary[(dim(summary)[1] + 1), ] <- c("output.folder", toString(output.folder),
                                          "output folder")
  }  
  if (exists("final.file.name")) {
    summary[(dim(summary)[1] + 1), ] <- c("final.file.name", toString(final.file.name),
                                          "final file name")
  }  
  if (exists("fcs.file.names")) {
    summary[(dim(summary)[1] + 1), ] <- c("fcs.file.names", toString(fcs.file.names),
                                          "selected files")
  }  
  if (exists("num.files")) {
    summary[(dim(summary)[1] + 1), ] <- c("num.files", num.files,
                                          "number of selected files")
  }  
  if (exists("starting.files")) {
    summary[(dim(summary)[1] + 1), ] <- c("starting.files", toString(starting.files),
                                          "selected starting file types")
  }
  if (exists("var.annotate")) {
    summary[(dim(summary)[1] + 1), ] <- c("var.annotate", toString(var.annotate),
                                          "markers included in this analysis")
  }
  if (exists("var.annotate")) {
    panel <- PrintPanel(var.annotate)
    summary[(dim(summary)[1] + 1), ] <- c("panel", toString(panel),
                                          "full panel including metals and corresponding marker")
  }
  if (exists("var.remove")) {
    summary[(dim(summary)[1] + 1), ] <- c("var.remove", toString(var.remove),
                                          "removed markers")
  }
  if (exists("clustering.var")) {
    summary[(dim(summary)[1] + 1), ] <- c("clustering.var", toString(clustering.var),
                                          "markers used for clustering and distance calculation")
  } 
  if (exists("distance.metric")) {
    summary[(dim(summary)[1] + 1), ] <- c("distance.metric", toString(distance.metric),
                                          "distance metric")
  }
  if (exists("exclude.pctile")) {
    if (is.null(exclude.pctile)) {
      summary[(dim(summary)[1] + 1), ] <- c("exclude.pctile", "NULL",
                                            "exclude percentile for downsampling")
    } else {
      summary[(dim(summary)[1] + 1), ] <- c("exclude.pctile", exclude.pctile,
                                            "exclude percentile for downsampling")
    }
  }
  if (exists("target.pctile")) {
    if (is.null(target.pctile)) {
      summary[(dim(summary)[1] + 1), ] <- c("target.pctile", "NULL",
                                            "target percentile for downsampling")
    } else {
      summary[(dim(summary)[1] + 1), ] <- c("target.pctile", target.pctile,
                                            "target percentile for downsampling")
    }
  }
  if (exists("target.number")) {
    if (is.null(target.number)) {
      summary[(dim(summary)[1] + 1), ] <- c("target.number", "NULL",
                                            "target number for downsampling")
    } else {
      summary[(dim(summary)[1] + 1), ] <- c("target.number", target.number,
                                            "target number for downsampling")
    }
  }
  if (exists("target.percent")) {
    if (is.null(target.percent)) {
      summary[(dim(summary)[1] + 1), ] <- c("target.percent", "NULL",
                                            "target percent for downsampling")
    } else {
      summary[(dim(summary)[1] + 1), ] <- c("target.percent", target.percent,
                                            "target percent for downsampling")
    }
  }
  if (exists("per")) {
    summary[(dim(summary)[1] + 1), ] <- c("per", per,
                                          "distance for calculated density (n percent)")
  } 
  if (exists("minimum")) {
    summary[(dim(summary)[1] + 1), ] <- c("minimum", minimum,
                                          "min number of edges")
  } 
  if (exists("maximum")) {
    summary[(dim(summary)[1] + 1), ] <- c("maximum", maximum,
                                          "max number of edges")
  } 
  if (exists("subsamples")) {
    summary[(dim(summary)[1] + 1), ] <- c("subsamples", subsamples,
                                          "subsamples for all FCS files")
  } 
  if (exists("cluster.numbers")) {
    summary[(dim(summary)[1] + 1), ] <- c("cluster.numbers", cluster.numbers,
                                          "number of clusters for all FCS files")
  } 
  if (exists("seed.X")) {
    summary[(dim(summary)[1] + 1), ] <- c("seed.X", seed.X,
                                          "set seed value")
  }
  file.name <- gsub(":", ".", gsub(" ", "_", Sys.time(), fixed = TRUE), fixed = TRUE)
  file.name <- paste(file.name, "FLOW-MAPR_run_settings_summary", sep = "_")
  file.name <- paste(file.name, ".txt", sep = "")
  cat("file.name", file.name, "\n")
  write.table(summary, file = file.name, row.names = FALSE, na = "")
}
