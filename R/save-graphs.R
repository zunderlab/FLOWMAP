
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

ChangeTimes <- function(num, time.convert) {
  # ind <- match(num, as.numeric(names(time.convert)))
  convert.num <- time.convert[[num]]
  return(convert.num)
}

ConvertOrigTime <- function(graph, orig.times) {
  which.ind <- grep(pattern = "Time", list.vertex.attributes(graph))
  if (length(which.ind) > 1) {
    stop("Too many time attributes in graph!")
  }
  time.attr <- list.vertex.attributes(graph)[which.ind]
  to.fix <- get.vertex.attribute(graph, name = time.attr, index = V(graph))
  if (length(unique(to.fix)) != length(orig.times)) {
    stop("Number of original time points does not match number of distinct time points in final graphml!")
  }
  times.to.fix <- sort(unique(to.fix))
  time.convert <- list()
  for (i in 1:length(times.to.fix)) {
    time.convert[times.to.fix[i]] <- orig.times[i]
  }
  print("time.convert")
  print(time.convert)
  fixed.times <- c()
  for (i in 1:length(to.fix)) {
    fixed.times <- c(fixed.times, ChangeTimes(to.fix[i], time.convert))
  }
  fixed.graph <- set.vertex.attribute(graph, name = time.attr, index = V(graph), fixed.times)
  return(fixed.graph)
}

ConvertToPDF <- function(graphml.file, scale = NULL,
                         which.palette = "bluered") {
  pctile.color <- c(0.2, 0.98)
  node.size.scale <- 2
  min.node.size <- 12
  max.node.size <- 24
  pdf.width <- 100
  pdf.height <- 100
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
    } else {
      warning("Unknown attribute type!")
    }
  }
  colnames(all.attributes) <- attrs.colnames
  rownames(all.attributes) <- V(graph)
  # get x-y graph layout
  graph.l <- matrix(data = c(V(graph)$x, V(graph)$y), nrow = length(V(graph)$x), ncol = 2)
  # set up color scale
  if (which.palette == "jet") {
    my.palette <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    time.palette <- my.palette
  } else if (which.palette == "bluered") {
    my.palette <- colorRampPalette(c("#1500FB", "#C3C3C7", "#D10100"))
    time.palette <- colorRampPalette(rev(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")))
  } else if (which.palette == "CB") {
    my.palette <- colorRampPalette(c("#0072B2", "#C3C3C7", "#E69F00"))
    time.palette <- my.palette
  } else {
    stop("Unknown color palette!")
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
    if (name == "name") {
      next
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
        next
      }
      grad <- seq(boundary[1], boundary[2], length.out = length(color.scale))
      color <- color.scale[findInterval(attribute, grad, all.inside = TRUE)]
      color[is.na(attribute)] <- "grey"
    }
    fill.color <- color
    is.na(fill.color) <- is.na(attribute)
    frame.color <- color
    pdf(file = paste(name, ".pdf", sep = ""),
        width = pdf.width, height = pdf.height, pointsize = 12,
        bg = "transparent")
    graph.aspect <- ((max(graph.l[, 2]) - min(graph.l[, 2])) / (max(graph.l[, 1]) - min(graph.l[, 1])))
    par(mar = c(1.5, 1.5, 1.5, 1.5))
    plot(graph, layout = graph.l, vertex.shape = "circle", 
         vertex.color = fill.color, vertex.frame.color = frame.color, 
         edge.color = "#FF000000", vertex.size = vsize, edge.label = NA, 
         vertex.label = NA, edge.arrow.size = 0.25, edge.arrow.width = 1, 
         asp = graph.aspect)
    pnts <- cbind(x = c(0.85, 0.9, 0.9, 0.85), y = c(1.2, 1.2, 0.9, 0.9))
    legend.gradient(pnts = pnts, cols = my.palette(20), title = name, round(c(min(attribute), max(attribute)), 4), cex = 10)
    dev.off()
  }
  remember.attr <- setdiff(remember.attr, "id")
  if (length(remember.attr) > 0) {
    for (name in remember.attr) {
      if (grepl(name, pattern = "Time")) {
        palette.use <- time.palette
      } else {
        palette.use <- my.palette
      }
      # get attribute name and data
      attribute <- get.vertex.attribute(graph, name, index = V(graph))
      cat("options are", unique(attribute), "\n")
      num.unique <- length(unique(attribute))
      color.scale <- palette.use(num.unique)
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
      par(mar = c(1.5, 1.5, 1.5, 1.5))
      plot(graph, layout = graph.l, vertex.shape = "circle", 
           vertex.color = fill.color, vertex.frame.color = frame.color, 
           edge.color = "#FF000000", vertex.size = vsize, edge.label = NA, 
           vertex.label = NA, edge.arrow.size = 0.25, edge.arrow.width = 1, 
           asp = graph.aspect)
      if (grepl(name, pattern = "Time")) {
        legend(0.875, 1.25, legend = unique(attribute),
               fill = color.scale, cex = 10)
      } else {
        legend(0.8, 1.25, legend = unique(attribute),
               fill = color.scale, cex = 10)
      }
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

PrintSummary <- function(env = parent.frame()) {
  summary <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Variable", "Value", "Description"))
  cat("Printing summary.", "\n")
  # starting.files = c("FCS", "cluster_matrix")
  summary[(dim(summary)[1] + 1), ] <- c("mode", env$mode,
                                        "selected FLOW-MAP mode")
  summary[(dim(summary)[1] + 1), ] <- c("files", toString(env$files),
                                        "files")
  summary[(dim(summary)[1] + 1), ] <- c("var.annotate", toString(env$var.annotate),
                                        "markers included in this analysis")
  panel <- PrintPanel(env$var.annotate)
  summary[(dim(summary)[1] + 1), ] <- c("panel", toString(env$panel),
                                        "full panel including metals and corresponding marker")
  summary[(dim(summary)[1] + 1), ] <- c("var.remove", toString(env$var.remove),
                                        "removed markers")
  summary[(dim(summary)[1] + 1), ] <- c("clustering.var", toString(env$clustering.var),
                                        "markers used for clustering and distance calculation")
  summary[(dim(summary)[1] + 1), ] <- c("distance.metric", toString(env$distance.metric),
                                        "distance metric")
  summary[(dim(summary)[1] + 1), ] <- c("per", env$per,
                                        "distance for calculated density (n percent)")
  summary[(dim(summary)[1] + 1), ] <- c("minimum", env$minimum,
                                        "min number of edges")
  summary[(dim(summary)[1] + 1), ] <- c("maximum", env$maximum,
                                        "max number of edges")
  summary[(dim(summary)[1] + 1), ] <- c("subsamples", env$subsamples,
                                        "subsamples for all FCS files")
  summary[(dim(summary)[1] + 1), ] <- c("cluster.numbers", env$cluster.numbers,
                                        "number of clusters for all FCS files")
  summary[(dim(summary)[1] + 1), ] <- c("seed.X", env$seed.X,
                                        "set seed value")
  file.name <- gsub(":", ".", gsub(" ", "_", Sys.time(), fixed = TRUE), fixed = TRUE)
  file.name <- paste(file.name, "FLOW-MAPR_run_settings_summary", sep = "_")
  file.name <- paste(file.name, ".txt", sep = "")
  cat("file.name", file.name, "\n")
  write.table(summary, file = file.name, row.names = FALSE, na = "")
}
