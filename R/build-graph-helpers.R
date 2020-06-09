#############################################################################
##FUNCTIONS GRAPH BUILDING====
#############################################################################

InitializeGraph <- function(FLOWMAP.clusters) {
  # This section initializes the graph
  total.nodes <- length(FLOWMAP.clusters$full.clusters[, 1])
  # create empty graph with the right number of nodes - will fill in the edges later
  initial.graph <- graph.empty(n = total.nodes, directed = FALSE)
  return(initial.graph)
}



AnnotateGraph <- function(output.graph, FLOWMAP.clusters) {
  # This section annotates the graph
  anno.cat <- c()
  anno <- list()
  # iterate through all times and annotate
  for (f in 1:length(FLOWMAP.clusters$cluster.medians)) {
    cat("Annotating graph for file", f, "\n")
    # get medians for all parameters and counts for all clusters
    counts <- FLOWMAP.clusters$cluster.counts[[f]]$Counts
    anno$count <- counts
    total.cell <- sum(counts)
    anno$percent.total <- data.frame(percent.total = c(counts / total.cell))
    anno$medians  <- FLOWMAP.clusters$cluster.medians[[f]]
    # add time information column
    time.matrix <- matrix(f, nrow = length(anno$count))
    colnames(time.matrix) <- c("Timepoint")
    anno$medians <- cbind(anno$medians, time.matrix)
    # add median and percent values
    for (col in c("medians", "percent.total")) {
      anno.cat[[col]] <- rbind(anno.cat[[col]], anno[[col]])
    }
  }
  # combine anno.cat matrices
  output.anno <- cbind(anno.cat[[1]], anno.cat[[2]])
  print("1")
  for (c in colnames(output.anno)) {
    output.graph <- set.vertex.attribute(output.graph, c,
                                         index = as.numeric(rownames(output.anno)),
                                         value = output.anno[, c])
  }
  print("2")
  # add name attribute
  V(output.graph)$name <- 1:length(V(output.graph))
  return(output.graph)
}

AnnotateSpecialGraph <- function(output.graph, FLOWMAP.clusters,
                                 label.key.special) {
  # This section annotates the graph
  anno.cat <- c()
  anno <- list()
  # iterate through all times and annotate
  for (f in 1:length(FLOWMAP.clusters$cluster.medians)) {
    cat("Annotating graph for file", f, "\n")
    # get medians for all parameters and counts for all clusters
    counts <- FLOWMAP.clusters$cluster.counts[[f]]$Counts
    anno$count <- counts
    total.cell <- sum(counts)
    anno$percent.total <- data.frame(percent.total = c(counts / total.cell))
    anno$medians  <- FLOWMAP.clusters$cluster.medians[[f]]
    # add time information column
    time.matrix <- matrix(f, nrow = length(anno$count))
    colnames(time.matrix) <- c("Timepoint")
    anno$medians <- cbind(anno$medians, time.matrix)
    # add median and percent values
    for (col in c("medians", "percent.total")) {
      anno.cat[[col]] <- rbind(anno.cat[[col]], anno[[col]])
    }
  }
  # combine anno.cat matrices
  output.anno <- cbind(anno.cat[[1]], anno.cat[[2]])
  print("label.key.special")
  print(label.key.special)
  output.anno <- ConvertCharacterLabelSpecial(output.anno, label.key.special)
  for (c in colnames(output.anno)) {
    if (c == "Condition") {
      output.graph <- set.vertex.attribute(output.graph, c,
                                           index = as.numeric(1:dim(output.anno)[1]),
                                           value = as.character(output.anno[, c]))
    } else {
      output.graph <- set.vertex.attribute(output.graph, c,
                                           index = as.numeric(1:dim(output.anno)[1]),
                                           value = output.anno[, c])
    }
  }
  # add name attribute
  V(output.graph)$name <- 1:length(V(output.graph))
  return(output.graph)
}
