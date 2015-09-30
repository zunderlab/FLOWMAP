library(igraph)
library(RForceAtlas2)
library(proxy)


forceDirectedXY <- function(graph) {
  force.graph1 <- layout.forceatlas2(graph, iter = 10000, stopping_tolerance = 0.001,
                                    prevent.overlap = FALSE)
  graph_with_xy <- graph
  cat("x:", head(force.graph1$lay[, 1]), "\n")
  cat("y:", head(force.graph1$lay[, 2]), "\n")
  V(graph_with_xy)$x <- force.graph1$lay[, 1]
  V(graph_with_xy)$y <- force.graph1$lay[, 2]
  force.graph2 <- layout.forceatlas2(graph_with_xy, iter = 1000, stopping_tolerance = 0.001,
                                    prevent.overlap = TRUE)
  cat("x:", head(force.graph2$lay[, 1]), "\n")
  cat("y:", head(force.graph2$lay[, 2]), "\n")
  V(graph_with_xy)$x <- force.graph2$lay[, 1]
  V(graph_with_xy)$y <- force.graph2$lay[, 2]
  return(graph_with_xy)
}

