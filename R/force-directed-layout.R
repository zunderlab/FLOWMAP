#' @importFrom vite layout_forceatlas2
ForceDirectedXY <- function(graph) {
  V(graph)$size <- rep(20, vcount(graph))
  force.graph1 <- vite::layout_forceatlas2(graph, iter = 10000,
                                                stopping.tolerance = 0.001,
                                                prevent.overlap = FALSE)
  graph.with.xy <- graph
  V(graph.with.xy)$x <- force.graph1$lay[, 1]
  V(graph.with.xy)$y <- force.graph1$lay[, 2]
  force.graph2 <- vite::layout_forceatlas2(graph.with.xy, iter = 1000,
                                                stopping.tolerance = 0.001,
                                                prevent.overlap = TRUE)
  V(graph.with.xy)$x <- force.graph2$lay[, 1]
  V(graph.with.xy)$y <- force.graph2$lay[, 2]
  return(graph.with.xy)
}
