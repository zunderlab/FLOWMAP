
#' @export
FLOWMAPfromDF <- function(mode = c("single", "multi", "one"), df, project.name,
                          clustering.var, distance.metric = "manhattan",
                          minimum = 2, maximum = 5, per = 1, save.folder = getwd(),
                          time.col.label = "Time", condition.col.label = NULL,
                          name.sort = TRUE, downsample = FALSE, clustering = FALSE,
                          seed.X = 1, savePDFs = TRUE, which.palette = "bluered",
                          cluster.numbers = NULL, exclude.pctile = NULL, target.pctile = NULL,
                          target.number = NULL, target.percent = NULL, ...) {
  set.seed(seed.X)
  cat("Seed set to", seed.X, "\n")
  setwd(save.folder)
  if (mode == "single") {
    check <- CheckDFModeSingle(df)
    cat("check", check, "\n")
    if (check) {
      stop("Unknown 'df' format provided for specified mode!")
    }
    runtype <- "SingleFLOWMAP"
    output.folder <- MakeOutFolder(runtype = runtype)
    setwd(output.folder)
    orig.times <- GetOrigTimesfromDF(df, time.col.label, name.sort = name.sort)
    df <- StripTimesfromDFList(df, time.col.label)
    print("removed times")
    df <- ShuffleCells(df)
    if (clustering) {
      file.clusters <- ClusterFCS(fcs.files = df, clustering.var = clustering.var,
                                  numcluster = cluster.numbers, distance.metric = distance.metric)
    } else {
      file.clusters <- ConstructSingleFLOWMAPCluster(df)
    }
    results <- BuildFLOWMAP(FLOWMAP.clusters = file.clusters, per = per, min = minimum,
                            max = maximum, distance.metric = distance.metric,
                            clustering.var = clustering.var)
    graph <- results$output.graph
  } else if (mode == "multi") {
    check <- CheckDFModeMulti(df)
    cat("check", check, "\n")
    if (check) {
      stop("Unknown 'df' format provided for specified mode!")
    }
    runtype <- "MultiFLOWMAP"
    output.folder <- MakeOutFolder(runtype = runtype)
    setwd(output.folder)
    df <- MultiShuffleCells(df)
    label.key <- GetLabelKeyfromDF(df, time.col.label, condition.col.label)
    if (clustering) {
      file.clusters <- MultiClusterFCS(fixed.fcs.files, clustering.var = clustering.var, numcluster = cluster.numbers,
                                       distance.metric = distance.metric)
    } else {
      file.clusters <- ConstructMultiFLOWMAPCluster(df)
    }
    graph <- BuildMultiFLOWMAP(file.clusters, per = per, min = minimum,
                               max = maximum, distance.metric = distance.metric,
                               label.key = label.key, clustering.var = clustering.var)
  } else if (mode == "one") {
    check <- CheckDFModeOne(df)
    cat("check", check, "\n")
    if (check) {
      stop("Unknown 'df' format provided for specified mode!")
    }
    runtype <- "OneTimepoint"
    output.folder <- MakeOutFolder(runtype = runtype)
    setwd(output.folder)
    if (clustering) {
      file.clusters <- ClusterFCS(fcs.files = df, clustering.var = clustering.var,
                                  numcluster = cluster.numbers, distance.metric = distance.metric)
    } else {
      file.clusters <- ConstructOneFLOWMAPCluster(df)
    }
    first.results <- BuildFirstFLOWMAP(FLOWMAP.clusters = file.clusters,
                                       per = per, min = minimum, max = maximum,
                                       distance.metric = distance.metric,
                                       clustering.var = clustering.var)
    output.graph <- first.results$output.graph
    output.graph <- AnnotateGraph(output.graph = output.graph,
                                  FLOWMAP.clusters = file.clusters)
    graph <- output.graph
  } else {
    stop("Unknown mode!")
  }
  file.name <- paste(project.name, "FLOW-MAP", sep = "_")
  ConvertToGraphML(output.graph = graph, file.name = file.name)
  graph.xy <- ForceDirectedXY(graph = graph)
  file.name.xy <- paste(file.name, "xy", sep = "_")
  final.file.name <- ConvertToGraphML(output.graph = graph.xy, file.name = file.name.xy)
  fixed.file.name <- paste(file.name.xy, "orig_time", sep = "_")
  if (mode != "one") {
    fixed.graph <- ConvertOrigTime(graph.xy, orig.times)
  } else {
    fixed.graph <- graph.xy
  }
  fixed.file <- ConvertToGraphML(output.graph = fixed.graph, file.name = fixed.file.name)
  # PrintSummary(env = parent.frame())
  # MakeFLOWMAPRFile(env = parent.frame())
  if (savePDFs) {
    cat("Printing pdfs.", "\n")
    ConvertToPDF(graphml.file = fixed.file, which.palette = which.palette)
  }
  return(graph.xy)
}
