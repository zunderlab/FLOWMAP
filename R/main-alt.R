
#' FLOWMAPfromDF - generate FLOWMAPR analysis results from
#' dataframes in R
#'
#' \code{FLOWMAP} generates FLOWMAPR analysis results from dataframe(s). For more information,
#' as well as a guide for how to choose the best settings for your analysis, go to
#' our GitHub repo \url{https://github.com/zunderlab/FLOWMAP/}.
#' 
#' @param mode FLOWMAPR mode to use in analysis based on starting input,
#' available options include \code{c("single", "multi", "one", "one-special")}
#' @param df single dataframe, list of dataframes with each member belonging to single-cell
#' data from a different timepoint, or a list of lists of dataframes belonging to the same timepoint,
#' but coming from different conditions, to be used in analysis
#' @param project.name Character string to label file output of FLOWMAPR analysis
#' @param clustering.var Vector naming channels to be used to calculate distances/differences
#' between cells for clustering (if requested) and edge-drawing steps
#' @param distance.metric Character specifying which metric to use to calculate between-node distances, valid options include \code{c("manhattan", "euclidean")}
#' @param minimum Numeric value specifying the minimum number of edges that will be allotted
#' during each density-dependent edge-building step of the FLOW-MAP graph, default value is
#' set to \code{2}, no less than 2 is recommended
#' @param maximum Numeric value specifying the maximum number of edges that will be allotted
#' during each density-dependent edge-building step of the FLOW-MAP graph, default value is
#' set to \code{5}, no less than 3 is recommended
#' @param save.folder Directory where all results generated should be saved
#' @param time.col.label Character specifying the name of the channel with the time labels for each cell
#' @param condition.col.label Character specifying the name of the channel with the condition labels for each cell
#' @param name.sort Logical specifying whether to sort FCS file path names alphanumerically or use
#' them in the order supplied by the user
#' @param per Percentile of edge length in fully connected graph that determines hypersphere radius
#' in density calculation
#' @param clustering Logical specifying whether to cluster single-cell data
#' @param seed.X Numeric value for the seed to set for reproducible FLOWMAPR runs
#' @param savePDFs Logical specifying whether to generate PDFs for the resolved graph with
#' nodes colored by each parameter
#' @param which.palette Optional variable, character specifying which color palette to use
#' in generated PDFs, valid options include \code{c("bluered", "jet", "CB")}, where \code{"CB"}
#' is a colorblind-friendly option
#' @param cluster.numbers Optional variable, single numeric or a vector of numerics specifying
#' how many clusters to generate from each separate dataframe
#' @param cluster.mode Character specifying which clustering algorithm to use, valid options include \code{c("hclust", "kmeans")}
#' @return the force-directed layout resolved igraph graph object
#' @export
FLOWMAPfromDF <- function(mode = c("single", "multi", "one"), df, project.name,
                          clustering.var, distance.metric = "manhattan",
                          minimum = 2, maximum = 5, save.folder = getwd(),
                          time.col.label = "Time", condition.col.label = NULL,
                          name.sort = TRUE, per = 10, clustering = FALSE, seed.X = 1,
                          savePDFs = TRUE, which.palette = "bluered", cluster.numbers = NULL,
                          cluster.mode) {
  set.seed(seed.X)
  cat("Seed set to", seed.X, "\n")
  setwd(save.folder)
  df <- RemoveRowNames(df)
  PrintSummaryfromDF(env = parent.frame())

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
    PrintSummaryfromDF(env = parent.frame())
    if (clustering) {
      file.clusters <- ClusterFCS(fcs.files = df, clustering.var = clustering.var,
                                  numcluster = cluster.numbers, distance.metric = distance.metric,
                                  cluster.mode = cluster.mode)
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
    label.key <- GetLabelKeyfromDF(df, time.col.label, condition.col.label)
    PrintSummaryfromDF(env = parent.frame())
    if (clustering) {
      file.clusters <- MultiClusterFCS(.data$fixed.fcs.files, clustering.var = clustering.var, numcluster = cluster.numbers,
                                       distance.metric = distance.metric, cluster.mode = cluster.mode)
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
    PrintSummaryfromDF(env = parent.frame())
    fcs.files <- list()
    fcs.files[[1]] <- df
    if (clustering) {
      file.clusters <- ClusterFCS(fcs.files = fcs.files, clustering.var = clustering.var,
                                  numcluster = cluster.numbers, distance.metric = distance.metric,
                                  cluster.mode = cluster.mode)
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
  } else if (mode == "one-special") {
    check <- CheckDFModeMulti(df)
    cat("check", check, "\n")
    if (check) {
      stop("Unknown 'df' format provided for specified mode!")
    }
    runtype <- "OneTimepoint-MultipleConditions"
    output.folder <- MakeOutFolder(runtype = runtype)
    setwd(output.folder)
    process.results <- GetConditionsfromDF(df.list = df, condition.col.label)
    fixed.df <- process.results$new.df.list
    label.key.special <- process.results$label.key.special
    PrintSummaryfromDF(env = parent.frame())
    if (clustering) {
      file.clusters <- MultiClusterFCS(list.of.files = fixed.df, clustering.var = clustering.var,
                                       numcluster = cluster.numbers, distance.metric = distance.metric,
                                       cluster.mode = cluster.mode)
    } else {
      file.clusters <- ConstructMultiFLOWMAPCluster(fixed.df)
    }
    remodel.FLOWMAP.clusters <- RemodelFLOWMAPClusterList(file.clusters)
    output.graph <- BuildFirstMultiFLOWMAP(list.of.FLOWMAP.clusters = remodel.FLOWMAP.clusters,
                                           per = per, min = minimum, max = maximum,
                                           distance.metric = distance.metric,
                                           clustering.var = clustering.var)
    output.graph <- AnnotateSpecialGraph(output.graph, remodel.FLOWMAP.clusters,
                                         label.key.special)
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
  if (mode == "single") {
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
