
#' FLOWMAP - generate FLOWMAPR analysis results from FCS files in R
#'
#' \code{FLOWMAP} generates FLOWMAPR analysis results from FCS files. For more information,
#' as well as a guide for how to choose the best settings for your analysis, go to
#' our GitHub repo \url{https://github.com/zunderlab/FLOWMAP/}.
#' 
#' @param mode FLOWMAPR mode to use in analysis based on starting input,
#' available options include \code{c("single", "multi", "one")}
#' @param files File paths for FCS files to be used or a folder containing
#' the FCS files to be used in analysis
#' @param var.remove Vector naming channels to be removed from all downstream analysis
#' @param var.annotate List mapping channel names to user-specified names to properly
#' annotate all FCS file data
#' @param clustering.var Vector naming channels to be used to calculate distances/differences
#' between cells for clustering (if requested) and edge-drawing steps
#' @param cluster.numbers A single numeric or a vector of numerics specifying how many clusters
#' to generate from each separate FCS file
#' @param distance.metric Character \code{c("manhattan", "euclidean")}
#' @param minimum Numeric value specifying the minimum number of edges that will be allotted
#' during each density-dependent edge-building step of the FLOW-MAP graph, default value is 
#' set to \code{2}, no less than 2 is recommended
#' @param maximum Numeric value specifying the maximum number of edges that will be allotted
#' during each density-dependent edge-building step of the FLOW-MAP graph, default value is
#' set to \code{5}, no less than 3 is recommended
#' @param per Numeric value specifying the top n% of edges by strength that will be used to assess
#' density and allot edges during the density-dependent edge-building step of the FLOW-MAP
#' graph, default value is set to \code{1}, changing this value is not advised
#' @param save.folder Directory where all results generated should be saved
#' @param subsamples A single numeric or a vector of numerics specifying how many cells to
#' subsample from each FCS file
#' @param name.sort Logical specifying whether to sort FCS file path names alphanumerically or use
#' them in the order supplied by the user
#' @param downsample Logical specifying whether to use SPADE density-dependent downsampling
#' @param seed.X Numeric value for the seed to set for reproducible FLOWMAPR runs
#' @param savePDFs Logical specifying whether to generate PDFs for the resolved graph with
#' nodes colored by each parameter
#' @param which.palette Optional variable, character specifying which color palette to use
#' in generated PDFs, valid options include \code{c("bluered", "jet", "CB")}, where \code{"CB"}
#' is a colorblind-friendly option
#' @param exclude.pctile Optional variable, numeric value for the downsampling_exclude_pctile variable
#' used as described in the SPADE driver function, see the documentation for the spade package at
#' \url{https://github.com/nolanlab/spade}
#' @param target.pctile Optional variable, numeric value for the downsampling_target_pctile variable
#' used as described in the SPADE driver function, see the documentation for the spade package at
#' \url{https://github.com/nolanlab/spade}
#' @param target.number Optional variable, numeric value for the downsampling_target_number variable
#' used as described in the SPADE driver function, see the documentation for the spade package at
#' \url{https://github.com/nolanlab/spade}
#' @param target.percent Optional variable, numeric value for the downsampling_target_percent variable
#' used as described in the SPADE driver function, see the documentation for the spade package at
#' \url{https://github.com/nolanlab/spade}
#' @return the force-directed layout resolved igraph graph object
#' @export
FLOWMAP <- function(mode = c("single", "multi", "one"), files, var.remove,
                    var.annotate, clustering.var, cluster.numbers = 100,
                    distance.metric = "manhattan", minimum = 2, maximum = 5,
                    per = 1, save.folder = getwd(), subsamples = 200, name.sort = TRUE,
                    downsample = FALSE, seed.X = 1, savePDFs = TRUE,
                    which.palette = "bluered", exclude.pctile = NULL, target.pctile = NULL,
                    target.number = NULL, target.percent = NULL, ...) {
  set.seed(seed.X)
  cat("Seed set to", seed.X, "\n")
  setwd(save.folder)
  if (mode == "single") {
    check <- CheckModeSingle(files)
    cat("check", check, "\n")
    if (check[1]) {
      stop("Unknown 'files' format provided for specified mode!")
    }
    runtype <- "SingleFLOWMAP"
    output.folder <- MakeOutFolder(runtype = runtype)
    setwd(output.folder)
    if (check[2] == "FCS") {
      fcs.file.names <- files
    }
    if (check[2] == "folder") {
      fcs.file.names <- GetFCSNames(folder = files, sort = name.sort)
    }
    orig.times <- ParseTimes(fcs.file.names, name.sort = name.sort)
    file.name <- fcs.file.names[1]
    if (downsample) {
      cat("Downsampling all files using SPADE downsampling", "\n")
      fcs.files <- DownsampleFCS(fcs.file.names, clustering.var, channel.annotate = var.annotate,
                                 channel.remove = var.remove, exclude.pctile = exclude.pctile,
                                 target.pctile = target.pctile, target.number = target.number,
                                 target.percent = target.percent, transform = TRUE)
    } else {
      fcs.files <- LoadCleanFCS(fcs.file.names = fcs.file.names, channel.remove = var.remove,
                                channel.annotate = var.annotate, subsamples = subsamples)
    }
    file.clusters <- ClusterFCS(fcs.files = fcs.files, clustering.var = clustering.var,
                                numcluster = cluster.numbers, distance.metric = distance.metric)
    if (downsample) {
      cat("Upsampling all clusters to reflect Counts prior to SPADE downsampling", "\n")
      file.clusters <- Upsample(fcs.file.names, file.clusters, fcs.files, var.remove, var.annotate, clustering.var)
    }
    results <- BuildFLOWMAP(FLOWMAP.clusters = file.clusters, per = per, min = minimum,
                            max = maximum, distance.metric = distance.metric,
                            clustering.var = clustering.var)
    graph <- results$output.graph
  } else if (mode == "multi") {
    check <- CheckModeMulti(files)
    cat("check", check, "\n")
    if (check[1]) {
      stop("Unknown 'files' format provided for specified mode!")
    }
    runtype <- "MultiFLOWMAP"
    output.folder <- MakeOutFolder(runtype = runtype)
    setwd(output.folder)
    if (check[2] == "list") {
      fcs.file.names <- files
      orig.times <- MultiListParseTimes(fcs.file.names, name.sort)
    }
    if (check[2] == "subfolder") {
      fcs.file.names <- GetMultiFCSNames(folder = files, sort = name.sort)
      orig.times <- MultiFolderParseTimes(files, fcs.file.names, name.sort = name.sort)
    }
    file.name <- fcs.file.names[[1]][1]
    if (downsample) {
      cat("Downsampling all files using SPADE downsampling", "\n")
      fcs.files <- MultiDownsampleFCS(fcs.file.names, clustering.var, channel.annotate = var.annotate,
                                      channel.remove = var.remove, exclude.pctile = exclude.pctile,
                                      target.pctile = target.pctile, target.number = target.number,
                                      target.percent = target.percent, transform = TRUE)
    } else {
      fcs.files <- LoadMultiCleanFCS(fcs.file.names, var.remove, var.annotate,
                                     subsamples = subsamples)
    }
    fcs.files.conversion <- ConvertNumericLabel(fcs.files)
    fixed.fcs.files <- fcs.files.conversion$fixed.list.FCS.files
    label.key <- fcs.files.conversion$label.key
    file.clusters <- MultiClusterFCS(fixed.fcs.files, clustering.var = clustering.var, numcluster = cluster.numbers,
                                     distance.metric = distance.metric)
    if (downsample) {
      cat("Upsampling all clusters to reflect Counts prior to SPADE downsampling", "\n")
      file.clusters <- MultiUpsample(fcs.file.names, file.clusters, fcs.files, var.remove, var.annotate, clustering.var)
    }
    graph <- BuildMultiFLOWMAP(file.clusters, per = per, min = minimum,
                               max = maximum, distance.metric = distance.metric,
                               label.key = label.key, clustering.var = clustering.var)
  } else if (mode == "one") {
    check <- CheckModeOne(files)
    cat("check", check, "\n")
    if (check[1]) {
      stop("Unknown 'files' format provided for specified mode!")
    }
    if (check[2] == "FCS") {
      fcs.file.names <- files
    }
    if (check[2] == "folder") {
      fcs.file.names <- GetFCSNames(folder = files, sort = name.sort)
    }
    file.name <- fcs.file.names
    runtype <- "OneTimepoint"
    output.folder <- MakeOutFolder(runtype = runtype)
    setwd(output.folder)
    
    if (downsample) {
      cat("Downsampling all files using SPADE downsampling", "\n")
      fcs.file <- DownsampleFCS(file.name, clustering.var, channel.annotate = var.annotate,
                                channel.remove = var.remove, exclude.pctile = exclude.pctile,
                                target.pctile = target.pctile, target.number = target.number,
                                target.percent = target.percent, transform = TRUE)
    } else {
      fcs.file <- LoadCleanFCS(fcs.file.names = file.name, channel.remove = var.remove,
                               channel.annotate = var.annotate, subsamples = subsamples)
    }
    file.clusters <- ClusterFCS(fcs.files = fcs.file, clustering.var = clustering.var,
                                numcluster = cluster.numbers, distance.metric = distance.metric)
    if (downsample) {
      cat("Upsampling all clusters to reflect Counts prior to SPADE downsampling", "\n")
      file.clusters <- Upsample(file.name, file.clusters, fcs.files, var.remove, var.annotate, clustering.var)
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
  file.name <- paste(unlist(strsplit(basename(file.name), "\\."))[1], "FLOW-MAP", sep = "_")
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
  PrintSummary(env = parent.frame())
  MakeFLOWMAPRFile(env = parent.frame())
  if (savePDFs) {
    cat("Printing pdfs.", "\n")
    ConvertToPDF(graphml.file = fixed.file, which.palette = which.palette)
  }
  return(graph.xy)
}
