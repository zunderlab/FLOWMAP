
#' FLOWMAP - generate FLOWMAPR analysis results from FCS files in R
#'
#' \code{FLOWMAP} generates FLOWMAPR analysis results from FCS files. For more information,
#' as well as a guide for how to choose the best settings for your analysis, go to
#' our GitHub repo \url{https://github.com/zunderlab/FLOWMAP/}.
#'
#' @param mode FLOWMAPR mode to use in analysis based on starting input,
#' available options include \code{c("single", "multi", "one", "static-multi")}
#' @param files File paths for FCS files to be used or a folder containing
#' the FCS files to be used in analysis
#' @param var.remove Vector naming channels to be removed from all downstream analysis, default
#' value is an empty vector, meaning that no channels will be removed
#' @param var.annotate List mapping channel names to user-specified names to properly
#' annotate all FCS file data, default value is \code{NULL}, which will then autogenerate the
#' \code{var.annotate} values using the \code{ConstructVarAnnotate()} function
#' @param clustering.var Vector naming channels to be used to calculate distances/differences
#' between cells for clustering (if requested) and edge-drawing steps
#' @param cluster.numbers A single numeric or a vector of numerics specifying how many clusters
#' to generate from each separate FCS file
#' @param cluster.mode Character specifying which clustering algorithm to use, valid options include \code{c("hclust", "kmeans")}
#' @param distance.metric Character specifying which metric to use to calculate between-node distances, valid options include \code{c("manhattan", "euclidean")}
#' @param density.metric Character string specifying which method to use for local density estimation for edge assignment in graph building
#' @param minimum Numeric value specifying the minimum number of edges that will be allotted
#' during each density-dependent edge-building step of the FLOW-MAP graph, default value is
#' set to \code{2}, no less than 2 is recommended
#' @param maximum Numeric value specifying the maximum number of edges that will be allotted
#' during each density-dependent edge-building step of the FLOW-MAP graph, default value is
#' set to \code{5}, no less than 3 is recommended
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
#' @param umap.settings
#' @param umap.k
#' @param graph.out
#' @return the force-directed layout resolved igraph graph object
#' @export
FLOWMAP <- function(mode = c("single", "multi", "one", "static-multi"),
                    var.remove = c(), var.annotate = NULL, clustering.var, cluster.numbers = 100,
                    cluster.mode = "hclust", distance.metric = "manhattan", umap.settings,
                    files, density.metric = c("kNN", "radius"), minimum = 2, maximum = 5,
                    save.folder = getwd(), subsamples = 200, name.sort = TRUE,
                    downsample = FALSE, seed.X = 1, savePDFs = TRUE, graph.out = c("ForceDirected"),
                    which.palette = "bluered", exclude.pctile = NULL, target.pctile = NULL,
                    target.number = NULL, target.percent = NULL, k = 10, per = 5, ...) {
  set.seed(seed.X)
  cat("Seed set to", seed.X, "\n")
  cat("Mode set to", mode, "\n")
  CheckSettings(mode, save.folder, var.remove, var.annotate,
                clustering.var, cluster.numbers, cluster.mode,
                distance.metric, minimum, maximum,
                subsamples, which.palette)

  if (downsample) {
    CheckDownsampleSettings(exclude.pctile, target.pctile,
                            target.number, target.percent)
  }
  setwd(save.folder)
#SINGLE====
  if (mode == "single") {
    check <- CheckModeSingle(files)
    cat("check", check, "\n")
    if (check[1]) {
      stop("Unknown 'files' format provided for specified mode!")
    }
    runtype <- "SingleFLOWMAP"
    if (density.metric == "radius") {
      output.folder <- MakeOutFolder(runtype = runtype, per = per, maximum = maximum, minimum = minimum)
    } else if (density.metric == "kNN") {
      output.folder <- MakeOutFolder(runtype = runtype, k = k, maximum = maximum, minimum = minimum)
    }

    setwd(output.folder)
    PrintSummary(env = parent.frame())
    if (check[2] == "FCS") {
      fcs.file.names <- files
    }
    if (check[2] == "folder") {
      fcs.file.names <- GetFCSNames(folder = files, sort = name.sort)
    }
    orig.times <- ParseTimes(fcs.file.names, name.sort = name.sort)
    file.name <- fcs.file.names[1]
    if (is.null(var.annotate)) {
      var.annotate <- ConstructVarAnnotate(file.name)
      assign("var.annotate", var.annotate, envir = .GlobalEnv)
    }
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
                                numcluster = cluster.numbers, distance.metric = distance.metric,
                                cluster.mode = cluster.mode)
    if (cluster.mode != "none") {
      cat("Upsampling all clusters to reflect Counts of entire file", "\n")
      file.clusters <- Upsample(fcs.file.names, file.clusters, fcs.files, var.remove, var.annotate, clustering.var)
    }

    #Build FLOWMAP single====
    if (density.metric == "radius") {
      results <- BuildFLOWMAP(FLOWMAP.clusters = file.clusters, per = per, min = minimum,
                                 max = maximum, distance.metric = distance.metric,
                                 clustering.var = clustering.var)
    } else if (density.metric == "kNN") {
      results <- BuildFLOWMAPkNN(FLOWMAP.clusters = file.clusters, k = k, min = minimum,
                              max = maximum, distance.metric = distance.metric,
                              clustering.var = clustering.var)
    }
    graph <- results$output.graph
    knn.indexes <- results$knn.indexes
    knn.distances <- results$knn.distances

#MULTI====
  } else if (mode == "multi") {
    check <- CheckModeMulti(files)
    cat("check", check, "\n")
    if (check[1]) {
      stop("Unknown 'files' format provided for specified mode!")
    }
    runtype <- "MultiFLOWMAP"
    if (density.metric == "radius") {
      output.folder <- MakeOutFolder(runtype = runtype, per = per, maximum = maximum, minimum = minimum)
    } else if (density.metric == "kNN") {
      output.folder <- MakeOutFolder(runtype = runtype, k = k, maximum = maximum, minimum = minimum)
    }
    setwd(output.folder)
    PrintSummary(env = parent.frame())
    if (check[2] == "list") {
      fcs.file.names <- files
      orig.times <- MultiListParseTimes(fcs.file.names, name.sort)
    }
    if (check[2] == "subfolder") {
      fcs.file.names <- GetMultiFCSNames(folder = files, sort = name.sort)
      orig.times <- MultiFolderParseTimes(files, fcs.file.names, name.sort = name.sort)
    }
    file.name <- fcs.file.names[[1]][1]
    if (is.null(var.annotate)) {
      var.annotate <- ConstructVarAnnotate(file.name)
      assign("var.annotate", var.annotate, envir = .GlobalEnv)
    }
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
                                     distance.metric = distance.metric, cluster.mode = cluster.mode)
    cat("Upsampling all clusters to reflect Counts of entire file", "\n")
    file.clusters <- MultiUpsample(fcs.file.names, file.clusters, fcs.files, var.remove, var.annotate, clustering.var)
    #Build FLOWMAP multi====
    if (density.metric == "radius") {
      graph <- BuildMultiFLOWMAP(file.clusters, per = per, min = minimum,
                                 max = maximum, distance.metric = distance.metric,
                                 label.key = label.key, clustering.var = clustering.var)
    }
    ##kNN does not yet work with multi
    # else if (density.metric == "kNN") {
    #   graph <- BuildMultiFLOWMAPkNN(file.clusters, k = maximum, min = minimum,
    #                              max = maximum, distance.metric = distance.metric,
    #                              label.key = label.key, clustering.var = clustering.var)
    # }
#ONE====
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
    if (density.metric == "radius") {
      output.folder <- MakeOutFolder(runtype = runtype, per = per, maximum = maximum, minimum = minimum)
    } else if (density.metric == "kNN") {
      output.folder <- MakeOutFolder(runtype = runtype, k = k, maximum = maximum, minimum = minimum)
    }
    setwd(output.folder)
    PrintSummary(env = parent.frame())
    if (is.null(var.annotate)) {
      var.annotate <- ConstructVarAnnotate(file.name)
      assign("var.annotate", var.annotate, envir = .GlobalEnv)
    }
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
                                numcluster = cluster.numbers, distance.metric = distance.metric,
                                cluster.mode = cluster.mode)
    cat("Upsampling all clusters to reflect Counts of entire file", "\n")
    file.clusters <- Upsample(file.name, file.clusters, fcs.file, var.remove, var.annotate, clustering.var)
    ##Build FLOWMAP one====
    if (density.metric == "radius") {
      first.results <- BuildFirstFLOWMAP(FLOWMAP.clusters = file.clusters,
                                         per = per, min = minimum, max = maximum,
                                         distance.metric = distance.metric,
                                         clustering.var = clustering.var)
    } else if (density.metric == "kNN") {
      first.results <- BuildFirstFLOWMAPkNN(FLOWMAP.clusters = file.clusters,
                                         k = k, min = minimum, max = maximum,
                                         distance.metric = distance.metric,
                                         clustering.var = clustering.var)
    }
    output.graph <- first.results$output.graph
    output.graph <- AnnotateGraph(output.graph = output.graph,
                                  FLOWMAP.clusters = file.clusters)
    graph <- output.graph
    knn.indexes <- results$knn.indexes
    knn.distances <- results$knn.distances

#STATIC-MULTI====
  } else if (mode == "static-multi") {
    check <- CheckModeSingle(files) # static-multi mode will look the same as single
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
    file.name <- fcs.file.names[1]
    runtype <- "OneTimepoint-MultipleConditions"
    if (density.metric == "radius") {
      output.folder <- MakeOutFolder(runtype = runtype, per = per, maximum = maximum, minimum = minimum)
    } else if (density.metric == "kNN") {
      output.folder <- MakeOutFolder(runtype = runtype, k = k, maximum = maximum, minimum = minimum)
    }
    setwd(output.folder)
    PrintSummary(env = parent.frame())
    if (is.null(var.annotate)) {
      var.annotate <- ConstructVarAnnotate(file.name)
      assign("var.annotate", var.annotate, envir = .GlobalEnv)
    }
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
    fcs.files.old <- fcs.files
    fcs.files.list <- list()
    process.results <- ProcessConditions(fcs.files.old, fcs.file.names)
    label.key.special <- process.results$label.key.special
    fcs.files.list[[1]] <- process.results$fixed.files
    file.clusters <- MultiClusterFCS(list.of.files = fcs.files.list, clustering.var = clustering.var,
                                     numcluster = cluster.numbers, distance.metric = distance.metric,
                                     cluster.mode = cluster.mode)
    cat("Upsampling all clusters to reflect Counts of entire file", "\n")
    temp.files.list <- list()
    temp.files.list[[1]] <- fcs.files
    temp.name.list <- list()
    temp.name.list[[1]] <- fcs.file.names
    file.clusters <- MultiUpsample(temp.name.list, file.clusters, temp.files.list, var.remove, var.annotate, clustering.var)
    remodel.FLOWMAP.clusters <- RemodelFLOWMAPClusterList(file.clusters)
    ##Build FLOWMAP one-special ====
    if (density.metric == "radius") {
      output.graph <- BuildFirstMultiFLOWMAP(list.of.FLOWMAP.clusters = remodel.FLOWMAP.clusters,
                                             k = maximum, min = minimum, max = maximum,
                                             distance.metric = distance.metric,
                                             clustering.var = clustering.var)
    }
    output.graph <- AnnotateSpecialGraph(output.graph, remodel.FLOWMAP.clusters,
                                         label.key.special)
    graph <- output.graph

  } else {
    stop("Unknown mode!")
  }
  #Make layout output ====
  if ("ForceDirected" %in% graph.out) {
    graph.xy <- RunForceDirectedLayout(mode=mode, file.name=file.name, graph=graph,
                                       orig.times=orig.times, which.palette=which.palette)
  }
  if ("UMAP"  %in% graph.out) {
    if (umap.k > minimum) {
      stop("Too few edges in input graph to generate UMAP layout with specified umap.k value!")
    } else {
      RunUMAPlayout(graph=graph, file.clusters=file.clusters, file.name=file.name,umap.settings=umap.settings)
    }
  }

}#end FLOWMAP function
