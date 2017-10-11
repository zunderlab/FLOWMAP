
#' @export
FLOWMAP <- function(mode = c("single", "multi", "one"), files, var.remove,
                    var.annotate, clustering.var, cluster.numbers = 100,
                    distance.metric = "manhattan", minimum = 2, maximum = 5,
                    per = 1, save.folder = getwd(), subsamples = 200, name.sort = TRUE,
                    downsample = FALSE, seed.X = 1, savePDFs = TRUE,
                    which.palette = "bluered", exclude.pctile = NULL, target.pctile = NULL,
                    target.number = NULL, target.percent = NULL, ...) {
  # optional variables
  # starting.files = c("FCS", "cluster_matrix")
  cat("exclude_pctile is", exclude.pctile, "\n")
  cat("target_pctile is", target.pctile, "\n")
  cat("target_number is", target.number, "\n")
  cat("target_percent is", target.percent, "\n")
  
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
      # fcs.files <- DownsampleFCS(fcs.file.names, clustering.var, channel.annotate = var.annotate,
      #                            channel.remove = var.remove, exclude.pctile = exclude.pctile,
      #                            target.pctile = target.pctile, target.number = target.number,
      #                            target.percent = target.percent)
      
      cat("exclude_pctile is", exclude.pctile, "\n")
      cat("target_pctile is", target.pctile, "\n")
      cat("target_number is", target.number, "\n")
      cat("target_percent is", target.percent, "\n")
      fcs.file.names <- DownsampleFCS(fcs.file.names, clustering.var,
                                      distance.metric,  exclude.pctile = exclude.pctile,
                                      target.pctile = target.pctile, target.number = target.number,
                                      target.percent = target.percent)
      subsamples <- FALSE
      fcs.files <- LoadCleanFCS(fcs.file.names = fcs.file.names, channel.remove = var.remove,
                                channel.annotate = var.annotate, subsamples = subsamples)
    } else {
      fcs.files <- LoadCleanFCS(fcs.file.names = fcs.file.names, channel.remove = var.remove,
                                channel.annotate = var.annotate, subsamples = subsamples)
    }
    file.clusters <- ClusterFCS(fcs.files = fcs.files, clustering.var = clustering.var,
                                numcluster = cluster.numbers, distance.metric = distance.metric)
    if (downsample) {
      cat("Upsampling all clusters to reflect Counts prior to SPADE downsampling", "\n")
      file.clusters <- Upsample(file.clusters)
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
      # fcs.files <- MultiDownsampleFCS(fcs.file.names, clustering.var, channel.annotate = var.annotate,
      #                                 channel.remove = var.remove, exclude.pctile = exclude.pctile,
      #                                 target.pctile = target.pctile, target.number = target.number,
      #                                 target.percent = target.percent)
      cat("exclude_pctile is", exclude.pctile, "\n")
      cat("target_pctile is", target.pctile, "\n")
      cat("target_number is", target.number, "\n")
      cat("target_percent is", target.percent, "\n")
      fcs.file.names <- MultiDownsampleFCS(fcs.file.names, clustering.var,
                                           distance.metric,  exclude.pctile = exclude.pctile,
                                           target.pctile = target.pctile, target.number = target.number,
                                           target.percent = target.percent)
      subsamples <- FALSE
      fcs.files <- LoadMultiCleanFCS(fcs.file.names, var.remove, var.annotate,
                                     subsamples = subsamples)
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
      file.clusters <- MultiUpsample(file.clusters)
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
    fcs.file <- LoadCleanFCS(fcs.file.names = fcs.file.names, channel.remove = var.remove,
                             channel.annotate = var.annotate, subsamples = subsamples)
    file.clusters <- ClusterFCS(fcs.files = fcs.file, clustering.var = clustering.var,
                                numcluster = cluster.numbers, distance.metric = distance.metric)
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
