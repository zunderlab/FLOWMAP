
MakePrintVarAnnotate <- function(var.annotate) {
  if (length(var.annotate) == 0 && is.null(var.annotate)) {
    p.var.annotate <- paste("var.annotate", " <- ", "NULL", sep = "")
  } else if (length(var.annotate) == 0 && !is.null(var.annotate)) {
    p.var.annotate <- paste("var.annotate", " <- ", "list()", sep = "")
  } else {
    for (i in 1:length(var.annotate)) {
      marker.i <- paste("'", names(var.annotate)[i], "'", " = ",
                        "'", unname(unlist(var.annotate))[i], "'", sep = "")
      if (i == 1) {
        p.var.annotate <- marker.i
      } else {
        p.var.annotate <- paste(p.var.annotate, marker.i, sep = ", ")
      }
    }
    p.var.annotate <- paste("list(", p.var.annotate, ")", sep = "")
    p.var.annotate <- paste("var.annotate", " <- ", p.var.annotate, sep = "")
  }
  return(p.var.annotate)
}

MakePrintVarRemove <- function(var.remove) {
  if (is.null(var.remove)) {
    p.var.remove <- "c()"
  } else {
    p.var.remove <- paste("'", var.remove, "'", sep = "")
    p.var.remove <- paste(p.var.remove, collapse = ", ")
    p.var.remove <- paste("c(", p.var.remove, ")", sep = "")
  }
  p.var.remove <- paste("var.remove", " <- ", p.var.remove, sep = "")
  return(p.var.remove)
}

MakePrintClusteringVar <- function(clustering.var) {
  p.clustering.var <- paste("'", clustering.var, "'", sep = "")
  p.clustering.var <- paste(p.clustering.var, collapse = ", ")
  p.clustering.var <- paste("c(", p.clustering.var, ")", sep = "")
  p.clustering.var <- paste("clustering.var", " <- ", p.clustering.var, sep = "")
  return(p.clustering.var)
}

MakePrintFiles <- function(files) {
  if (length(files) == 1) {
    if (exists("globe.raw.FCS.dir")) {
      files <- paste(.data$globe.raw.FCS.dir, files, sep = "/")
    }
    p.files <- MakePrintChar("files", files)
  } else if (length(files) > 1 & !is.list(files)) {
    p.files <- paste("'", files, "'", sep = "")
    p.files <- paste(p.files, collapse = ", ")
    p.files <- paste("c(", p.files, ")", sep = "")
    p.files <- paste("files", " <- ", p.files, sep = "")
  } else if (length(files) > 1 & is.list(files)){
    for (i in 1:length(files)) {
      files.in.i <- paste("'", files[[i]], "'", sep = "")
      files.in.i <- paste(files.in.i, collapse = ", ")
      files.in.i <- paste("c(", files.in.i, ")", sep = "")
      if (i == 1) {
        p.files <- files.in.i
      } else {
        p.files <- paste(p.files, files.in.i, sep = ", ")
      }
    }
    p.files = gsub("\\\\", "/", p.files)
    print(p.files)
    p.files <- paste("list(", p.files, ")", sep = "")
    p.files <- paste("files", " <- ", p.files, sep = "")
  } else {
    warning("Files format not recognizable for printing!")
    p.files <- MakePrintChar("files", files)
  }
  return(p.files)
}

MakePrintNum <- function(var.assgn, num) {
  p.num <- paste(var.assgn, " <- ", as.character(num), sep = "")
  return(p.num)
}

MakePrintChar <- function(var.assgn, char) {
  p.char <- paste(var.assgn, " <- ", "'", char, "'", sep = "")
  return(p.char)
}

MakeFLOWMAPRFile <- function(env = parent.frame()) {
  cat("Generating .R file of FLOWMAPR run", "\n")
  p.files <- MakePrintFiles(env$files)
  p.mode <- MakePrintChar("mode", env$mode)
  p.save.folder <- MakePrintChar("save.folder", env$save.folder)
  p.minimum <- MakePrintNum("minimum", env$minimum)
  p.maximum <- MakePrintNum("maximum", env$maximum)
  p.distance.metric <- MakePrintChar("distance.metric", env$distance.metric)
  p.cluster.numbers <- MakePrintNum("cluster.numbers", env$cluster.numbers)
  p.var.annotate <- MakePrintVarAnnotate(env$var.annotate)
  p.var.remove <- MakePrintVarRemove(env$var.remove)
  p.clustering.var <- MakePrintClusteringVar(env$clustering.var)
  p.seed.X <- MakePrintNum("seed.X", env$seed.X)
  p.set.seed.X <- "set.seed(seed.X)"
  p.name.sort = MakePrintNum("name.sort", env$name.sort)
  p.downsample = MakePrintNum("downsample", env$downsample)
  p.savePDFs = MakePrintNum("savePDFs", env$savePDFs)
  p.which.palette = MakePrintChar("which.palette", env$which.palette)
  
  if (env$downsample) {
    p.exclude.pctile <- MakePrintNum("exclude.pctile", env$exclude.pctile)
    p.target.pctile <- MakePrintNum("target.pctile", env$target.pctile)
    if (is.null(env$target.number)) {
      p.target.number <- "target.number <- NULL"
    } else {
      p.target.number <- MakePrintNum("target.number", env$target.number)
    }
    if (is.null(env$target.percent)) {
      p.target.percent <- "target.percent <- NULL"
    } else {
      p.target.percent <- MakePrintNum("target.percent", env$target.percent)
    }
    p.subsamples <- "subsamples <- FALSE"
  } else {
    p.subsamples <- MakePrintNum("subsamples", env$subsamples)
  }
  if (env$downsample) {
    p.FLOWMAP <- "FLOWMAP(seed.X = seed.X, files = files, var.remove = var.remove,
    var.annotate = var.annotate, clustering.var = clustering.var,
    cluster.numbers = cluster.numbers, subsamples = subsamples,
    distance.metric = distance.metric, minimum = minimum,
    maximum = maximum, save.folder = save.folder,
    mode = mode, name.sort = name.sort, downsample = downsample,
    savePDFs = savePDFs, which.palette = which.palette,
    exclude.pctile = exclude.pctile, target.pctile = target.pctile,
    target.number = target.number, target.percent = target.percent)"
  } else {
    p.FLOWMAP <- "FLOWMAP(seed.X = seed.X, files = files, var.remove = var.remove,
    var.annotate = var.annotate, clustering.var = clustering.var,
    cluster.numbers = cluster.numbers, subsamples = subsamples,
    distance.metric = distance.metric, minimum = minimum,
    maximum = maximum, save.folder = save.folder,
    mode = mode, name.sort = name.sort, downsample = downsample,
    savePDFs = savePDFs, which.palette = which.palette)"
  }
  output.file <- "run_FLOWMAPR.R"
  output <- file(output.file)
  p.setup1 <- "rm(list = ls())"
  p.setup2 <- "library(FLOWMAPR)"
  if (env$downsample) {
    writeLines(c(p.setup1, p.setup2, p.files, p.mode, p.save.folder, p.minimum,
                 p.maximum, p.distance.metric, p.cluster.numbers,
                 p.var.annotate, p.var.remove, p.clustering.var,
                 p.seed.X, p.set.seed.X, p.subsamples, p.name.sort,
                 p.downsample, p.savePDFs, p.which.palette,
                 p.exclude.pctile, p.target.pctile, p.target.number,
                 p.target.percent, p.FLOWMAP), output)
  } else {
    writeLines(c(p.setup1, p.setup2, p.files, p.mode, p.save.folder, p.minimum,
                 p.maximum, p.distance.metric, p.cluster.numbers,
                 p.var.annotate, p.var.remove, p.clustering.var,
                 p.seed.X, p.set.seed.X, p.subsamples, p.name.sort,
                 p.downsample, p.savePDFs, p.which.palette, p.FLOWMAP), output)
  }
  close(output)
}
