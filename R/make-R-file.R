
MakePrintVarAnnotate <- function(var.annotate) {
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

MakePrintNum <- function(var.assgn, num) {
  p.num <- paste(var.assgn, " <- ", as.character(num), sep = "")
  return(p.num)
}

MakePrintChar <- function(var.assgn, char) {
  p.char <- paste(var.assgn, " <- ", "'", char, "'", sep = "")
  return(p.char)
}

MakeFLOWMAPRFile <- function(env = parent.frame()) {
  p.files <- MakePrintChar("files", files)
  p.mode <- MakePrintChar("mode", mode)
  p.save.folder <- MakePrintChar("save.folder", save.folder)
  p.per <- MakePrintNum("per", per)
  p.minimum <- MakePrintNum("minimum", minimum)
  p.maximum <- MakePrintNum("maximum", maximum)
  p.distance.metric <- MakePrintChar("distance.metric", distance.metric)
  p.cluster.numbers <- MakePrintNum("cluster.numbers", cluster.numbers)
  p.var.annotate <- MakePrintVarAnnotate(var.annotate)
  p.var.remove <- MakePrintVarRemove(var.remove)
  p.clustering.var <- MakePrintClusteringVar(clustering.var)
  p.seed.X <- MakePrintNum("seed.X", seed.X)
  p.set.seed.X <- "set.seed(seed.X)"
  p.name.sort = MakePrintNum("name.sort", name.sort)
  p.downsample = MakePrintNum("downsample", downsample)
  p.savePDFs = MakePrintNum("savePDFs", savePDFs)
  p.which.palette = MakePrintChar("which.palette", which.palette)
  
  if (downsample) {
    p.exclude.pctile <- MakePrintNum("exclude.pctile", exclude.pctile)
    p.target.pctile <- MakePrintNum("target.pctile", target.pctile)
    if (is.null(target.number)) {
      p.target.number <- "target.number <- NULL"
    } else {
      p.target.number <- MakePrintNum("target.number", target.number)
    }
    if (is.null(target.percent)) {
      p.target.percent <- "target.percent <- NULL"
    } else {
      p.target.percent <- MakePrintNum("target.percent", target.percent)
    }
    p.subsamples <- "subsamples <- FALSE"
  } else {
    p.subsamples <- MakePrintNum("subsamples", subsamples)
  }
  if (downsample) {
    p.FLOWMAP <- "FLOWMAP(seed.X = seed.X, files = files, var.remove = var.remove,
    var.annotate = var.annotate, clustering.var = clustering.var,
    cluster.numbers = cluster.numbers, subsamples = subsamples,
    distance.metric = distance.metric, minimum = minimum,
    maximum = maximum, per = per, save.folder = save.folder,
    mode = mode, name.sort = name.sort, downsample = downsample,
    savePDFs = savePDFs, which.palette = which.palette,
    exclude.pctile = exclude.pctile, target.pctile = target.pctile,
    target.number = target.number, target.percent = target.percent)"
  } else {
    p.FLOWMAP <- "FLOWMAP(seed.X = seed.X, files = files, var.remove = var.remove,
    var.annotate = var.annotate, clustering.var = clustering.var,
    cluster.numbers = cluster.numbers, subsamples = subsamples,
    distance.metric = distance.metric, minimum = minimum,
    maximum = maximum, per = per, save.folder = save.folder,
    mode = mode, name.sort = name.sort, downsample = downsample,
    savePDFs = savePDFs, which.palette = which.palette)"
  }
  output.file <- "run_FLOWMAPR.R"
  output <- file(output.file)
  p.setup1 <- "rm(list = ls())"
  p.setup2 <- "library(FLOWMAPR)"
  if (downsample) {
    writeLines(c(p.setup1, p.setup2, p.files, p.mode, p.save.folder, p.per, p.minimum,
                 p.maximum, p.distance.metric, p.cluster.numbers,
                 p.var.annotate, p.var.remove, p.clustering.var,
                 p.seed.X, p.set.seed.X, p.subsamples, p.name.sort,
                 p.downsample, p.savePDFs, p.which.palette,
                 p.exclude.pctile, p.target.pctile, p.target.number,
                 p.target.percent, p.FLOWMAP), output)
  } else {
    writeLines(c(p.setup1, p.setup2, p.files, p.mode, p.save.folder, p.per, p.minimum,
                 p.maximum, p.distance.metric, p.cluster.numbers,
                 p.var.annotate, p.var.remove, p.clustering.var,
                 p.seed.X, p.set.seed.X, p.subsamples, p.name.sort,
                 p.downsample, p.savePDFs, p.which.palette, p.FLOWMAP), output)
  }
  close(output)
}
