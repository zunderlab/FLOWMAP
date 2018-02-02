
CheckSettings <- function(mode, var.remove, var.annotate,
                          clustering.var, cluster.numbers,
                          distance.metric, minimum, maximum,
                          subsamples, which.palette) {
  check.whole.number <- function(x) { return (x %% 1 == 0) }
  check.pos.number <- function(x) { return (x > 0) }
  
  # need to add to checksettings
  # - maximum edges can’t be more than number of clusters in each time point
  # - clustering.var are all numerical values in files
  
  if (!dir.exists(save.folder)) {
    stop("save.folder does not exist!")
  }
  if (!(mode %in% c("single", "multi", "one", "one-special"))) {
    stop("mode is not a recognized value!")
  }
  if (!(distance.metric %in% c("manhattan", "euclidean"))) {
    stop("distance.metric is not a recognized value!")
  }
  if (!check.pos.number(minimum) || !check.whole.number(minimum)) {
    stop("minimum must be positive whole number!")
  }
  if (!check.pos.number(maximum) || !check.whole.number(maximum)) {
    stop("maximum must be positive whole number!")
  }
  if (minimum > maximum) {
    stop("maximum must be greater than or equal to minimum!")
  }
  if (is.numeric(cluster.numbers)) {
    if (!check.pos.number(cluster.numbers) || !check.whole.number(cluster.numbers)) {
      stop("cluster.numbers must be positive whole number!")
    }
  }
  if (is.numeric(subsamples)) {
    if (!check.pos.number(subsamples) || !check.whole.number(subsamples)) {
      stop("subsamples must be positive whole number!")
    }
  }
  if (is.numeric(subsamples) && is.numeric(cluster.numbers)) {
    if (subsamples < cluster.numbers) {
      stop("subsamples must be greater than cluster.numbers!")
    }
  }
  if (!(which.palette %in% c("bluered", "jet", "CB"))) {
    stop("which.palette is not a recognized value!")
  }
  if (!is.null(var.annotate)) {
    temp <- unname(unlist(var.annotate))
    temp2 <- setdiff(temp, var.remove)
    if (length(setdiff(clustering.var, temp2)) != 0) {
      stop("clustering.var contains a marker that is not in var.annotate!")
    }
    if (length(setdiff(var.remove, temp)) != 0) {
      warning("var.remove contains a marker that is not in var.annotate!")
    }
  }
  return()
}

CheckDownsampleSettings <- function(exclude.pctile, target.pctile, target.number, target.percent) {
  check.percent <- function(x) { return (x < 1 && x > 0) }
  if (is.null(exclude.pctile)) {
    stop("exclude.pctile not provided!")
  } else {
    if (!check.percent(exclude.pctile)) {
      stop("exclude.pctile must be value between 0 and 1!")
    }
  }
  if (is.null(target.pctile)) {
    stop("target.pctile not provided!")
  } else {
    if (!check.percent(target.pctile)) {
      stop("target.pctile must be value between 0 and 1!")
    }
  }
  if (is.null(target.number) && is.null(target.percent)) {
    warning("both target.number and target.percent are not provided!")
  }
  return() 
}
