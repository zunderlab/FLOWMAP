
#' @export
RestructureDF <- function(df, time.col.label = "Time", condition.col.label = NULL) {
  results <- list()
  all.times <- unique(df[, time.col.label])
  results[["all.times"]] <- all.times
  if (is.null(condition.col.label)) {
    list.of.df <- list()
    for (i in 1:length(all.times)) {
      inds <- which(df[, time.col.label] == all.times[i])
      list.of.df[[i]] <- df[inds, ]
      rm(inds)
    }
    new.df <- list.of.df
  } else {
    all.conditions <- unique(df[, condition.col.label])
    results[["all.conditions"]] <- all.conditions
    multi.list.df <- list()
    for (i in 1:length(all.conditions)) {
      list.of.df <- list()
      for (i in 1:length(all.times)) {
        inds <- which(df[, time.col.label] == all.times[i])
        list.of.df[[i]] <- df[inds, ]
        rm(inds)
      }
      multi.list.df[[i]] <- list.of.df
    }
    new.df <- multi.list.df
  }
  results[["new.df"]] <- new.df
  return(results) 
}

GetOrigTimesfromDF <- function(list.of.df, time.col.label = "Time", name.sort = TRUE) {
  if (class(list.of.df) == "list") {
    orig.times <- c()
    for (i in 1:length(list.of.df)) {
      orig.times <- c(orig.times, as.character(unique(list.of.df[[i]][, time.col.label])))
    }
  } else {
    stop("'list.of.df' is of unknown type!")
  }
  if (name.sort) {
    orig.times <- orig.times[order(as.numeric(orig.times))]
  }
  return(orig.times)
}

StripTimesfromDF <- function(df, time.col.label = "Time") {
  keep <- setdiff(colnames(df), time.col.label)
  fixed.df <- df[, keep]
  return(fixed.df)
}

StripTimesfromDFList <- function(list.of.df, time.col.label = "Time") {
  fixed.df.list <- list()
  for (i in 1:length(list.of.df)) {
    fixed.df.list[[i]] <- StripTimesfromDF(list.of.df[[i]], time.col.label)
  }
  return(fixed.df.list)
}

RemoveRowNames <- function(df) {
  fixed.df <- df
  if (is.list(df) && !is.data.frame(df)) {
    for (i in 1:length(df)) {
      rownames(fixed.df[[i]]) <- 1:nrow(df[[i]])
    }
  } else if (is.data.frame(df)) {
    rownames(fixed.df) <- 1:nrow(df)
  }
  return(fixed.df)
}

GetLabelKeyfromDF <- function(multi.list.df, time.col.label, condition.col.label) {
  if (class(multi.list.df) == "list" & class(multi.list.df[[1]]) == "list") {
    label.key <- list()
    all.times <- c()
    for (i in 1:length(multi.list.df)) {
      label.key[[i]] <- c()
      this.time <- c()
      for (j in 1:length(multi.list.df[[i]])) {
        label.key[[i]] <- c(label.key[[i]], as.character(unique(multi.list.df[[i]][[j]][, condition.col.label])))
        this.time <- c(this.time, as.character(unique(multi.list.df[[i]][[j]][, time.col.label])))
      }
      if (unique(this.time) > 1) {
        warning("More than one timepoint detected in sublist of multi.list.df! Using first timepoint found for label.")
        all.times <- c(all.times, this.time[[1]])
      } else {
        all.times <- c(all.times, unique(this.time))
      }
    }
  } else {
    stop("'multi.list.df' is of unknown type!")
  }
  names(label.key) <- all.times
  return(label.key)
}

CheckDFModeOne <- function(df) {
  # "df" variable should be one of the following:
  # dataframe with multiple rows/columns, all from one timepoint
  fail.flag <- TRUE
  if (class(df) == "data.frame") {
    if (dim(df)[1] >= 1 & dim(df)[2] >= 1) {
      fail.flag <- FALSE
    } else {
      stop("Dataframe does not contain at least one row and/or column!")
    } 
  } else {
    stop("Dataframe provided is not of type dataframe!")
  }
  return(fail.flag)
}

CheckDFModeSingle <- function(list.of.df) {
  fail.flag <- TRUE
  # "list.of.df" variable should be one of the following:
  # list of dataframes, each dataframe is from one timepoint
  # list goes in order of timepoints
  if (class(list.of.df) == "list") {
    for (i in 1:length(list.of.df)) {
      temp.fail.flag <- CheckDFModeOne(list.of.df[[i]])
      if (temp.fail.flag) {
        break
      }
    }
  } else {
    stop("List of dataframes provided is not of type list!")
  }
  fail.flag <- FALSE
  return(fail.flag)
}

CheckDFModeMulti <- function(multi.list.df) {
  fail.flag <- TRUE
  # "multi.list.df" variable should be one of the following:
  # list of list of dataframes, each sublist is from one timepoint
  # each dataframe in sublist is a different condition
  # in that specific timepoint
  # top-level list goes in order of timepoints
  if (class(multi.list.df) == "list") {
    for (i in 1:length(multi.list.df)) {
      temp.fail.flag <- CheckDFModeSingle(multi.list.df[[i]])
      if (temp.fail.flag) {
        break
      }
    }
  } else {
    stop("List of dataframes provided is not of type list!")
  }
  fail.flag <- FALSE
  return(fail.flag)
}
