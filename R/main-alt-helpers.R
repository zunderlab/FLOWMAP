
#' Restructure a dataframe into a list of dataframes
#'
#' \code{RestructureDF} generates the list of dataframes containing single-cell data
#' where each consecutive member is a dataframe from a different timepoint, from a single
#' dataframe where a column designates the timepoints and if relevant, the conditions, of
#' each individual cell. This function is a pre-processing step before starting a FLOWMAPR
#' analysis using FLOWMAPfromDF as it only accepts single-cell data in a certain format.
#'
#' @param df a single dataframe containing all single-cell data from all treatments
#' and conditions (if relevant) to be restructured
#' @param time.col.label Character specifying the channel name that should
#' be used to find the timepoint label for each cell, default is set to \code{"Time"}
#' @param condition.col.label Character specifying the channel name that should
#' be used to find the condition label for each cell, default is set to \code{NULL}
#' @return a list of dataframes containing the single-cell data or if \code{condition.col.label}
#' is provided (not \code{NULL}), a list of lists containing dataframes of single-cell data
#' @examples
#' \dontrun{df <- read.csv(file = "/single-cell-data.csv", header = TRUE, sep = ",")}
#' \dontrun{RestructureDF(df, time.col.label = "Timepoint", condition.col.label = "Treatment")}
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
    for (i in 1:length(all.times)) {
      list.of.df <- list()
      inds <- which(df[, time.col.label] == all.times[i])
      sub.df <- df[inds, ]
      cond.in.this.time <- unique(sub.df[, condition.col.label])
      for (j in 1:length(cond.in.this.time)) {
        cond.inds <- which(sub.df[, condition.col.label] == cond.in.this.time[j])
        list.of.df[[cond.in.this.time[j]]] <- sub.df[cond.inds, ]
      }
      multi.list.df[[i]] <- list.of.df
      rm(list.of.df)
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
      if (is.list(df[[i]]) && !is.data.frame(df[[i]])) {
        for (j in 1:length(df[[i]])) {
          rownames(fixed.df[[i]][[j]]) <- 1:nrow(df[[i]][[j]])
        }
      } else {
        rownames(fixed.df[[i]]) <- 1:nrow(df[[i]])
      }
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
      condition.labels <- c()
      this.time <- c()
      for (j in 1:length(multi.list.df[[i]])) {
        condition.labels <- c(condition.labels, as.character(unique(multi.list.df[[i]][[j]][, condition.col.label])))
        this.time <- c(this.time, as.character(unique(multi.list.df[[i]][[j]][, time.col.label])))
      }
      label.key[[i]] <- condition.labels
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
  print("label.key")
  print(label.key)
  return(label.key)
}

GetConditionsfromDF <- function(df.list, condition.col.label) {
  label.key.special <- c()
  new.df.list <- df.list
  for (i in 1:length(df.list[[1]])) {
    this.condition <- as.character(unique(df.list[[1]][[i]][, condition.col.label]))
    label.key.special <- c(label.key.special, this.condition)
    new.df.list[[1]][[i]][, condition.col.label] <- i
  }
  print("label.key.special")
  print(label.key.special)
  results <- list(new.df.list = new.df.list,
                  label.key.special = label.key.special)
  return(results)
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
