Asinh <- function(value) {
  # Define a function for arcsinh transformation (based on definition from Nikolay).
  value <- value - 1
  for(i in 1:length(value)) {
    if ((value[i] < 0) | is.na(value[i])) {
      value[i] <- rnorm(1, mean = 0, sd = 0.01)
    }
    value <- value / 5
    value <- asinh(value) # value <- log(value + sqrt(value^2 + 1))  
    return(value)
  }
}


GetFCSNames <- function(folder, file.format, sort = TRUE) {
  # get FCS files
  fcs.files = list.files(path = folder, pattern = file.format,
                         recursive = FALSE, full.names = TRUE)
  if (sort) {
    # sort to organize by hour
    fcs.files <- sort(fcs.files)
  }
  return(fcs.files)
}

######################
#### FIX THIS???? ####
######################
GetMultiFCSNamesOLD <- function(folder, file.format, sort = TRUE) {
  # get FCS files
  subfolders <- list.files(folder, full.names = TRUE)
  list.of.condition.file.names <- list()
  for (folder in subfolders) {
    fcs.files = list.files(path = folder, pattern = file.format,
                           recursive = FALSE, full.names = TRUE)
    condition <- basename(folder)
    list.of.condition.file.names[[condition]] <- fcs.files
  }
  return(list.of.condition.file.names)
}


GetMultiFCSNames <- function(folder, file.format, sort = TRUE) {
  # get FCS files
  subfolders <- list.files(folder, full.names = TRUE)
  list.of.time.file.names <- list()
  for (folder in subfolders) {
    fcs.files = list.files(path = folder, pattern = file.format,
                           recursive = FALSE, full.names = TRUE)
    time <- basename(folder)
    list.of.time.file.names[[time]] <- fcs.files
  }
  return(list.of.time.file.names)
}


LoadCleanFCS <- function(fcs.file.names, channel.remove, channel.annotate,
                         subsamples = 10000, subsample.rand = FALSE,
                         transform = TRUE, scale = FALSE) {
  clean.fcs.files <- list()
  if (length(subsamples) == 1 & subsamples != FALSE) {
    cat("Subsampling all files to:", subsamples, "\n")
    subsample.new <- rep(subsamples, times = length(fcs.file.names))
    subsamples <- subsample.new
  }
  for (i in 1:length(fcs.file.names)) {
    # print FCS file that is being read
    current.file <- tail(strsplit(fcs.file.names[i],"/")[[1]],n=1)
    cat("Reading FCS file data from:", current.file, "\n")
    # store currently read FCS file
    if (!subsamples) {
      tmp.FCS1 <- read.FCS(fcs.file.names[i], transformation = "linearize", which.lines = NULL,
                           alter.names = FALSE, column.pattern = NULL, invert.pattern = FALSE,
                           decades = 0, ncdf = FALSE, min.limit = NULL, dataset = NULL,
                           emptyValue = TRUE)  
      tmp.FCS1 <- head(tmp.FCS1, n = unname(dim(tmp.FCS1)[1]))
    }
    else { 
      cat("Subsampling", current.file, "to", subsamples[i], "cells\n")
      if (subsample.rand) {
        subsamp.FCS1 <- read.FCS(fcs.file.names[i], transformation = "linearize",
                                 which.lines = subsamples[i],
                                 alter.names = FALSE, column.pattern = NULL,
                                 invert.pattern = FALSE, decades = 0, ncdf = FALSE,
                                 min.limit = NULL, dataset = NULL, emptyValue = TRUE)
        tmp.FCS1 <- head(subsamp.FCS1, n = subsamples[i])
      }
      else {
        full.FCS1 <- read.FCS(fcs.file.names[i], transformation = "linearize", which.lines = NULL,
                              alter.names = FALSE, column.pattern = NULL,
                              invert.pattern = FALSE, decades = 0, ncdf = FALSE,
                              min.limit = NULL, dataset = NULL, emptyValue = TRUE)
        tmp.FCS1 <- head(full.FCS1, n = subsamples[i])
      }
    }
    # rename variables with protein marker measured instead of metal channel
    cat("Fixing channel names from:", current.file, "\n")
    for (x in 1:length(colnames(tmp.FCS1))) {
      if (exists(colnames(tmp.FCS1)[x], where = channel.annotate)) {
        colnames(tmp.FCS1)[x] <- channel.annotate[[colnames(tmp.FCS1)[x]]]
      }
    }
    # remove unneeded variables
    cat("Removing unnecessary channel names from:", current.file, "\n")
    tmp.FCS2 <- subset(tmp.FCS1, select = colnames(tmp.FCS1)[!colnames(tmp.FCS1) %in% channel.remove])
    if (transform) {
      cat("Transforming data from:",current.file,"\n")
      tmp.FCS3 <- apply(tmp.FCS2, 2, Asinh) 
    }
    if (scale) {
      print("no scale implemented yet")
    }
    tmp.FCS4 <- as.data.frame(tmp.FCS3)
    clean.fcs.files[[i]] <- tmp.FCS4
    rm(tmp.FCS1, tmp.FCS2, tmp.FCS3, tmp.FCS4)
  }  
  return(clean.fcs.files)
}


LoadMultiCleanFCSOLD <- function(list.of.condition.file.names, channel.remove, channel.annotate,
                                 subsamples, subsample.rand = FALSE, transform = TRUE, scale = FALSE) {
  list.of.condition.clean.FCS.files <- list()
  subsamp.orig <- subsamples
  # print("subsamp.orig")
  # print(subsamp.orig)
  for (condition in names(list.of.condition.file.names)) {
    cat("condition is", condition, "\n")
    fcs.file.names <- list.of.condition.file.names[[condition]]
    cat("fcs.file.names are", basename(fcs.file.names), "\n")
    if (length(subsamp.orig) == 1 & subsamp.orig != FALSE) {
      cat("Subsampling all files to:", subsamp.orig, "\n")
      subsample.new <- rep(subsamp.orig, times = length(fcs.file.names))
      subsamples <- subsample.new
    } else {
      subsamples <- subsamp.orig[[condition]]
    }
    list.of.condition.clean.FCS.files[[condition]] <- LoadCleanFCS(fcs.file.names,
                                                                   channel.remove,
                                                                   channel.annotate,
                                                                   subsamples,
                                                                   subsample.rand = FALSE,
                                                                   transform, scale)
    for (i in 1:length(list.of.condition.clean.FCS.files[[condition]])) {
      Condition <- rep(condition, times = dim(list.of.condition.clean.FCS.files[[condition]][[i]])[1])
      list.of.condition.clean.FCS.files[[condition]][[i]] <- cbind(list.of.condition.clean.FCS.files[[condition]][[i]],
                                                                   Condition)
      rm(Condition)
    }
  }
  return(list.of.condition.clean.FCS.files)
}


LoadMultiCleanFCS <- function(list.of.time.file.names, channel.remove, channel.annotate,
                              subsamples, subsample.rand = FALSE, transform = TRUE, scale = FALSE) {
  list.of.time.clean.FCS.files <- list()
  subsamp.orig <- subsamples
  # print("subsamp.orig")
  # print(subsamp.orig)
  for (time in names(list.of.time.file.names)) {
    cat("time is", time, "\n")
    fcs.file.names <- list.of.time.file.names[[time]]
    cat("fcs.file.names are", basename(fcs.file.names), "\n")
    if (length(subsamp.orig) == 1 & subsamp.orig != FALSE) {
      cat("Subsampling all files to:", subsamp.orig, "\n")
      subsample.new <- rep(subsamp.orig, times = length(fcs.file.names))
      subsamples <- subsample.new
    } else {
      subsamples <- subsamp.orig[[time]]
    }
    list.of.time.clean.FCS.files[[time]] <- LoadCleanFCS(fcs.file.names,
                                                         channel.remove,
                                                         channel.annotate,
                                                         subsamples,
                                                         subsample.rand = FALSE,
                                                         transform, scale)
    f.names <- c() 
    for (i in 1:length(list.of.time.clean.FCS.files[[time]])) {
      Time <- rep(as.numeric(time), times = dim(list.of.time.clean.FCS.files[[time]][[i]])[1])
      list.of.time.clean.FCS.files[[time]][[i]] <- cbind.data.frame(list.of.time.clean.FCS.files[[time]][[i]],
                                                                    Time, stringsAsFactors = FALSE)
      this.name <- basename(fcs.file.names[i])
      this.name <- gsub(".fcs", "", this.name)
      if (grepl("-", this.name)) {
        this.name <- unlist(strsplit(this.name, split = "-"))[1]
      }
      if (grepl("\\.", this.name)) {
        this.name <- unlist(strsplit(this.name, split = "\\."))[1]
      }
      Condition <- rep(this.name, times = dim(list.of.time.clean.FCS.files[[time]][[i]])[1])
      # print("class(Condition)")
      # print(class(Condition))
      # print("class(Time)")
      # print(class(Time))
      list.of.time.clean.FCS.files[[time]][[i]] <- cbind.data.frame(list.of.time.clean.FCS.files[[time]][[i]],
                                                                    Condition, stringsAsFactors = FALSE)
      f.names <- c(f.names, this.name)
      rm(Time, Condition)
    }
    names(list.of.time.clean.FCS.files[[time]]) <- f.names
  }
  return(list.of.time.clean.FCS.files)
}



ConvertNumericLabelOLD <- function(list.of.clean.FCS.files.with.labels) {
  label.key <- list()
  fixed.list.FCS.files <- list()
  for (i in 1:length(list.of.clean.FCS.files.with.labels)) {
    label.key[[i]] <- names(list.of.clean.FCS.files.with.labels)[i]
    tmp.label.fcs <- list.of.clean.FCS.files.with.labels[[i]]
    for (n in 1:length(tmp.label.fcs)) {
      inds <- which(tmp.label.fcs[[n]] == label.key[[i]], arr.ind = TRUE)
      col.ind <- unname(inds[1, 2])
      tmp.label.fcs[[n]][, col.ind] <- i
    }
    fixed.list.FCS.files[[i]] <- tmp.label.fcs
  }
  names(fixed.list.FCS.files) <- names(list.of.clean.FCS.files.with.labels)
  return(list(fixed.list.FCS.files = fixed.list.FCS.files,
              label.key = label.key))
}

ConvertNumericLabel <- function(list.of.clean.FCS.files.with.labels) {
  label.key <- list()
  fixed.list.FCS.files <- list()
  for (i in 1:length(list.of.clean.FCS.files.with.labels)) {
    label.key[[i]] <- names(list.of.clean.FCS.files.with.labels[[i]])
    tmp.label.fcs <- list.of.clean.FCS.files.with.labels[[i]]
    for (n in 1:length(tmp.label.fcs)) {
      inds <- which(tmp.label.fcs[[n]] == label.key[[i]], arr.ind = TRUE)
      col.ind <- unname(inds[1, 2])
      # tmp.label.fcs[[n]][, col.ind] <- i
      tmp.label.fcs[[n]][, col.ind] <- as.numeric(n)
      # print(head())
      # for (i in 1:ncol(tmp.label.fcs[[n]])) {
      #   print(colnames(tmp.label.fcs[[n]])[i])
      #   print(head(tmp.label.fcs[[n]][, i]))
      #   print("class(tmp.label.fcs[[n]][, i])")
      #   print(class(tmp.label.fcs[[n]][, i]))
      # }
    }
    fixed.list.FCS.files[[i]] <- tmp.label.fcs
  }
  names(fixed.list.FCS.files) <- names(list.of.clean.FCS.files.with.labels)
  return(list(fixed.list.FCS.files = fixed.list.FCS.files,
              label.key = label.key))
}



ConvertCharacterLabel <- function(data.frame.with.numeric.labels, label.key) {
  # print("head(data.frame.with.numeric.labels)")
  # print(head(data.frame.with.numeric.labels))
  data.frame.with.character.labels <- data.frame.with.numeric.labels
  # print("label.key")
  # print(label.key)
  times <- unique(data.frame.with.numeric.labels[, "Time"])
  # cat("times are", times, "\n")
  for (t in times) {
    this.label <- label.key[[t]]
    this.ind <- which(data.frame.with.character.labels[, "Time"] == t)
    # cat("this.label is", this.label, "\n")
    for (i in length(this.label)) {
      # cat("i is", i, "and this.label[i] is", this.label[i], "\n")
      fix.ind <- which(data.frame.with.character.labels[, "Condition"] == i)
      fix.ind <- union(this.ind, fix.ind)
      # print("fix.ind")
      # print(fix.ind)
      data.frame.with.character.labels[fix.ind, "Condition"] <- this.label[i]
    }
  }
  # print("head(data.frame.with.character.labels)")
  # print(head(data.frame.with.character.labels))
  # print("tail(data.frame.with.character.labels)")
  # print(tail(data.frame.with.character.labels))
  return(data.frame.with.character.labels)
}


ConvertNumericLabel2 <- function(list.of.clean.FCS.files.with.labels) {
  label.key <- list()
  fixed.list.FCS.files <- list()
  for (i in 1:length(list.of.clean.FCS.files.with.labels)) {
    label.key[[i]] <- names(list.of.clean.FCS.files.with.labels)[i]
    # print("label.key[[i]]")
    # print(label.key[[i]])
    tmp.label.fcs <- list.of.clean.FCS.files.with.labels[[i]]
    # print("names(tmp.label.fcs)")
    # print(names(tmp.label.fcs))
    for (n in 1:length(tmp.label.fcs)) {
      inds <- which(tmp.label.fcs[[n]] == label.key[[i]], arr.ind = TRUE)
      col.ind <- unname(inds[1, 2])
      tmp.label.fcs[[n]][, col.ind] <- i
      # print("head(tmp.label.fcs[[n]])")
      # print(head(tmp.label.fcs[[n]]))
    }
    # print("names(tmp.label.fcs)")
    # print(names(tmp.label.fcs))
    fixed.list.FCS.files[[i]] <- tmp.label.fcs
    # print("names(fixed.list.FCS.files[[i]])")
    # print(names(fixed.list.FCS.files[[i]]))
  }
  names(fixed.list.FCS.files) <- names(list.of.clean.FCS.files.with.labels)
  # print("names(fixed.list.FCS.files)")
  # print(names(fixed.list.FCS.files))
  # print("label.key")
  # print(label.key)
  return(list(fixed.list.FCS.files = fixed.list.FCS.files,
              label.key = label.key))
}





# LoadCleanFCS.OLD <- function(fcs.file.names, channel.remove, channel.annotate,
#                              subsample = 10000, subsample.rand = FALSE,
#                              transform = TRUE, scale = FALSE) {
#   clean.fcs.files <- list()
#   if (length(subsample) == 1 & subsample != FALSE) {
#     cat("Subsampling all files to:", subsample, "\n")
#     subsample.new <- rep(subsample, times = length(fcs.file.names))
#     subsample <- subsample.new
#   }
#   for (i in 1:length(fcs.file.names)) {
#     # print FCS file that is being read
#     current.file <- tail(strsplit(fcs.file.names[i],"/")[[1]],n=1)
#     cat("Reading FCS file data from:", current.file, "\n")
#     # store currently read FCS file
#     if (!subsample) {
#       tmp.FCS1 <- read.FCS(fcs.file.names[i], transformation = "linearize", which.lines = NULL,
#                            alter.names = FALSE, column.pattern = NULL, invert.pattern = FALSE,
#                            decades = 0, ncdf = FALSE, min.limit = NULL, dataset = NULL,
#                            emptyValue = TRUE)  
#       tmp.FCS1 <- head(tmp.FCS1, n = unname(dim(tmp.FCS1)[1]))
#     }
#     else { 
#       cat("Subsampling", current.file, "to", subsample[i], "cells\n")
#       if (subsample.rand) {
#         subsamp.FCS1 <- read.FCS(fcs.file.names[i], transformation = "linearize",
#                                  which.lines = subsample[i],
#                                  alter.names = FALSE, column.pattern = NULL,
#                                  invert.pattern = FALSE, decades = 0, ncdf = FALSE,
#                                  min.limit = NULL, dataset = NULL, emptyValue = TRUE)
#         tmp.FCS1 <- head(subsamp.FCS1, n = subsample[i])
#       }
#       else {
#         full.FCS1 <- read.FCS(fcs.file.names[i], transformation = "linearize", which.lines = NULL,
#                               alter.names = FALSE, column.pattern = NULL,
#                               invert.pattern = FALSE, decades = 0, ncdf = FALSE,
#                               min.limit = NULL, dataset = NULL, emptyValue = TRUE)
#         tmp.FCS1 <- head(full.FCS1, n = subsample[i])
#       }
#     }
#     # rename variables with protein marker measured instead of metal channel
#     cat("Fixing channel names from:", current.file, "\n")
#     for (x in 1:length(colnames(tmp.FCS1))) {
#       if (exists(colnames(tmp.FCS1)[x], where = channel.annotate)) {
#         colnames(tmp.FCS1)[x] <- channel.annotate[[colnames(tmp.FCS1)[x]]]
#       }
#     }
#     # remove unneeded variables
#     cat("Removing unnecessary channel names from:", current.file, "\n")
#     tmp.FCS2 <- subset(tmp.FCS1, select = colnames(tmp.FCS1)[!colnames(tmp.FCS1) %in% channel.remove])
#     if (transform) {
#       cat("Transforming data from:",current.file,"\n")
#       tmp.FCS3 <- apply(tmp.FCS2, 2, Asinh) 
#     }
#     if (scale) {
#       print("no scale implemented yet")
#     }
#     tmp.FCS4 <- as.data.frame(tmp.FCS3)
#     clean.fcs.files[[i]] <- tmp.FCS4
#     rm(tmp.FCS1, tmp.FCS2, tmp.FCS3, tmp.FCS4)
#   }  
#   return(clean.fcs.files)
# }



# LoadMultiCleanFCS.OLD <- function(list.of.condition.file.names, channel.remove, channel.annotate,
#                                   subsamples = 10000, subsample.rand = FALSE, transform = TRUE, scale = FALSE) {
#   list.of.condition.clean.FCS.files <- list()
#   subsamp.orig <- subsample
#   print("subsamp.orig")
#   print(subsamp.orig)
#   for (condition in names(list.of.condition.file.names)) {
#     cat("condition is", condition, "\n")
#     fcs.file.names <- list.of.condition.file.names[[condition]]
#     cat("fcs.file.names are", basename(fcs.file.names), "\n")
#     if (length(subsamp.orig) == 1 & subsamp.orig != FALSE) {
#       cat("Subsampling all files to:", subsamp.orig, "\n")
#       subsample.new <- rep(subsamp.orig, times = length(fcs.file.names))
#       subsample <- subsample.new
#     } else {
#       print("boop")
#       subsample <- subsamp.orig[[condition]]
#     }
#     list.of.condition.clean.FCS.files[[condition]] <- LoadCleanFCS(fcs.file.names,
#                                                                    channel.remove,
#                                                                    channel.annotate,
#                                                                    subsample,
#                                                                    subsample.rand = FALSE,
#                                                                    transform, scale)
#     for (i in 1:length(list.of.condition.clean.FCS.files[[condition]])) {
#       Condition <- rep(condition, times = dim(list.of.condition.clean.FCS.files[[condition]][[i]])[1])
#       list.of.condition.clean.FCS.files[[condition]][[i]] <- cbind(list.of.condition.clean.FCS.files[[condition]][[i]],
#                                                                    Condition)
#       rm(Condition)
#     }
#   }
#   return(list.of.condition.clean.FCS.files)
# }
