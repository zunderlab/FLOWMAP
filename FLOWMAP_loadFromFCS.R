asinhNik <- function(value) {
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


getFCSNames <- function(folder, file_format, sort = TRUE) {
  # get FCS files
  fcs_files = list.files(path = folder, pattern = file_format,
                         recursive = FALSE, full.names = TRUE)
  if (sort) {
    # sort to organize by hour
    fcs_files <- sort(fcs_files)
  }
  return(fcs_files)
}


getMultiFCSNames <- function(listOfTreatments, folder, file_format, sort = TRUE) {
  # get FCS files
  listOfTreatFileNames <- list()
  for (treat in listOfTreatments) {
    folder_name <- paste(folder, treat, sep = "/")
    listOfTreatFileNames[[treat]] <- getFCSNames(folder_name, file_format, sort = TRUE)
  }
  return(listOfTreatFileNames)
}


loadCleanFCS <- function(fcs_file_names, channel_remove, channel_annotate,
                         subsample = 10000, subsampleRand = FALSE,
                         transform = TRUE, scale = FALSE) {
  clean_fcs_files <- list()
  for (i in 1:length(fcs_file_names)) {
    # print FCS file that is being read
    currentfile <- tail(strsplit(fcs_file_names[i],"/")[[1]],n=1)
    cat("Reading FCS file data from:", currentfile, "\n")
    # store currently read FCS file
    if (!subsample) {
      tmp_FCS1 <- read.FCS(fcs_file_names[i], transformation = "linearize", which.lines = NULL,
                           alter.names = FALSE, column.pattern = NULL, invert.pattern = FALSE,
                           decades = 0, ncdf = FALSE, min.limit = NULL, dataset = NULL,
                           emptyValue = TRUE)  
      tmp_FCS1 <- head(tmp_FCS1, n = unname(dim(tmp_FCS1)[1]))
    }
    else { 
      cat("Subsampling", currentfile, "to", subsample, "cells\n")
      if (subsampleRand) {
        subsamp_FCS1 <- read.FCS(fcs_file_names[i], transformation = "linearize", which.lines = subsample,
                                 alter.names = FALSE, column.pattern = NULL,
                                 invert.pattern = FALSE, decades = 0, ncdf = FALSE,
                                 min.limit = NULL, dataset = NULL, emptyValue = TRUE)
        tmp_FCS1 <- head(subsamp_FCS1, n = subsample)
      }
      else {
        full_FCS1 <- read.FCS(fcs_file_names[i], transformation = "linearize", which.lines = NULL,
                              alter.names = FALSE, column.pattern = NULL,
                              invert.pattern = FALSE, decades = 0, ncdf = FALSE,
                              min.limit = NULL, dataset = NULL, emptyValue = TRUE)
        tmp_FCS1 <- head(full_FCS1, n = subsample)
      }
    }
    # rename variables with protein marker measured instead of metal channel
    cat("Fixing channel names from:", currentfile, "\n")
    for (x in 1:length(colnames(tmp_FCS1))) {
      if (exists(colnames(tmp_FCS1)[x], where = channel_annotate)) {
        colnames(tmp_FCS1)[x] <- channel_annotate[[colnames(tmp_FCS1)[x]]]
      }
    }
    # remove unneeded variables
    cat("Removing unnecessary channel names from:", currentfile, "\n")
    tmp_FCS2 <- subset(tmp_FCS1, select = colnames(tmp_FCS1)[!colnames(tmp_FCS1) %in% channel_remove])
    if (transform) {
      cat("Transforming data from:",currentfile,"\n")
      tmp_FCS3 <- apply(tmp_FCS2, 2, asinhNik) 
    }
    if (scale) {
      print("no scale implemented yet")
    }
    tmp_FCS4 <- as.data.frame(tmp_FCS3)
    clean_fcs_files[[i]] <- tmp_FCS4
    rm(tmp_FCS1, tmp_FCS2, tmp_FCS3, tmp_FCS4)
  }  
  return(clean_fcs_files)
}


loadMultiCleanFCS <- function(listOfTreatFileNames, channel_remove, channel_annotate,
                              subsample = 10000, subsampleRand = FALSE, transform = TRUE, scale = FALSE) {
  listOfTreatCleanFCSFiles <- list()
  for (treat in names(listOfTreatFileNames)) {
    fcs_file_names <- listOfTreatFileNames[[treat]]
    listOfTreatCleanFCSFiles[[treat]] <- loadCleanFCS(fcs_file_names,
                                                      channel_remove,
                                                      channel_annotate,
                                                      subsample,
                                                      subsampleRand = FALSE,
                                                      transform, scale)
    for (i in 1:length(listOfTreatCleanFCSFiles[[treat]])) {
      Treat <- rep(treat, times = dim(listOfTreatCleanFCSFiles[[treat]][[i]])[1])
      listOfTreatCleanFCSFiles[[treat]][[i]] <- cbind(listOfTreatCleanFCSFiles[[treat]][[i]],
                                                 Treat)
#       cat("Treatment for file", treat, i, "is:", Treat, "\n")
#       print(head(Treat))
#       print(colnames(listOfTreatCleanFCSFiles[[treat]][[i]]))
      rm(Treat)
    }
  }
  return(listOfTreatCleanFCSFiles)
}

