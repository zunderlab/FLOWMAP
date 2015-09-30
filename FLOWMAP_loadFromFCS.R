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

# gateCisplatin <- function(fcs_file, VAR_ANNOTATE) {
#   cisplatin_channel <- names(VAR_ANNOTATE)[[which(VAR_ANNOTATE == "Cisplatin")]]
#   fcs_file[, cisplatin_channel]
#   return(gated_fcs_file)
# }

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


upsampleFCS <- function(clean_fcs_files, channel_upsample,
                        percent_upsample_limit, subsample) {
  # add in low, mid, high functionality for which population to upsample
  # somehow sample so it's a uniform distribution of that marker?
  upsampled_fcs_files <- list()
  all <- data.frame()
  for (i in 1:length(clean_fcs_files)) {
    all <- rbind(all, clean_fcs_files[[i]])
  }
  dist <- hist(all[, channel_upsample], xlab = channel_upsample, main = NULL)
  local_minima <- which(diff(sign(diff(dist$counts))) == 2) + 1
  if (length(local_minima) > 1) {
    local_minima <- sort(local_minima)
    local_maxima <- which(diff(sign(diff(dist$counts))) == -2) + 1
    local_maxima <- sort(local_maxima)
    last <- tail(local_maxima, n = 1)
    min <- local_minima[which(local_minima < last)]
  }    
  else {
    min <- local_minima
  }
  # min <- local_minima[which.min(dist$counts[local_minima])]
  leeway <- round(length(dist$breaks) / 10)
  leeway <- max(leeway, 1)
  window <- c(dist$breaks[(min - leeway)], dist$breaks[(min + leeway)])
  upsample_num <- sum(dist$counts[(min - leeway):(min + leeway)])  
  cat("Number of cells in upsampled population across files is", upsample_num, "\n")
  #   window <- dist$breaks[window_breaks]
  #   max <- local_maxima[which.max(dist$counts[local_maxima])]
  #   max_val <- local_maxima[which.max(dist$counts[local_maxima])]
  for (i in 1:length(clean_fcs_files)) {
    currentfile <- clean_fcs_files[[i]]
    cat("Upsampling on: File", i, "\n")
    upsamplepop <- currentfile[which(currentfile[, channel_upsample] < window[2]),]
    upsamplepop <- upsamplepop[which(upsamplepop[, channel_upsample] > window[1]),]
    if (dim(upsamplepop)[1] > 0) {
      upsample_limit <- percent_upsample_limit * subsample
      upsample_limit <- min(upsample_limit, dim(upsamplepop)[1])
      upsamplepop <- head(upsamplepop, n = upsample_limit)
      # subsample <- min(dim(currentfile)[1], subsample)
      #           if (dim(upsamplepop)[1] > upsample_limit) {
      #             upsamplepop <- head(upsamplepop, n = upsample_limit)
      #           }
      allminuspop <- currentfile[which(currentfile[, channel_upsample] > window[2]),]
      allminuspop <- rbind(allminuspop, currentfile[which(currentfile[, channel_upsample] < window[1]),])
      tmp_FCS <- head(allminuspop, n = (subsample - upsample_limit))
      tmp_FCS <- rbind(tmp_FCS, upsamplepop)
    }
    else {
      tmp_FCS <- head(currentfile, n = subsample)
    }
    upsampled_fcs_files[[i]] <- tmp_FCS
    rm(tmp_FCS, allminuspop, upsamplepop)
  }   
  cat("Window of values for upsampled variable", channel_upsample, ": ", window, "\n")
  return(upsampled_fcs_files)
}


upsampleMultiFCS <- function(listOfTreatCleanFCSFiles, channel_upsample,
                             percent_upsample_limit, subsample) {
  listOfTreatUpsampleFCSFiles <- list()
  for (treat in names(listOfTreatCleanFCSFiles)) {
    clean_fcs_files <- listOfTreatCleanFCSFiles[[treat]]
    listOfTreatUpsampleFCSFiles[[treat]] <- upsampleFCS(clean_fcs_files,
                                                        channel_upsample,
                                                        percent_upsample_limit,
                                                        subsample)
  }
  return(listOfTreatUpsampleFCSFiles)
}
