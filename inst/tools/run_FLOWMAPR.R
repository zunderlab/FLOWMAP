rm(list = ls())

# EXAMPLE FLOW-MAP RUN

library(FLOWMAPR)
files <- '/Users/mesako/Desktop/SingleFLOWMAP' # Change to where your FCS files are located
mode <- 'single'
save.folder <- '/Users/mesako/Desktop' # Change to where you want your results saved
minimum <- 2
maximum <- 5
distance.metric <- 'manhattan'
cluster.numbers <- 100
var.annotate <- list('marker1' = 'marker1', 'marker2' = 'marker2')
var.remove <- c()
clustering.var <- c('marker1', 'marker2')
seed.X <- 1
set.seed(seed.X)
subsamples <- 200
name.sort <- FALSE
downsample <- FALSE
savePDFs <- TRUE
which.palette <- 'bluered'

FLOWMAP(mode = mode, files = files, var.remove = var.remove, var.annotate = var.annotate,
        clustering.var = clustering.var, cluster.numbers = cluster.numbers,
        distance.metric = distance.metric, minimum = minimum, maximum = maximum,
        save.folder = save.folder, subsamples = subsamples,
        name.sort = name.sort, downsample = downsample, seed.X = seed.X,
        savePDFs = savePDFs, which.palette = which.palette)
