rm(list = ls())

# EXAMPLE FLOW-MAP RUN

library(FLOWMAPR)
files <- '/Users/mesako/Desktop/SingleFLOWMAP' # Change to where your FCS files are located
mode <- 'single'
save.folder <- '/Users/mesako/Desktop' # Change to where you want your results saved
per <- 1
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

FLOWMAP(seed.X = seed.X, files = files, var.remove = var.remove,
    var.annotate = var.annotate, clustering.var = clustering.var,
    cluster.numbers = cluster.numbers, subsamples = subsamples,
    distance.metric = distance.metric, minimum = minimum,
    maximum = maximum, per = per, save.folder = save.folder,
    mode = mode, name.sort = name.sort, downsample = downsample,
    savePDFs = savePDFs, which.palette = which.palette)
