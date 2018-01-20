rm(list = ls())

# EXAMPLE FLOW-MAP RUN
# with SPADE downsampling instead of random subsampling

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
subsamples <- FALSE
name.sort <- FALSE
downsample <- TRUE
exclude.pctile <- 0.01
target.pctile <- 0.99
target.number <- 200
target.percent <- NULL
savePDFs <- TRUE
which.palette <- 'bluered'

FLOWMAP(mode = mode, files = files, var.remove = var.remove, var.annotate = var.annotate,
        clustering.var = clustering.var, cluster.numbers = cluster.numbers,
        distance.metric = distance.metric, minimum = minimum, maximum = maximum,
        save.folder = save.folder, subsamples = subsamples,
        name.sort = name.sort, downsample = downsample, seed.X = seed.X,
        savePDFs = savePDFs, which.palette = which.palette,
        exclude.pctile = exclude.pctile, target.pctile = target.pctile,
        target.number = target.number, target.percent = target.percent)

