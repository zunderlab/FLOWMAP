rm(list = ls(all = T))

library("MASS")
library("Matrix")
library("lattice")
library("flowCore")

setwd("/Users/mesako/Downloads/NewSyntheticData")

template.desc <- read.FCS("synthetic_use_description.fcs")
template.desc <- description(template.desc)
keep.these <- names(template.desc)[!grepl(names(template.desc), pattern = "flowCore")]
template.desc <- template.desc[keep.these[grepl(keep.these, pattern = "P")]]

fix.these <- names(template.desc)[grepl(names(template.desc), pattern = "R")]
fix.these <- fix.these[grepl("\\d", fix.these)]

epb = 500000 # events per branch
outnum = 50000 # events per fcs file

t1a <- mvrnorm(n = outnum, mu = c(1, 1),
               Sigma = matrix(c(0.2, 0.1, 0.1, 0.2), 2, 2))

t2a <- rbind(mvrnorm(n = (epb * 0.03), mu = c(1, 1),
                     Sigma = matrix(c(0.2, 0.1, 0.1, 0.2), 2, 2)), 
             mvrnorm(n = (epb * 0.3), mu = c(2, 2),
                     Sigma = matrix(c(0.2, 0.1, 0.1, 0.2), 2, 2)))
t2a <- t2a[sample(1:nrow(t2a), outnum, replace = FALSE), ]

t3b <- rbind(mvrnorm(n = (epb * 0.1), mu = c(2, 2),
                     Sigma = matrix(c(0.2, 0.1, 0.1, 0.2), 2, 2)), 
             mvrnorm(n = (epb * 0.4), mu = c(3.25, 3.25),
                     Sigma = matrix(c(0.2, 0.1, 0.1, 0.2), 2, 2)))
t3b <- t3b[sample(1:nrow(t3b), outnum, replace = FALSE), ]

t3c <- rbind(mvrnorm(n = (epb * 0.1), mu = c(1.5, 1.5),
                     Sigma = matrix(c(0.2, 0.1, 0.1, 0.2), 2, 2)), 
             mvrnorm(n = (epb * 0.4), mu = c(3.5, 2),
                     Sigma = matrix(c(0.5, 0, 0, 0.1), 2, 2)))
t3c <- t3c[sample(1:nrow(t3c), outnum, replace = FALSE), ]

t3d <- rbind(mvrnorm(n = (epb * 0.1), mu = c(1.5, 1.5),
                     Sigma = matrix(c(0.2, 0.1, 0.1, 0.2), 2, 2)), 
             mvrnorm(n = (epb * 0.4), mu = c(2, 3.5),
                     Sigma = matrix(c(0.1, 0, 0, 0.5), 2, 2)))
t3d <- t3d[sample(1:nrow(t3d), outnum, replace = FALSE), ]

t4b <- rbind(mvrnorm(n = (epb * 0.03), mu = c(2.75, 2.75),
                     Sigma = matrix(c(0.2, 0.1, 0.1, 0.2), 2, 2)), 
             mvrnorm(n = (epb * 0.05), mu = c(3.5, 3.5),
                     Sigma = matrix(c(0.2, 0.1, 0.1, 0.2), 2, 2)), 
             mvrnorm(n = (epb * 0.025), mu = c(4, 4),
                     Sigma = matrix(c(0.2, 0.1, 0.1, 0.2), 2, 2)))
t4b <- t4b[sample(1:nrow(t4b), outnum, replace = FALSE), ]

t4c <- rbind(mvrnorm(n = (epb * 0.04), mu = c(3, 2),
                     Sigma = matrix(c(0.5, 0, 0, 0.1), 2, 2)), 
             mvrnorm(n = (epb * 0.075), mu = c(4.5, 2),
                     Sigma = matrix(c(0.5, 0, 0, 0.1), 2, 2)))
t4c <- t4c[sample(1:nrow(t4c), outnum, replace = FALSE), ]

t4d <- rbind(mvrnorm(n = (epb * 0.04), mu = c(2, 3),
                     Sigma = matrix(c(0.1, 0, 0, 0.5), 2, 2)), 
             mvrnorm(n = (epb * 0.075), mu = c(2, 4.5),
                     Sigma = matrix(c(0.1, 0, 0, 0.5), 2, 2)))
t4d <- t4d[sample(1:nrow(t4d), outnum, replace = FALSE), ]

t5b <- rbind(mvrnorm(n = (epb * 0.04), mu = c(3.5, 3.5),
                     Sigma = matrix(c(0.2, 0.1, 0.1, 0.2), 2, 2)), 
             mvrnorm(n = (epb * 0.05), mu = c(3.75, 3.75),
                     Sigma = matrix(c(0.5, 0.4, 0.4, 0.5), 2, 2)), 
             mvrnorm(n = (epb * 0.2), mu = c(4.5, 4.5),
                     Sigma = matrix(c(0.3, 0.15, 0.15, 0.3), 2, 2)))
t5b <- t5b[sample(1:nrow(t5b), outnum, replace = FALSE), ]

t5c <- rbind(mvrnorm(n = (epb * 0.02), mu = c(3, 2),
                     Sigma = matrix(c(0.5, 0, 0, 0.1), 2, 2)), 
             mvrnorm(n = (epb * 0.075), mu = c(4, 2),
                     Sigma = matrix(c(0.5, 0, 0, 0.1), 2, 2)), 
             mvrnorm(n = (epb * 0.175), mu = c(5.5, 2),
                     Sigma = matrix(c(0.4, 0, 0, 0.1), 2, 2)))
t5c <- t5c[sample(1:nrow(t5c), outnum, replace = FALSE), ]

t5d <- rbind(mvrnorm(n = (epb * 0.02), mu = c(2, 3),
                     Sigma = matrix(c(0.1, 0, 0, 0.5), 2, 2)), 
             mvrnorm(n = (epb * 0.075), mu = c(2, 4),
                     Sigma = matrix(c(0.1, 0, 0, 0.5), 2, 2)), 
             mvrnorm(n = (epb * 0.175), mu = c(2, 5.5),
                     Sigma = matrix(c(0.1, 0, 0, 0.4), 2, 2)))
t5d <- t5d[sample(1:nrow(t5d), outnum, replace = FALSE), ]

t6b <- rbind(mvrnorm(n = (epb * 0.1), mu = c(4.4, 4.4),
                     Sigma = matrix(c(0.3, 0.15, 0.15, 0.3), 2, 2)), 
             mvrnorm(n = (epb * 0.2), mu = c(5.5, 5.5),
                     Sigma = matrix(c(0.3, 0.15, 0.15, 0.3), 2, 2)))
t6b <- t6b[sample(1:nrow(t6b), outnum, replace = FALSE), ]

t6c <- rbind(mvrnorm(n = (epb * 0.075), mu = c(5, 2),
                     Sigma = matrix(c(0.4, 0, 0, 0.1), 2, 2)), 
             mvrnorm(n = (epb * 0.175), mu = c(6.5, 2),
                     Sigma = matrix(c(0.4, 0, 0, 0.1), 2, 2)))
t6c <- t6c[sample(1:nrow(t6c), outnum, replace = FALSE), ]

t6d <- rbind(mvrnorm(n = (epb * 0.075), mu = c(2, 5),
                     Sigma = matrix(c(0.1, 0, 0, 0.4), 2, 2)), 
             mvrnorm(n = (epb * 0.175), mu = c(2, 6.5),
                     Sigma = matrix(c(0.1, 0, 0, 0.4), 2, 2)))
t6d <- t6d[sample(1:nrow(t6d), outnum, replace = FALSE), ]

palette <- colorRampPalette(c("blue", "blue", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red"))

outlist <- list(t1a, t2a, t3b, t4b, t5b, t6b, t3c, t4c, t5c, t6c, t3d, t4d, t5d, t6d)
outnames <- c("t1a", "t2a", "t3b", "t4b", "t5b", "t6b", "t3c", "t4c", "t5c", "t6c", "t3d", "t4d", "t5d", "t6d")

setwd("/Users/mesako/Downloads/NewSyntheticData/")

for (i in 1:length(outlist)) {
  d <- outlist[[i]]
  pointsize <- 1
  OVERLAY_COLOR <- densCols(x = d[, 1], y = d[, 2], colramp = palette)
  plot(x = d[, 1], y = d[, 2], col = OVERLAY_COLOR, pch = 15, cex = 0.2,
       ylim = c(-1, 10), xlim = c(-1, 10), main = outnames[i])
  d <- 5 * sinh(d)
  plot(x = d[, 1], y = d[, 2], col = OVERLAY_COLOR, pch = 15, cex = 0.2,
       ylim = c(-1, 10), xlim = c(-1, 10), main = paste(outnames[i], "_sinh"))
  colnames(d) <- c("marker1", "marker2")
  this.desc <- template.desc
  for (x in 1:length(fix.these)) {
    this.desc[[fix.these[x]]] <- max(d[, x])
    print(max(d[, x]))
  }
  print(this.desc)
  outframe <- flowFrame(exprs = d, description = this.desc)
  file.name <- paste(outnames[i], ".fcs", sep = "")
  write.FCS(outframe, filename = file.name)
}

# Make singleFLOWMAP data
names(outlist) <- outnames
for (i in 1:6) {
  use.these <- outnames[grepl(outnames, pattern = i)]
  use.these <- outlist[use.these]
  d <- c()
  for (x in 1:length(use.these)) {
    d <- rbind(d, use.these[[x]])
  }
  d <- d[sample(1:nrow(d), nrow(d), replace = FALSE), ]
  pointsize <- 1
  OVERLAY_COLOR <- densCols(x = d[, 1], y = d[, 2], colramp = palette)
  plot(x = d[, 1], y = d[, 2], col = OVERLAY_COLOR, pch = 15, cex = 0.2,
       ylim = c(-1, 10), xlim = c(-1, 10), main = paste("time", i, sep = ""))
  d <- 5 * sinh(d)
  plot(x = d[, 1], y = d[, 2], col = OVERLAY_COLOR, pch = 15, cex = 0.2,
       ylim = c(-1, 10), xlim = c(-1, 10), main = paste("time", i, "_sinh", sep = ""))
  colnames(d) <- c("marker1", "marker2")
  this.desc <- template.desc
  for (x in 1:length(fix.these)) {
    this.desc[[fix.these[x]]] <- max(d[, x])
    print(max(d[, x]))
  }
  print(this.desc)
  outframe <- flowFrame(exprs = d, description = this.desc)
  file.name <- paste("t", i, ".fcs", sep = "")
  write.FCS(outframe, filename = file.name)
}

# Make one-specialFLOWMAP data
all.conditions <- unique(unlist(lapply(strsplit(outnames, split = ""), tail, 1)))
for (i in 1:length(all.conditions)) {
  use.these <- outnames[grepl(outnames, pattern = all.conditions[i])]
  print(use.these)
  use.these <- outlist[use.these]
  d <- c()
  for (x in 1:length(use.these)) {
    d <- rbind(d, use.these[[x]])
  }
  d <- d[sample(1:nrow(d), nrow(d), replace = FALSE), ]
  pointsize <- 1
  OVERLAY_COLOR <- densCols(x = d[, 1], y = d[, 2], colramp = palette)
  plot(x = d[, 1], y = d[, 2], col = OVERLAY_COLOR, pch = 15, cex = 0.2,
       ylim = c(-1, 10), xlim = c(-1, 10), main = paste("condition", all.conditions[i], sep = ""))
  d <- 5 * sinh(d)
  plot(x = d[, 1], y = d[, 2], col = OVERLAY_COLOR, pch = 15, cex = 0.2,
       ylim = c(-1, 10), xlim = c(-1, 10), main = paste("condition", all.conditions[i], "_sinh", sep = ""))
  colnames(d) <- c("marker1", "marker2")
  this.desc <- template.desc
  for (x in 1:length(fix.these)) {
    this.desc[[fix.these[x]]] <- max(d[, x])
    print(max(d[, x]))
  }
  print(this.desc)
  outframe <- flowFrame(exprs = d, description = this.desc)
  file.name <- paste("condition", all.conditions[i], ".fcs", sep = "")
  write.FCS(outframe, filename = file.name)
}

# Make oneFLOWMAP data
d <- c()
for (i in 1:length(outlist)) {
  d <- rbind(d, (outlist)[[i]])
}
pointsize <- 1
OVERLAY_COLOR <- densCols(x = d[, 1], y = d[, 2], colramp = palette)
plot(x = d[, 1], y = d[, 2], col = OVERLAY_COLOR, pch = 15, cex = 0.2,
     ylim = c(-1, 10), xlim = c(-1, 10), main = paste("all_cells", sep = ""))
d <- 5 * sinh(d)
plot(x = d[, 1], y = d[, 2], col = OVERLAY_COLOR, pch = 15, cex = 0.2,
     ylim = c(-1, 10), xlim = c(-1, 10), main = paste("all_cells", "_sinh", sep = ""))

temp.d <- asinh(d / 5)
OVERLAY_COLOR <- densCols(x = temp.d[, 1], y = temp.d[, 2], colramp = palette)
plot(x = temp.d[, 1], y = temp.d[, 2], col = OVERLAY_COLOR, pch = 15, cex = 0.2,
     ylim = c(-1, 10), xlim = c(-1, 10), main = paste("all_cells", "_asinh", sep = ""))

d <- d[sample(1:nrow(d), nrow(d), replace = FALSE), ]
colnames(d) <- c("marker1", "marker2")
this.desc <- template.desc
for (x in 1:length(fix.these)) {
  this.desc[[fix.these[x]]] <- max(d[, x])
  print(max(d[, x]))
}
print(this.desc)
outframe <- flowFrame(exprs = d, description = this.desc)
file.name <- "SyntheticOne.fcs"
write.FCS(outframe, filename = file.name)


library(Biobase)
library(devtools)
install_github("nolanlab/cytofCore")
library(flowCore)

files <- "/Users/mesako/Downloads/NewSyntheticData/NewSingleFLOWMAP"
files <- list.files(files, pattern = "\\.fcs", recursive = TRUE, full.names = TRUE)
file.names <- basename(files)
file.names <- gsub(file.names, pattern = "\\.fcs", replacement = "")

FCS.frames <- lapply(files, read.FCS)
names(FCS.frames) <- file.names

new.temp.frames <- list()
for (i in 1:length(FCS.frames)) {
  temp.exprs <- exprs(FCS.frames[[i]])
  temp.exprs <- cbind(temp.exprs, Timepoint = rep(i, times = nrow(temp.exprs)))
  temp.params <- as(parameters(FCS.frames[[i]]), "data.frame")
  temp.params <- rbind(temp.params, c("Timepoint", "Timepoint", 0, i, i))
  rownames(temp.params)[nrow(temp.params)] <- "$P3"
  temp.desc <- description(FCS.frames[[i]])
  temp.desc[["$P3B"]] <- 32.0
  temp.desc[["$P3E"]] <- "0,0"
  temp.desc[["$P3N"]] <- "Timepoint"
  temp.desc[["$P3R"]] <- 0
  temp.desc[["$P3S"]] <- " "
  new.temp.frames[[i]] <- flowFrame(exprs = temp.exprs, parameters = AnnotatedDataFrame(temp.params),
                                    description = temp.desc)
  new.file.name <- gsub(files[i], pattern = ".fcs", replacement = "")
  new.file.name <- paste(new.file.name, "_with_time_param.fcs", sep = "")
  write.FCS(new.temp.frames[[i]], filename = new.file.name)
}

folder.name <- "/Users/mesako/Downloads/NewSyntheticData/NewSingleFLOWMAP/merged_files"
these.files <- list.files(folder.name, full.names = TRUE)
merged.file.name <- "Single_Synthetic_merged_by_time.fcs"
merged.exprs <- c()
for (i in 1:length(these.files)) {
  this.file <- these.files[i]
  this.file <- read.FCS(this.file)
  temp.exprs <- exprs(this.file)
  merged.exprs <- rbind(merged.exprs, temp.exprs)
}
merged.params <- as(parameters(this.file), "data.frame")
new.range <- ceiling(max(merged.exprs[, "Timepoint"]) - min(merged.exprs[, "Timepoint"]))
new.min.max <- c(min(merged.exprs[, "Timepoint"]), max(merged.exprs[, "Timepoint"]))
merged.params[nrow(merged.params), ] <- c("Timepoint", "Timepoint",
                                          new.range, new.min.max)
merged.params <- AnnotatedDataFrame(merged.params)
merged.desc <- description(this.file)
merged.desc[["$P3R"]] <- new.range
merged.file <- flowFrame(exprs = merged.exprs, parameters = merged.params,
                         description = merged.desc)
write.FCS(merged.file, filename = merged.file.name)  
rm(merged.file, merged.file.name, merged.desc, merged.params, merged.exprs)
