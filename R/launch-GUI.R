#' GUI for launching shiny Application
#' 
#' This code was adapted from the cytofkit R package,
#' available here: https://github.com/JinmiaoChenLab/cytofkit/.
#' 
#' Specifically, the code in this file is adapted from cytofkit_GUI.R, but
#' modified to take user-specified inputs for FLOW-MAP runs and
#' then interface with the FLOWMAPR library back-end.
#' 
#' This code, authored and maintained by Jinmiao Chen and Hao Chen,
#' was used under the Artistic License 2.0 or Artistic-2.0 available
#' here: https://opensource.org/licenses/Artistic-2.0.
#' 
#' @param message A message to determine if open the shiny Application
#' @param dir Result direcroty.
#' 
#' @export
#' @examples
#' # LaunchGUI()

LaunchGUI <- function() {
  require(shiny)
  # parameter initialization
  current.dir <- getwd()
  distance.metric <- c("manhattan", "euclidean")
  mode <- c("multi", "single", "one")
  color.palette <- c("bluered", "jet", "CB")
  init.raw.FCS.dir <- tclVar(current.dir)
  init.result.dir <- tclVar(current.dir)
  downsample.toggle <- tclVar("0")
  savePDFs.toggle <- tclVar("0")
  subsample.num <- tclVar("200")
  cluster.num <- tclVar("100")
  seed.num <- tclVar("1")
  edge.pct.num <- tclVar("1")
  edge.num.min <- tclVar("2")
  edge.num.max <- tclVar("5")
  ret.var <- tclVar("")
  quit.var <- TRUE
  
  # button functions
  SetRawFCSDir <- function() {
    raw.FCS.dir <- tclvalue(tkchooseDirectory(title = "Choose your raw FCS files directory ..."))
    tclvalue(init.raw.FCS.dir) <- raw.FCS.dir
  }
  SetResultDir <- function() {
    result.dir <- tclvalue(tkchooseDirectory(title = "Choose your result directory ..."))
    tclvalue(init.result.dir) <- result.dir
  }
  RawFCSDirHelp <- function() {
    tkmessageBox(title = "raw.FCS.dir", message = "The directory that contains the raw FCS files.", 
                 icon = "info", type = "ok")
  }
  DistanceMetricHelp <- function() {
    tkmessageBox(title = "distance.metric", message = "Select the appropriate distance metric.", 
                 icon = "info", type = "ok")
  }
  ModeHelp <- function() {
    tkmessageBox(title = "mode", message = "Method of analysis via multiple FCS files or via one.", 
                 icon = "info", type = "ok")
  }
  ColorHelp <- function() {
    tkmessageBox(title = "color.palette", message = "Color palette to use for saved PDFs of output graph.", 
                 icon = "info", type = "ok")
  }
  ResultDirHelp <- function() {
    tkmessageBox(title = "result.dir", message = "The directory where result files will be generated.", 
                 icon = "info", type = "ok")
  }
  DownsampleHelp <- function() {
    tkmessageBox(title = "downsample.toggle", message = "Turn on/off SPADE downsampling.",
                 icon = "info", type = "ok")
  }
  SavePDFsHelp <- function() {
    tkmessageBox(title = "savePDFs.toggle", message = "Turn on/off saving PDFs of output graph.",
                 icon = "info", type = "ok")
  }
  SubsampleNumHelp <- function() {
    tkmessageBox(title = "subsample.num", message = "The number of cells to sample from each FCS file. If SPADE downsampling is selected, this specifies the target number.",
                 icon = "info", type = "ok")
  }
  ClusterNumHelp <- function() {
    tkmessageBox(title = "cluster.num", message = "The number of clusters from each FCS file.",
                 icon = "info", type = "ok")
  }
  SeedNumHelp <- function() {
    tkmessageBox(title = "seed.num", message = "The seed number for reproducible analysis.",
                 icon = "info", type = "ok")
  }
  EdgePctNumHelp <- function() {
    tkmessageBox(title = "edge.pct.num", message = "The edge percentile number",
                 icon = "info", type = "ok")
  }
  EdgeNumMinHelp <- function() {
    tkmessageBox(title = "edge.num.min", message = "The lower bound for number of edges during density-dependent edge drawing steps.",
                 icon = "info", type = "ok")
  }
  EdgeNumMaxHelp <- function() {
    tkmessageBox(title = "edge.num.max", message = "The upper bound for number of edges during density-dependent edge drawing steps.",
                 icon = "info", type = "ok")
  }
  Reset <- function() {
    tclvalue(raw.FCS.dir) <- current.dir
    tclvalue(result.dir) <- current.dir
    tclvalue(distance.metric) <- distance.metric[1]
    tclvalue(mode) <- mode[2]
    tclvalue(color.palette) <- color.palette[1]
    tclvalue(downsample.toggle) <- "0"
    tclvalue(savePDFs.toggle) <- "0"
    tclvalue(subsample.num) <- "200"
    tclvalue(cluster.num) <- "100"
    tclvalue(seed.num) <- "1"
    tclvalue(edge.pct.num) <- "1"
    tclvalue(edge.num.min) <- "2"
    tclvalue(edge.num.max) <- "5"
  }
  Submit <- function() {
    has.error <- FALSE
    if (has.error == FALSE) {
      quit.var <<- FALSE
      tclvalue(ret.var) <- "OK"
      tkdestroy(tt)
    }
  }
  Quit <- function() {
    tkdestroy(tt)
    stop("Exiting FLOWMAPR GUI.")
  }
  
  # build the GUI
  # head line
  tt <- tktoplevel(borderwidth = 20)
  tkwm.title(tt, "FLOWMAPR")
  if (.Platform$OS.type == "windows") {
    box.length <- 63
  } else {
    box.length <- 55 
  }
  cell.width <- 3
  bt.width <- 8
  
  imgfile <- system.file("extdata", "help.gif", package = "cytofkit")
  image1 <- tclVar()
  tkimage.create("photo", image1, file = imgfile)
  image2 <- tclVar()
  tkimage.create("photo", image2)
  tcl(image2, "copy", image1, subsample = 6)
  
  # raw.FCS.dir
  raw.FCS.dir.label <- tklabel(tt, text = "Raw FCS Files Directory:")
  raw.FCS.dir.entry <- tkentry(tt, textvariable = init.raw.FCS.dir, width = box.length)
  raw.FCS.dir.button <- tkbutton(tt, text = " Choose... ", width = bt.width, command = SetRawFCSDir)
  raw.FCS.dir.hBut <- tkbutton(tt, image = image2, command = RawFCSDirHelp)
  
  # result.dir
  result.dir.label <- tklabel(tt, text = "Result Directory:")
  result.dir.entry <- tkentry(tt, textvariable = init.result.dir, width = box.length)
  result.dir.button <- tkbutton(tt, text = " Choose... ", width = bt.width, 
                                command = SetResultDir)
  result.dir.hBut <- tkbutton(tt, image = image2, command = ResultDirHelp)
  
  # downsample.toggle
  downsample.label <- tklabel(tt, text = "SPADE Downsampling:")
  downsample.hBut <- tkbutton(tt, image = image2, command = DownsampleHelp)
  downsample.rbuts <- tkframe(tt)
  tkpack(tklabel(downsample.rbuts, text = ""), side = "left")
  tkpack(tkcheckbutton(downsample.rbuts, variable = downsample.toggle), 
         side = "left")
  
  # savePDFs.toggle
  savePDFs.label <- tklabel(tt, text = "Save Graph PDFs:")
  savePDFs.hBut <- tkbutton(tt, image = image2, command = SavePDFsHelp)
  savePDFs.rbuts <- tkframe(tt)
  tkpack(tklabel(savePDFs.rbuts, text = ""), side = "left")
  tkpack(tkcheckbutton(savePDFs.rbuts, variable = savePDFs.toggle), 
         side = "left")
  
  # color.palette Method
  color.label <- tklabel(tt, text = "Color Palette:")
  color.hBut <- tkbutton(tt, image = image2,
                         command = ColorHelp)
  color.rbuts <- tkframe(tt)
  tkpack(tklabel(color.rbuts, text = ""), side = "left")
  tkpack(tkradiobutton(color.rbuts, text = color.palette[1], 
                       variable = color.palette, value = color.palette[1]), side = "left")
  tkpack(tkradiobutton(color.rbuts, text = color.palette[2],
                       variable = color.palette, value = color.palette[2]), side = "left")
  tkpack(tkradiobutton(color.rbuts, text = color.palette[3],
                       variable = color.palette, value = color.palette[3]), side = "left")
  
  # FLOWMAPtypeMethod
  mode.label <- tklabel(tt, text = "FLOW-MAP Mode:")
  mode.hBut <- tkbutton(tt, image = image2,
                        command = ModeHelp)
  mode.rbuts <- tkframe(tt)
  tkpack(tklabel(mode.rbuts, text = ""), side = "left")
  tkpack(tkradiobutton(mode.rbuts, text = mode[1], 
                       variable = mode, value = mode[1]), side = "left")
  tkpack(tkradiobutton(mode.rbuts, text = mode[2],
                       variable = mode, value = mode[2]), side = "left")
  tkpack(tkradiobutton(mode.rbuts, text = mode[3],
                       variable = mode, value = mode[3]), side = "left")
  
  # distance.metric
  distance.metric.label <- tklabel(tt, text = "Distance Metric:")
  distance.metric.hBut <- tkbutton(tt, image = image2, command = DistanceMetricHelp)
  distance.metric.rbuts <- tkframe(tt)
  tkpack(tklabel(distance.metric.rbuts, text = ""), side = "left")
  tkpack(tkradiobutton(distance.metric.rbuts, text = distance.metric[1], 
                       variable = distance.metric, value = distance.metric[1]), 
         side = "left")
  tkpack(tkradiobutton(distance.metric.rbuts, text = distance.metric[2],
                       variable = distance.metric, value = distance.metric[2]), 
         side = "left")
  
  # subsample.num
  subsample.num.label <- tklabel(tt, text = "Subsample Number (Downsample Target Number):")
  subsample.num.entry <- tkentry(tt, textvariable = subsample.num, width = 9)
  subsample.num.hBut <- tkbutton(tt, image = image2, command = SubsampleNumHelp)
  
  # cluster.num
  cluster.num.label <- tklabel(tt, text = "Number of Clusters:")
  cluster.num.entry <- tkentry(tt, textvariable = cluster.num, width = 9)
  cluster.num.hBut <- tkbutton(tt, image = image2, command = ClusterNumHelp)
  
  # seed.num
  seed.num.label <- tklabel(tt, text = "Seed Number:")
  seed.num.entry <- tkentry(tt, textvariable = seed.num, width = 9)
  seed.num.hBut <- tkbutton(tt, image = image2, command = SeedNumHelp)
  
  # edge.pct.num
  edge.pct.num.label <- tklabel(tt, text = "Edge Percentile Number:")
  edge.pct.num.entry <- tkentry(tt, textvariable = edge.pct.num, width = 9)
  edge.pct.num.hBut <- tkbutton(tt, image = image2, command = EdgePctNumHelp)
  
  # edge.num.min
  edge.num.min.label <- tklabel(tt, text = "Minimum Number of Edges:")
  edge.num.min.entry <- tkentry(tt, textvariable = edge.num.min, width = 9)
  edge.num.min.hBut <- tkbutton(tt, image = image2, command = EdgeNumMinHelp)
  
  # edge.num.max
  edge.num.max.label <- tklabel(tt, text = "Maximum Number of Edges:")
  edge.num.max.entry <- tkentry(tt, textvariable = edge.num.max, width = 9)
  edge.num.max.hBut <- tkbutton(tt, image = image2, command = EdgeNumMaxHelp)
  
  # submit / reset / quit
  submit.button <- tkbutton(tt, text = "Submit", command = Submit)
  reset.button <- tkbutton(tt, text = "Reset", command = Reset)
  quit.button <- tkbutton(tt, text = "Quit", command = Quit)
  
  # display GUI
  tkgrid(raw.FCS.dir.label, raw.FCS.dir.hBut, raw.FCS.dir.entry, raw.FCS.dir.button, 
         padx = cell.width)
  tkgrid.configure(raw.FCS.dir.label, raw.FCS.dir.entry, raw.FCS.dir.button, 
                   sticky = "e")
  tkgrid.configure(raw.FCS.dir.hBut, sticky = "e")
  
  tkgrid(result.dir.label, result.dir.hBut, result.dir.entry, result.dir.button, 
         padx = cell.width)
  tkgrid.configure(result.dir.label, result.dir.entry, result.dir.button, 
                   sticky = "e")
  tkgrid.configure(result.dir.hBut, sticky = "e")
  
  tkgrid(downsample.label, downsample.hBut, downsample.rbuts, 
         padx = cell.width)
  tkgrid.configure(downsample.label, downsample.hBut, sticky = "e")
  tkgrid.configure(downsample.rbuts, sticky = "w")
  
  tkgrid(savePDFs.label, savePDFs.hBut, savePDFs.rbuts, 
         padx = cell.width)
  tkgrid.configure(savePDFs.label, savePDFs.hBut, sticky = "e")
  tkgrid.configure(savePDFs.rbuts, sticky = "w")
  
  tkgrid(distance.metric.label, distance.metric.hBut, distance.metric.rbuts, 
         padx = cell.width)
  tkgrid.configure(distance.metric.label, distance.metric.hBut, sticky = "e")
  tkgrid.configure(distance.metric.rbuts, sticky = "w")
  
  tkgrid(mode.label, mode.hBut, mode.rbuts,
         padx = cell.width)
  tkgrid.configure(mode.label, mode.hBut, sticky = "e")
  tkgrid.configure(mode.rbuts, sticky = "w")
  
  tkgrid(color.label, color.hBut, color.rbuts,
         padx = cell.width)
  tkgrid.configure(color.label, color.hBut, sticky = "e")
  tkgrid.configure(color.rbuts, sticky = "w")
  
  tkgrid(subsample.num.label, subsample.num.hBut, subsample.num.entry, padx = cell.width)
  tkgrid.configure(subsample.num.label, subsample.num.hBut, sticky = "e")
  tkgrid.configure(subsample.num.entry, sticky = "w")
  
  tkgrid(cluster.num.label, cluster.num.hBut, cluster.num.entry, padx = cell.width)
  tkgrid.configure(cluster.num.label, cluster.num.hBut, sticky = "e")
  tkgrid.configure(cluster.num.entry, sticky = "w")
  
  tkgrid(seed.num.label, seed.num.hBut, seed.num.entry, padx = cell.width)
  tkgrid.configure(seed.num.label, seed.num.hBut, sticky = "e")
  tkgrid.configure(seed.num.entry, sticky = "w")
  
  tkgrid(edge.pct.num.label, edge.pct.num.hBut, edge.pct.num.entry, padx = cell.width)
  tkgrid.configure(edge.pct.num.label, edge.pct.num.hBut, sticky = "e")
  tkgrid.configure(edge.pct.num.entry, sticky = "w")
  
  tkgrid(edge.num.min.label, edge.num.min.hBut, edge.num.min.entry, padx = cell.width)
  tkgrid.configure(edge.num.min.label, edge.num.min.hBut, sticky = "e")
  tkgrid.configure(edge.num.min.entry, sticky = "w")
  
  tkgrid(edge.num.max.label, edge.num.max.hBut, edge.num.max.entry, padx = cell.width)
  tkgrid.configure(edge.num.max.label, edge.num.max.hBut, sticky = "e")
  tkgrid.configure(edge.num.max.entry, sticky = "w")
  
  tkgrid(tklabel(tt, text = "\n"), padx = cell.width)  # leave blank line
  
  tkgrid(reset.button, tklabel(tt, text = ""), submit.button, 
         quit.button, padx = cell.width)
  tkgrid.configure(reset.button, sticky = "e")
  tkgrid.configure(quit.button, sticky = "w")
  
  tkwait.window(tt)
  
  # Return parameters
  if (tclvalue(ret.var) != "OK") {
    okMessage <- "Analysis is cancelled."
  } else {
    inputs <- list()
    inputs[["mode"]] <- tclvalue(mode)
    inputs[["color.palette"]] <- tclvalue(color.palette)
    inputs[["downsample.toggle"]] <- tclvalue(downsample.toggle)
    inputs[["savePDFs.toggle"]] <- tclvalue(savePDFs.toggle)
    inputs[["subsample.num"]] <- tclvalue(subsample.num)
    inputs[["distance.metric"]] <- tclvalue(distance.metric)
    inputs[["cluster.num"]] <- tclvalue(cluster.num)
    inputs[["seed.num"]] <- tclvalue(seed.num)
    inputs[["edge.pct.num"]] <- tclvalue(edge.pct.num)
    inputs[["edge.max.num"]] <- tclvalue(edge.num.max)
    inputs[["edge.min.num"]] <- tclvalue(edge.num.min)
    inputs[["quit"]] <- quit.var
    globe.inputs <<- inputs
    globe.raw.FCS.dir <<- tclvalue(init.raw.FCS.dir)
    timeNow <- Sys.time()
    globe.result.dir <<- tclvalue(init.result.dir)
    timeNow <- gsub("[:]","-", timeNow) 
    okMessage <- paste0("Analysis Done, results are saved under ",
                        inputs[["resultDir"]])
  }
  
  runApp(appDir = file.path(system.file(package = "FLOWMAPR"), "shinyGUI"))
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
  # runApp(appDir = file.path("C:\\Users\\Rohit\\Desktop\\FLOWMAP\\inst\\shinyGUI"))
=======
>>>>>>> origin/master
=======

>>>>>>> parent of cdf2a3d... quit button functionality
=======

>>>>>>> parent of cdf2a3d... quit button functionality
}
