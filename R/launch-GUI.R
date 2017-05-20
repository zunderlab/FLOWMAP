
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
  cur_dir <- getwd()
  distanceMetric <- c("manhattan", "euclidean")
  mode <- c("multi", "single", "one")
  colorpalette <- c("bluered", "jet", "CB")
  rawFCSdir <- tclVar(cur_dir)
  resDir <- tclVar(cur_dir)
  downsampleToggle <- tclVar("0")
  savePDFsToggle <- tclVar("0")
  subsampleNum <- tclVar("200")
  clusterNum <- tclVar("100")
  seedNum <- tclVar("1")
  edgepctNum <- tclVar("1")
  edgeNumMin <- tclVar("2")
  edgeNumMax <- tclVar("5")
  ret_var <- tclVar("")
  
  #  button functions
  reset_rawFCS_dir <- function() {
    rawFCS_dir <- ""
    rawFCS_dir <- tclvalue(tkchooseDirectory(title = "Choose your raw FCS files directory ..."))
    if (rawFCS_dir != "") {
      tclvalue(rawFCSdir) <- rawFCS_dir
      tclvalue(resDir) <- rawFCS_dir
    }
  }
  reset_res_dir <- function() {
    res_dir <- ""
    res_dir <- tclvalue(tkchooseDirectory(title = "Choose your result directory ..."))
    if (res_dir != "") {
      tclvalue(resDir) <- res_dir
    }
  }
  reset_fcs_data <- function() {
    fnames <- ""
    fnames <- tk_choose.files(default = paste(tclvalue(rawFCSdir), 
                                              "fcs", sep = .Platform$file.sep),
                              caption = "Select FCS files", 
                              multi = TRUE, filters = matrix(c("{fcs files}", "{.fcs}"), 1, 2), index = 1)
    if (length(fnames) >= 1) {
      fnames <- fnames[!(grepl(paste0(.Platform$file.sep, 
                                      "fcs$"), fnames))]  # remove empty .fcs files
    }
  }
  reset_num2null <- function() {
    tclvalue(fixedNum) <- "NULL"
  }
  reset_num2any <- function() {
    tclvalue(fixedNum) <- "1000"
  }
  rawFCSdir_help <- function() {
    tkmessageBox(title = "rawFCSdir", message = "The directory that contains the raw FCS files.", 
                 icon = "info", type = "ok")
  }
  distanceMetric_help <- function() {
    tkmessageBox(title = "distanceMetric", message = "Select the appropriate distance metric.", 
                 icon = "info", type = "ok")
  }
  mode_help <- function() {
    tkmessageBox(title = "mode", message = "Method of analysis via multiple FCS files or via one.", 
                 icon = "info", type = "ok")
  }
  color_help <- function() {
    tkmessageBox(title = "colorpalette", message = "Color palette to use for saved PDFs of output graph.", 
                 icon = "info", type = "ok")
  }
  resDir_help <- function() {
    tkmessageBox(title = "resDir", message = "The directory where result files will be generated.", 
                 icon = "info", type = "ok")
  }
  downsample_help = function() {
    tkmessageBox(title = "downsampleToggle", message = "Turn on/off SPADE downsampling.",
                 icon = "info", type = "ok")
  }
  savePDFs_help = function() {
    tkmessageBox(title = "savePDFsToggle", message = "Turn on/off saving PDFs of output graph.",
                 icon = "info", type = "ok")
  }
  subsampleNum_help = function() {
    tkmessageBox(title = "subsampleNum", message = "The number of cells to sample from each FCS file. If SPADE downsampling is selected, this specifies the target number.",
                 icon = "info", type = "ok")
  }
  clusterNum_help = function() {
    tkmessageBox(title = "clusterNum", message = "The number of clusters from each FCS file.",
                 icon = "info", type = "ok")
  }
  seedNum_help = function() {
    tkmessageBox(title = "seedNum", message = "The seed number for reproducible analysis.",
                 icon = "info", type = "ok")
  }
  edgepctNum_help = function() {
    tkmessageBox(title = "edgepctNum", message = "The edge percentile number",
                 icon = "info", type = "ok")
  }
  edgeNumMin_help = function() {
    tkmessageBox(title = "edgeNumMin", message = "The lower bound for number of edges during density-dependent edge drawing steps.",
                 icon = "info", type = "ok")
  }
  edgeNumMax_help = function() {
    tkmessageBox(title = "edgeNumMax", message = "The upper bound for number of edges during density-dependent edge drawing steps.",
                 icon = "info", type = "ok")
  }
  reset <- function() {
    tclvalue(rawFCSdir) <- cur_dir
    tclvalue(resDir) <- cur_dir
    tclvalue(distanceMetric) <- distanceMetric[1]
    tclvalue(mode) <- mode[2]
    tclvalue(colorpalette) <- colorpalette[1]
    tclvalue(downsampleToggle) <- "0"
    tclvalue(savePDFsToggle) <- "0"
    tclvalue(subsampleNum) <- "200"
    tclvalue(clusterNum) <- "100"
    tclvalue(seedNum) <- "1"
    tclvalue(edgepctNum) <- "1"
    tclvalue(edgeNumMin) <- "2"
    tclvalue(edgeNumMax) <- "5"
  }
  submit <- function() {
    has_error <- FALSE
    if (has_error == FALSE) {
      tclvalue(ret_var) <- "OK"
      tkdestroy(tt)
    }
  }
  quit <- function() {
    tkdestroy(tt)
  }
  # build the GUI
  # head line
  tt <- tktoplevel(borderwidth = 20)
  tkwm.title(tt, "FLOWMAPR")
  
  if (.Platform$OS.type == "windows") {
    box_length <- 63
  } else {
    box_length <- 55 
  }
  cell_width <- 3
  bt_width <- 8
  
  imgfile <- system.file("extdata", "help.gif", package = "cytofkit")
  image1 <- tclVar()
  tkimage.create("photo", image1, file = imgfile)
  image2 <- tclVar()
  tkimage.create("photo", image2)
  tcl(image2, "copy", image1, subsample = 6)
  
  # rawFCSdir
  rawFCSdir_label <- tklabel(tt, text = "Raw FCS Files Directory:")
  rawFCSdir_entry <- tkentry(tt, textvariable = rawFCSdir, width = box_length)
  rawFCSdir_button <- tkbutton(tt, text = " Choose... ", width = bt_width, command = reset_rawFCS_dir)
  rawFCSdir_hBut <- tkbutton(tt, image = image2, command = rawFCSdir_help)
  
  # resDir
  resDir_label <- tklabel(tt, text = "Result Directory:")
  resDir_entry <- tkentry(tt, textvariable = resDir, width = box_length)
  resDir_button <- tkbutton(tt, text = " Choose... ", width = bt_width, 
                            command = reset_res_dir)
  resDir_hBut <- tkbutton(tt, image = image2, command = resDir_help)
  
  # downsampleToggle
  downsample_label <- tklabel(tt, text = "SPADE Downsampling:")
  downsample_hBut <- tkbutton(tt, image = image2, command = downsample_help)
  downsample_rbuts <- tkframe(tt)
  tkpack(tklabel(downsample_rbuts, text = ""), side = "left")
  tkpack(tkcheckbutton(downsample_rbuts, variable = downsampleToggle), 
         side = "left")
  
  # savePDFsToggle
  savePDFs_label <- tklabel(tt, text = "Save Graph PDFs:")
  savePDFs_hBut <- tkbutton(tt, image = image2, command = savePDFs_help)
  savePDFs_rbuts <- tkframe(tt)
  tkpack(tklabel(savePDFs_rbuts, text = ""), side = "left")
  tkpack(tkcheckbutton(savePDFs_rbuts, variable = savePDFsToggle), 
         side = "left")
  
  # colorpalette Method
  color_label <- tklabel(tt, text = "Color Palette:")
  color_hBut <- tkbutton(tt, image = image2,
                         command = color_help)
  color_rbuts <- tkframe(tt)
  tkpack(tklabel(color_rbuts, text = ""), side = "left")
  tkpack(tkradiobutton(color_rbuts, text = colorpalette[1], 
                       variable = colorpalette, value = colorpalette[1]), side = "left")
  tkpack(tkradiobutton(color_rbuts, text = colorpalette[2],
                       variable = colorpalette, value = colorpalette[2]), side = "left")
  tkpack(tkradiobutton(color_rbuts, text = colorpalette[3],
                       variable = colorpalette, value = colorpalette[3]), side = "left")
  
  # FLOWMAPtypeMethod
  mode_label <- tklabel(tt, text = "FLOW-MAP Mode:")
  mode_hBut <- tkbutton(tt, image = image2,
                        command = mode_help)
  mode_rbuts <- tkframe(tt)
  tkpack(tklabel(mode_rbuts, text = ""), side = "left")
  tkpack(tkradiobutton(mode_rbuts, text = mode[1], 
                       variable = mode, value = mode[1]), side = "left")
  tkpack(tkradiobutton(mode_rbuts, text = mode[2],
                       variable = mode, value = mode[2]), side = "left")
  tkpack(tkradiobutton(mode_rbuts, text = mode[3],
                       variable = mode, value = mode[3]), side = "left")
  
  # distanceMetric
  distanceMetric_label <- tklabel(tt, text = "Distance Metric:")
  distanceMetric_hBut <- tkbutton(tt, image = image2, command = distanceMetric_help)
  distanceMetric_rbuts <- tkframe(tt)
  tkpack(tklabel(distanceMetric_rbuts, text = ""), side = "left")
  tkpack(tkradiobutton(distanceMetric_rbuts, text = distanceMetric[1], 
                       variable = distanceMetric, value = distanceMetric[1]), 
         side = "left")
  tkpack(tkradiobutton(distanceMetric_rbuts, text = distanceMetric[2],
                       variable = distanceMetric, value = distanceMetric[2]), 
         side = "left")
  
  # subsampleNum
  subsampleNum_label = tklabel(tt, text = "Subsample Number (Downsample Target Number):")
  subsampleNum_entry = tkentry(tt, textvariable = subsampleNum, width = 9)
  subsampleNum_hBut = tkbutton(tt, image = image2, command = subsampleNum_help)
  
  # clusterNum
  clusterNum_label = tklabel(tt, text = "Number of Clusters:")
  clusterNum_entry = tkentry(tt, textvariable = clusterNum, width = 9)
  clusterNum_hBut = tkbutton(tt, image = image2, command = clusterNum_help)
  
  # seedNum
  seedNum_label = tklabel(tt, text = "Seed Number:")
  seedNum_entry = tkentry(tt, textvariable = seedNum, width = 9)
  seedNum_hBut = tkbutton(tt, image = image2, command = seedNum_help)
  
  # edgepctNum
  edgepctNum_label = tklabel(tt, text = "Edge Percentile Number:")
  edgepctNum_entry = tkentry(tt, textvariable = edgepctNum, width = 9)
  edgepctNum_hBut = tkbutton(tt, image = image2, command = edgepctNum_help)
  
  # edgeNumMin
  edgeNumMin_label = tklabel(tt, text = "Minimum Number of Edges:")
  edgeNumMin_entry = tkentry(tt, textvariable = edgeNumMin, width = 9)
  edgeNumMin_hBut = tkbutton(tt, image = image2, command = edgeNumMin_help)
  
  # edgeNumMax
  edgeNumMax_label = tklabel(tt, text = "Maximum Number of Edges:")
  edgeNumMax_entry = tkentry(tt, textvariable = edgeNumMax, width = 9)
  edgeNumMax_hBut = tkbutton(tt, image = image2, command = edgeNumMax_help)
  
  # submit / reset / quit
  submit_button <- tkbutton(tt, text = "Submit", command = submit)
  reset_button <- tkbutton(tt, text = "Reset", command = reset)
  quit_button <- tkbutton(tt, text = "Quit", command = quit)
  
  # display GUI
  tkgrid(rawFCSdir_label, rawFCSdir_hBut, rawFCSdir_entry, rawFCSdir_button, 
         padx = cell_width)
  tkgrid.configure(rawFCSdir_label, rawFCSdir_entry, rawFCSdir_button, 
                   sticky = "e")
  tkgrid.configure(rawFCSdir_hBut, sticky = "e")
  
  tkgrid(resDir_label, resDir_hBut, resDir_entry, resDir_button, 
         padx = cell_width)
  tkgrid.configure(resDir_label, resDir_entry, resDir_button, 
                   sticky = "e")
  tkgrid.configure(resDir_hBut, sticky = "e")
  
  tkgrid(downsample_label, downsample_hBut, downsample_rbuts, 
         padx = cell_width)
  tkgrid.configure(downsample_label, sticky = "e")
  tkgrid.configure(downsample_hBut, sticky = "e")
  tkgrid.configure(downsample_rbuts, sticky = "w")
  
  tkgrid(savePDFs_label, savePDFs_hBut, savePDFs_rbuts, 
         padx = cell_width)
  tkgrid.configure(savePDFs_label, sticky = "e")
  tkgrid.configure(savePDFs_hBut, sticky = "e")
  tkgrid.configure(savePDFs_rbuts, sticky = "w")
  
  tkgrid(distanceMetric_label, distanceMetric_hBut, distanceMetric_rbuts, 
         padx = cell_width)
  tkgrid.configure(distanceMetric_label, sticky = "e")
  tkgrid.configure(distanceMetric_hBut, sticky = "e")
  tkgrid.configure(distanceMetric_rbuts, sticky = "w")
  
  tkgrid(mode_label, mode_hBut, mode_rbuts,
         padx = cell_width)
  tkgrid.configure(mode_label, sticky = "e")
  tkgrid.configure(mode_rbuts, sticky = "w")
  tkgrid.configure(mode_hBut, sticky = "e")
  
  tkgrid(color_label, color_hBut, color_rbuts,
         padx = cell_width)
  tkgrid.configure(color_label, sticky = "e")
  tkgrid.configure(color_rbuts, sticky = "w")
  tkgrid.configure(color_hBut, sticky = "e")
  
  tkgrid(subsampleNum_label, subsampleNum_hBut, subsampleNum_entry, padx = cell_width)
  tkgrid.configure(subsampleNum_label, sticky = "e")
  tkgrid.configure(subsampleNum_entry, sticky = "w")
  tkgrid.configure(subsampleNum_hBut, sticky = "e")
  
  tkgrid(clusterNum_label, clusterNum_hBut, clusterNum_entry, padx = cell_width)
  tkgrid.configure(clusterNum_label, sticky = "e")
  tkgrid.configure(clusterNum_entry, sticky = "w")
  tkgrid.configure(clusterNum_hBut, sticky = "e")
  
  tkgrid(seedNum_label, seedNum_hBut, seedNum_entry, padx = cell_width)
  tkgrid.configure(seedNum_label, sticky = "e")
  tkgrid.configure(seedNum_entry, sticky = "w")
  tkgrid.configure(seedNum_hBut, sticky = "e")

  tkgrid(edgepctNum_label, edgepctNum_hBut, edgepctNum_entry, padx = cell_width)
  tkgrid.configure(edgepctNum_label, sticky = "e")
  tkgrid.configure(edgepctNum_entry, sticky = "w")
  tkgrid.configure(edgepctNum_hBut, sticky = "e")
  
  tkgrid(edgeNumMin_label, edgeNumMin_hBut, edgeNumMin_entry, padx = cell_width)
  tkgrid.configure(edgeNumMin_label, sticky = "e")
  tkgrid.configure(edgeNumMin_entry, sticky = "w")
  tkgrid.configure(edgeNumMin_hBut, sticky = "e")
  
  tkgrid(edgeNumMax_label, edgeNumMax_hBut, edgeNumMax_entry, padx = cell_width)
  tkgrid.configure(edgeNumMax_label, sticky = "e")
  tkgrid.configure(edgeNumMax_entry, sticky = "w")
  tkgrid.configure(edgeNumMax_hBut, sticky = "e")
  
  tkgrid(tklabel(tt, text = "\n"), padx = cell_width)  # leave blank line
  
  tkgrid(reset_button, tklabel(tt, text = ""), submit_button, 
         quit_button, padx = cell_width)
  tkgrid.configure(reset_button, sticky = "e")
  tkgrid.configure(quit_button, sticky = "w")
  
  tkwait.window(tt)
  
  # Return parameters
  if (tclvalue(ret_var) != "OK") {
    okMessage <- "Analysis is cancelled."
  } else {
    inputs <- list()
    inputs[["mode"]] <- tclvalue(mode)
    inputs[["colorpalette"]] <- tclvalue(colorpalette)
    inputs[["downsampleToggle"]] <- tclvalue(downsampleToggle)
    inputs[["savePDFsToggle"]] <- tclvalue(savePDFsToggle)
    inputs[["subsampleNum"]] <- tclvalue(subsampleNum)
    inputs[["distanceMetric"]] <- tclvalue(distanceMetric)
    inputs[["clusterNum"]] <- tclvalue(clusterNum)
    inputs[["seedNum"]] <- tclvalue(seedNum)
    inputs[["edgepctNum"]] <- tclvalue(edgepctNum)
    inputs[["edgeMaxNum"]] <- tclvalue(edgeNumMax)
    inputs[["edgeminNum"]] <- tclvalue(edgeNumMin)
    inputs[["resultDir"]] <- tclvalue(resDir)
    
    print("inputs")
    print(inputs)
    
    globe_input <<- inputs
    globe_resdir <<- tclvalue(rawFCSdir)
    timeNow <- Sys.time()
    globe_resdir2 <<- tclvalue(resDir)
    timeNow <- gsub("[:]","-", timeNow) 
    okMessage <- paste0("Analysis Done, results are saved under ",
                        inputs[["resultDir"]])
  }
  runApp(appDir = file.path(system.file(package = "FLOWMAPR"), "shinyGUI"))
}

OpenDir <- function(dir = getwd()){
  if (.Platform['OS.type'] == "windows"){
    shell.exec(dir)
  } else {
    system(paste(Sys.getenv("R_BROWSER"), dir))
  }
}
