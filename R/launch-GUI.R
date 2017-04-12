
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
#' # launch_GUI()

launch_GUI <- function() {
  # parameter initialization
  cur_dir <- getwd()
  distanceMetrics <- c("manhattan", "euclidean")
  # multiSingles <- c("multi", "single", "one")
  multiSingles <- c("multi", "single")
  
  rawFCSdir <- tclVar(cur_dir)
  resDir <- tclVar(cur_dir)
  projectName <- tclVar("FLOWMAP")
  distanceMetric <- tclVar("manhattan")
  multiSingle = tclVar("single")
  subsampleNum = tclVar("200")
  clusterNum = tclVar("100")
  seedNum = tclVar("1")
  edgepctNum = tclVar("1")
  edgeNumMin = tclVar("2")
  edgeNumMax = tclVar("5")
  ret_var <- tclVar("")
  
  #  button functions
  reset_rawFCS_dir <- function() {
    rawFCS_dir <- ""
    rawFCS_dir <- tclvalue(tkchooseDirectory(title = "Choose your rawFCS dircetory ..."))
    if (rawFCS_dir != "") {
      tclvalue(rawFCSdir) <- rawFCS_dir
      tclvalue(resDir) <- rawFCS_dir
    }
  }
  
  reset_res_dir <- function() {
    res_dir <- ""
    res_dir <- tclvalue(tkchooseDirectory(title = "Choose your result dircetory ..."))
    if (res_dir != "") {
      tclvalue(resDir) <- res_dir
    }
  }
  
  reset_fcs_data <- function() {
    fnames <- ""
    fnames <- tk_choose.files(default = paste(tclvalue(rawFCSdir), 
                                              "fcs", sep = .Platform$file.sep), caption = "Select FCS files", 
                              multi = TRUE, filters = matrix(c("{fcs files}", "{.fcs}"), 
                                                             1, 2), index = 1)
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
    tkmessageBox(title = "rawFCSdir", message = "The directory that contains fcs files.", 
                 icon = "info", type = "ok")
  }
  
  projectName_help <- function() {
    tkmessageBox(title = "projectName", message = "A prefix that will be added to the names of result files.", 
                 icon = "info", type = "ok")
  }
  
  distanceMetric_help <- function() {
    tkmessageBox(title = "distanceMetric", message = "Select the appropriate distance metric.", 
                 icon = "info", type = "ok")
  }
  
  multiSingle_help <- function() {
    tkmessageBox(title = "multiSingle", message = "Method of analysis via multiple FCS files or via one.", 
                 icon = "info", type = "ok")
  }
  
  resDir_help <- function() {
    tkmessageBox(title = "resDir", message = "The directory where result files will be generated", 
                 icon = "info", type = "ok")
  }
  
  subsampleNum_help = function() {
    tkmessageBox(title = "subsampleNum", message = "The subsampling number",
                 icon = "info", type = "ok")
  }
  
  clusterNum_help = function() {
    tkmessageBox(title = "clusterNum", message = "The clustering number",
                 icon = "info", type = "ok")
  }
  
  seedNum_help = function() {
    tkmessageBox(title = "seedNum", message = "The seeding number",
                 icon = "info", type = "ok")
  }
  
  edgepctNum_help = function() {
    tkmessageBox(title = "edgepctNum", message = "The edge percentile number",
                 icon = "info", type = "ok")
  }
  
  edgeNumMin_help = function() {
    tkmessageBox(title = "edgeNumMin", message = "The lower bound for number of edges",
                 icon = "info", type = "ok")
  }
  
  edgeNumMax_help = function() {
    tkmessageBox(title = "edgeNumMax", message = "The upper bound for number of edges",
                 icon = "info", type = "ok")
  }
  
  reset <- function() {
    tclvalue(rawFCSdir) = cur_dir
    tclvalue(resDir) = cur_dir
    tclvalue(projectName) = "FLOWMAP"
    tclvalue(distanceMetric) = distanceMetrics[1]
    tclvalue(multiSingle) = multiSingles[2]
    tclvalue(subsampleNum) = "200"
    tclvalue(clusterNum) = "100"
    tclvalue(seedNum) = "1"
    tclvalue(edgepctNum) = "1"
    tclvalue(edgeNumMin) = "2"
    tclvalue(edgeNumMax) = "5"
  }
  
  submit <- function() {
    has_error = FALSE
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
  tkwm.title(tt, "FLOWMAP")
  
  if(.Platform$OS.type == "windows"){
    box_length <- 63
  }else{
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
  rawFCSdir_label <- tklabel(tt, text = "Raw FCS Directory :")
  rawFCSdir_entry <- tkentry(tt, textvariable = rawFCSdir, width = box_length)
  rawFCSdir_button <- tkbutton(tt, text = " Choose... ", width = bt_width, command = reset_rawFCS_dir)
  rawFCSdir_hBut <- tkbutton(tt, image = image2, command = rawFCSdir_help)
  
  # subsampleNum
  subsampleNum_label = tklabel(tt, text = "Subsample Number :")
  subsampleNum_entry = tkentry(tt, textvariable = subsampleNum, width = 9)
  subsampleNum_hBut = tkbutton(tt, image = image2, command = subsampleNum_help)
  
  # clusterNum
  clusterNum_label = tklabel(tt, text = "Clustering Number :")
  clusterNum_entry = tkentry(tt, textvariable = clusterNum, width = 9)
  clusterNum_hBut = tkbutton(tt, image = image2, command = clusterNum_help)
  
  # seedNum
  seedNum_label = tklabel(tt, text = "Seed number :")
  seedNum_entry = tkentry(tt, textvariable = seedNum, width = 9)
  seedNum_hBut = tkbutton(tt, image = image2, command = seedNum_help)
  
  # edgepctNum
  edgepctNum_label = tklabel(tt, text = "Edge Percentile Number :")
  edgepctNum_entry = tkentry(tt, textvariable = edgepctNum, width = 9)
  edgepctNum_hBut = tkbutton(tt, image = image2, command = edgepctNum_help)
  
  # edgeNumMin
  edgeNumMin_label = tklabel(tt, text = "Minimum Number of Edges :")
  edgeNumMin_entry = tkentry(tt, textvariable = edgeNumMin, width = 9)
  edgeNumMin_hBut = tkbutton(tt, image = image2, command = edgeNumMin_help)
  
  # edgeNumMax
  edgeNumMax_label = tklabel(tt, text = "Maximum Number of Edges :")
  edgeNumMax_entry = tkentry(tt, textvariable = edgeNumMax, width = 9)
  edgeNumMax_hBut = tkbutton(tt, image = image2, command = edgeNumMax_help)
  
  # resDir
  resDir_label <- tklabel(tt, text = "Result Directory :")
  resDir_entry <- tkentry(tt, textvariable = resDir, width = box_length)
  resDir_button <- tkbutton(tt, text = " Choose... ", width = bt_width, 
                            command = reset_res_dir)
  resDir_hBut <- tkbutton(tt, image = image2, command = resDir_help)
  
  # projectName
  projectName_label <- tklabel(tt, text = "Project Name :")
  projectName_entry <- tkentry(tt, textvariable = projectName, width = box_length)
  projectName_hBut <- tkbutton(tt, image = image2, command = projectName_help)
  
  # distanceMetric
  distanceMetric_label <- tklabel(tt, text = "Distance Metric :")
  distanceMetric_hBut <- tkbutton(tt, image = image2, command = distanceMetric_help)
  distanceMetric_rbuts <- tkframe(tt)
  tkpack(tklabel(distanceMetric_rbuts, text = ""), side = "left")
  tkpack(tkradiobutton(distanceMetric_rbuts, text = distanceMetrics[1], 
                       variable = distanceMetric, value = distanceMetrics[1]), 
         side = "left")
  tkpack(tkradiobutton(distanceMetric_rbuts, text = distanceMetrics[2],
                       variable = distanceMetric, value = distanceMetrics[2]), 
         side = "left")
  
  # transformMethod
  multiSingle_label <- tklabel(tt, text = "FLOW-MAP method: ")
  multiSingle_hBut <- tkbutton(tt, image = image2,
                               command = multiSingle_help)
  multiSingle_rbuts <- tkframe(tt)
  tkpack(tklabel(multiSingle_rbuts, text = ""), side = "left")
  tkpack(tkradiobutton(multiSingle_rbuts, text = multiSingles[1], 
                       variable = multiSingle, value = multiSingles[1]), side = "left")
  tkpack(tkradiobutton(multiSingle_rbuts, text = multiSingles[2],
                       variable = multiSingle, value = multiSingles[2]), side = "left")
  
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
  
  tkgrid(projectName_label, projectName_hBut, projectName_entry, padx = cell_width)
  tkgrid.configure(projectName_label, projectName_entry, sticky = "e")
  tkgrid.configure(projectName_hBut, sticky = "e")
  
  tkgrid(distanceMetric_label, distanceMetric_hBut, distanceMetric_rbuts, 
         padx = cell_width)
  tkgrid.configure(distanceMetric_label, sticky = "e")
  tkgrid.configure(distanceMetric_hBut, sticky = "e")
  tkgrid.configure(distanceMetric_rbuts, sticky = "w")
  
  tkgrid(multiSingle_label, multiSingle_hBut, multiSingle_rbuts,
         padx = cell_width)
  tkgrid.configure(multiSingle_label, sticky = "e")
  tkgrid.configure(multiSingle_rbuts, sticky = "w")
  tkgrid.configure(multiSingle_hBut, sticky = "e")
  
  tkgrid(subsampleNum_label, subsampleNum_hBut, subsampleNum_entry, padx = cell_width)
  tkgrid.configure(subsampleNum_hBut, sticky = "w")
  tkgrid.configure(subsampleNum_label, subsampleNum_entry, sticky = "w")
  
  tkgrid(clusterNum_label, clusterNum_hBut, clusterNum_entry, padx = cell_width)
  tkgrid.configure(clusterNum_hBut, sticky = "w")
  tkgrid.configure(clusterNum_label, clusterNum_entry, sticky = "w")
  
  tkgrid(seedNum_label, seedNum_hBut, seedNum_entry, padx = cell_width)
  tkgrid.configure(seedNum_hBut, sticky = "w")
  tkgrid.configure(seedNum_label, seedNum_entry, sticky = "w")
  
  tkgrid(edgepctNum_label, edgepctNum_hBut, edgepctNum_entry, padx = cell_width)
  tkgrid.configure(edgepctNum_hBut, sticky = "w")
  tkgrid.configure(edgepctNum_label, edgepctNum_entry, sticky = "w")
  
  tkgrid(edgeNumMin_label, edgeNumMin_hBut, edgeNumMin_entry, padx = cell_width)
  tkgrid.configure(edgeNumMin_hBut, sticky = "w")
  tkgrid.configure(edgeNumMin_label, edgeNumMin_entry, sticky = "w")
  
  tkgrid(edgeNumMax_label, edgeNumMax_hBut, edgeNumMax_entry, padx = cell_width)
  tkgrid.configure(edgeNumMax_hBut, sticky = "w")
  tkgrid.configure(edgeNumMax_label, edgeNumMax_entry, sticky = "w")
  
  tkgrid(tklabel(tt, text = "\n"), padx = cell_width)  # leave blank line
  
  tkgrid(reset_button, tklabel(tt, text = ""), submit_button, 
         quit_button, padx = cell_width)
  tkgrid.configure(reset_button, sticky = "e")
  tkgrid.configure(quit_button, sticky = "w")
  
  tkwait.window(tt)
  
  # Return parameters
  if (tclvalue(ret_var) != "OK") {
    okMessage <- "Analysis is cancelled."
  }else{
    inputs <- list()
    inputs[["multiSingle"]] = tclvalue(multiSingle)
    inputs[["subsampleNum"]] = tclvalue(subsampleNum)
    inputs[["distanceMetric"]] = tclvalue(distanceMetric)
    inputs[["clusterNum"]] = tclvalue(clusterNum)
    inputs[["seedNum"]] = tclvalue(seedNum)
    inputs[["edgepctNum"]] = tclvalue(edgepctNum)
    inputs[["edgeMaxNum"]] = tclvalue(edgeNumMax)
    inputs[["edgeminNum"]] = tclvalue(edgeNumMin)
    inputs[["projectName"]] <- tclvalue(projectName)
    inputs[["resultDir"]] <- tclvalue(resDir)
    
    globe_input <<- inputs
    globe_resdir <<- tclvalue(rawFCSdir)
    globe_proj <<- tclvalue(projectName)
    timeNow = Sys.time()
    globe_resdir2 <<- tclvalue(resDir)
    
    timeNow = gsub("[:]","-", timeNow) 
    
    write.csv(inputs, paste0(inputs[["resultDir"]], "/", inputs[["projectName"]], timeNow, ".csv"))
    
    okMessage <- paste0("Analysis Done, results are saved under ",
                        inputs[["resultDir"]])
  }

  runApp(appDir = file.path(system.file(package = "FLOWMAPR"), "shinyGUI"))
  print("in launch-GUI.R")
  print("files")
  print(files)
  print("var.remove")
  print(var.remove)
  print("var.annotate")
  print(var.annotate)
  print("clustering.var")
  print(clustering.var)
  print("cluster.numbers")
  print(cluster.numbers)
  print("subsamples")
  print(subsamples)
  print("distance.metric")
  print(distance.metric)
  print("minimum")
  print(minimum)
  print("maximum")
  print(maximum)
  print("per")
  print(per)
  print("save.folder")
  print(save.folder)
  print("mode")
  print(mode)
  FLOWMAPR::FLOWMAP(files = files, var.remove = var.remove, var.annotate = var.annotate,
                    clustering.var = clustering.var, cluster.numbers = cluster.numbers,
                    subsamples = subsamples, distance.metric = distance.metric,
                    minimum = minimum, maximum = maximum, per = per,
                    save.folder = save.folder, mode = mode,
                    shuffle = TRUE, name.sort = FALSE, downsample = FALSE)
  print("ending app")
  stopApp()
}


opendir <- function(dir = getwd()){
  if (.Platform['OS.type'] == "windows"){
    shell.exec(dir)
  } else {
    system(paste(Sys.getenv("R_BROWSER"), dir))
  }
}
