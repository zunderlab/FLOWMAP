require(flowCore)
require(shiny)
require(rhandsontable)

print("globe.inputs")
print(globe.inputs)
print("globe.raw.FCS.dir")
print(globe.raw.FCS.dir)
print("globe.result.dir")
print(globe.result.dir)

FileOrder <- function(dir.now) {
  file.names <- list.files(dir.now, pattern = "\\.fcs")
  len.filenames <- seq(1, length(file.names))
  return(list(len.filenames = len.filenames,
              file.names = file.names))
}

InitializePanel <- function() {
  panel.info <- data.frame(channels = c(NA), removal = c(NA), cluster = c(NA), annotate = c(NA))
  return(panel.info)
}

GetMarkerNameParam <- function(file.iter, order, folder.name) {
  fcs.list <- list()
  temp.list <- list()
  file.iter <- file.iter[!is.na(file.iter)]
  # Reads FCS Files, gets name and Description, add to a list of different FCS files
  setwd(folder.name)
  for (i in 1:length(file.iter)) {
    fcs.file <- read.FCS(file.iter[order[i]], emptyValue = FALSE)
    fcs.name <- as.vector(fcs.file@parameters@data[, 1])
    fcs.param <- as.vector(fcs.file@parameters@data[, 2])
    temp.list[[1]] <- unlist(fcs.name)
    temp.list[[2]] <- unlist(fcs.param)
    final <- paste(temp.list[[1]], temp.list[[2]], sep = "_")
    fcs.list[[i]] <- final
  }
  return(list(fcs.list = fcs.list,
              temp.list = temp.list))
}

BuildVarAnnotate <- function(fcs.file, flowfile) {
  var.annotate <- list()
  original.names <- read.FCS(fcs.file)
  original.names <- unname(original.names@parameters@data[, 1])
  flowfile$annotate <- as.character(flowfile$annotate)
  for (i in 1:nrow(flowfile)) {
    var.annotate[[original.names[i]]] <- flowfile$annotate[i]
  }
  return(var.annotate)
}

SelectClusteringVar <- function(flowfile, var.annotate) {
  clustering.var.ind <- which(flowfile$cluster == TRUE)
  all.var <- unname(unlist(var.annotate))
  clustering.var <- all.var[clustering.var.ind]
  return(clustering.var)
}

SelectVarRemove <- function(flowfile, var.annotate) {
  var.remove.ind <- which(flowfile$removal == TRUE)
  all.var <- unname(unlist(var.annotate))
  var.remove <- all.var[var.remove.ind]
  return(var.remove)
}

ComparePanels <- function(fcs.list) {
  same <- Reduce(intersect, fcs.list)
  every <- Reduce(union, fcs.list)
  diffs <- every[! every %in% same]
  return(list(same = same,
              diffs = diffs))
}

UpdatePanel <- function(final.new.same, final.new.diff) {
  if (length(final.new.diff) == 0) {
    panel.info <- data.frame(channels = c(final.new.same, final.new.diff),
                             removal = logical(length = length(final.new.same)),
                             cluster = logical(length = length(final.new.diff) + length(final.new.same)),
                             annotate = c(final.new.same, final.new.diff), stringsAsFactors = FALSE)
  } else {
    panel.info <- data.frame(channels = c(final.new.same, final.new.diff),
                             removal = c(logical(length = length(final.new.same)),
                                         !logical(length = length(final.new.diff))),
                             cluster = logical(length = length(final.new.diff) + length(final.new.same)),
                             annotate = c(final.new.same, final.new.diff), stringsAsFactors = FALSE)
  }
  return(panel.info)
}

MakePanelOneMode <- function(final.new.same) {
  panel.info <- data.frame(channels = c(final.new.same),
                           removal = logical(length = length(final.new.same)),
                           cluster = logical(length = length(final.new.same)),
                           annotate = c(final.new.same), stringsAsFactors = FALSE)
  return(panel.info)
}

GetMultiFilePaths <- function(multi.list.global) {
  fcs.file.path <- c()
  for (i in multi.list.global) {
    fcs.file.path <- c(fcs.file.path, i)
  }
  fcs.file.path <- unlist(fcs.file.path) 
  if (paste(levels(fcs.file.path), collapse = "") != "") {
    fcs.file.path <- levels(droplevels(fcs.file.path))
  }
  return(fcs.file.path)
}

GetFilePathsfromCSV <- function(csv.path) {
  print(csv.path)
  csv.data <- read.csv(csv.path, header = TRUE)
  temp.csv <- csv.data[, 2:ncol(csv.data)]
  multi.list <- list()
  for (i in 1:nrow(temp.csv)) {
    temp.vec <- as.vector(temp.csv[i, ])
    temp.vec <- temp.vec[temp.vec != ""] 
    multi.list[[i]] <- temp.vec
  }
  names(multi.list) <- csv.data[, 1]
  return(multi.list)
}


# build server based on FLOW-MAP mode
if (globe.inputs[["mode"]] == "single") {
  server <- function(input, output, session) {
    if (globe.inputs[["quit"]]) {
      print("Exiting FLOWMAP")
      stopApp()
    }
    options(shiny.maxRequestSize = 1000 * 1024^2)
    panel.info <- InitializePanel()
    final.new.same <- NULL
    final.new.diff <- NULL
    file.info <- FileOrder(globe.raw.FCS.dir)
    len.filenames <- file.info$len.filenames
    file.names <- file.info$file.names
    if (identical(file.names, character(0))) {
      choice <- "No FCS Files"
    } else {
      choice <- paste(len.filenames, file.names, sep = " ")
    }
    observe({
      updateSelectInput(session, "check.group.files",
                        choices = choice)
    })
    ChosenOrder <- eventReactive(input$gener.param.button, {
      print(paste(len.filenames, sep = "", collapse = ", "))
      actual.input <- input$file.order.input
      if ( actual.input == "") {
        actual.input <- paste(len.filenames, sep = "", collapse = ",")
      }
      actual.input
    })
    GetFCSinOrder <- eventReactive(input$gener.param.button, {
      order <- as.numeric(unlist(strsplit(ChosenOrder(), ",")))
      fcs.list <- file.names[order]
      fcs.list
    })
    ContentDiff <- eventReactive(input$gener.param.button, {
      order <- as.numeric(unlist(strsplit(ChosenOrder(), ",")))
      print(order)
      temp.result <- GetMarkerNameParam(file.iter = file.names, order = order, folder.name = globe.raw.FCS.dir)
      fcs.list <- temp.result$fcs.list
      temp.list <- temp.result$temp.list
      diffs <- ComparePanels(fcs.list)$diffs
      # Gets different parameters from the FCS files
      final.new.diff <<- diffs
      diffs
    })
    ContentSame <- eventReactive(input$gener.param.button, {
      rows <- length(FileOrder(globe.raw.FCS.dir))
      order <- as.numeric(unlist(strsplit(ChosenOrder(), ",")))
      temp.result <- GetMarkerNameParam(file.iter = file.names, order = order, folder.name = globe.raw.FCS.dir)
      fcs.list <- temp.result$fcs.list
      temp.list <- temp.result$temp.list
      same <- ComparePanels(fcs.list)$same
      final.new.same <<- same
      same
      # gives the same paramters
    })
    TableCreate <- eventReactive(input$gener.param.button, {
      panel.info <<- UpdatePanel(final.new.same, final.new.diff)
      output$table <- renderRHandsontable({
        rhandsontable(panel.info) %>%
          hot_col("channels", readOnly = TRUE)
      })
      panel.info.edit <<- panel.info
    })
    observeEvent(input$gener.param.button, {
      panel.info <<- UpdatePanel(final.new.same, final.new.diff)
    })
    observe({
      updateSelectInput(session, "check.group.sim", choices = ContentSame())
    })
    # updates the checkbox group to show same
    observe({
      updateSelectInput(session, "check.group.diff", choices = ContentDiff())
    })
    FileMergeDiff <- eventReactive(input$merge.button, {
      files.tbm <- input$check.group.diff
      merge.name <- input$file.merge
      new.diff <- ContentDiff() [! ContentDiff() %in% files.tbm]
      print(new.diff)
      new.diff
    })
    FileMergeSame <- eventReactive(input$merge.button, {
      new.same <- c(input$file.merge, ContentSame())
      print(new.same)
      new.same
    })
    FileMergeTable <- eventReactive(input$merge.button, {
      files.tbm <- input$check.group.diff
      new.panel.info <- panel.info.edit
      print(input$check.group.diff)
      for (i in files.tbm) {
        new.panel.info[new.panel.info$channels == i, "annotate"] <- input$file.merge
        panel.info.edit <<- new.panel.info
      }
      print("new.panel.info")
      print(new.panel.info)
      output$table <- renderRHandsontable({
        rhandsontable(new.panel.info) %>%
          hot_col("channels", readOnly = TRUE)
      })
    })
    observe({
      updateSelectInput(session, "check.group.sim", choices = FileMergeSame())
    })
    # updates the checkbox group to show same
    observe({
      updateSelectInput(session, "check.group.diff", choices = FileMergeDiff())
    })
    observe({
      FileMergeTable()
    })
    WriteFile <- eventReactive(input$start.button, {
      file.order <- as.numeric(unlist(strsplit(ChosenOrder(), split = ",")))
      flowfile <- (hot_to_r(input$table))
      files <- list.files(globe.raw.FCS.dir, full.names = TRUE, pattern = "\\.fcs")[file.order]
      mode <- globe.inputs[["mode"]]
      save.folder <- globe.result.dir
      var.annotate <- BuildVarAnnotate(files[1], flowfile)
      var.remove <- SelectVarRemove(flowfile, var.annotate)
      clustering.var <- SelectClusteringVar(flowfile, var.annotate)
      maximum <- as.numeric(globe.inputs[["edge.max.num"]])
      minimum <- as.numeric(globe.inputs[["edge.min.num"]])
      distance.metric <- globe.inputs[["distance.metric"]]
      subsamples <- as.numeric(globe.inputs[["subsample.num"]])
      cluster.numbers <- as.numeric(globe.inputs[["cluster.num"]])
      seed.X <- as.numeric(globe.inputs[["seed.num"]])
      savePDFs <- as.logical(as.numeric(globe.inputs[["savePDFs.toggle"]]))
      which.palette <- globe.inputs[["color.palette"]]
      name.sort <- FALSE
      downsample <- as.logical(as.numeric(globe.inputs[["downsample.toggle"]]))
      
      print("var.annotate")
      print(var.annotate)
      print("clustering.var")
      print(clustering.var)
      
      # Run FLOW-MAP
      if (downsample) {
        print("Downsampling")
        target.number <- subsamples
        subsamples <- FALSE
        target.percent <- NULL
        exclude.pctile <- input$exclude.pctile
        target.pctile <- input$target.pctile
        FLOWMAP(mode = mode, files = files, var.remove = var.remove, var.annotate = var.annotate,
                clustering.var = clustering.var, cluster.numbers = cluster.numbers,
                distance.metric = distance.metric, minimum = minimum, maximum = maximum,
                save.folder = save.folder, subsamples = subsamples,
                name.sort = name.sort, downsample = downsample, seed.X = seed.X,
                savePDFs = savePDFs, which.palette = which.palette,
                exclude.pctile = exclude.pctile, target.pctile = target.pctile,
                target.number = target.number, target.percent = target.percent)
      } else {
        print("No Downsampling")
        FLOWMAP(mode = mode, files = files, var.remove = var.remove, var.annotate = var.annotate,
                clustering.var = clustering.var, cluster.numbers = cluster.numbers,
                distance.metric = distance.metric, minimum = minimum, maximum = maximum,
                save.folder = save.folder, subsamples = subsamples,
                name.sort = name.sort, downsample = downsample, seed.X = seed.X,
                savePDFs = savePDFs, which.palette = which.palette)
      }
      stopApp()
    })
    output$writefile <- renderText({
      WriteFile()
      NULL
    })
    output$vartable <- renderText({
      TableCreate()
      NULL
    })
    output$ordering <- renderText({
      ChosenOrder()
      NULL
    })
    output$fcsorder <- renderText({
      GetFCSinOrder()
    })
  }
} else if (globe.inputs[["mode"]] == "multi") {
  server <- function(input, output, session) {
    if (globe.inputs[["quit"]]) {
      stopApp()
    }
    options(shiny.maxRequestSize = 1000 * 1024^2)
    panel.info <- InitializePanel()
    final.new.same <- NULL
    final.new.diff <- NULL
    file.info <- FileOrder(globe.raw.FCS.dir)
    len.filenames <- file.info$len.filenames
    file.names <- file.info$file.names
    csv.order <- list.files(globe.raw.FCS.dir, pattern = "\\.csv")
    if(identical(csv.order, character(0))){
      choice <- "No CSV Files in Provided Folder!"
    } else {
      choice <- csv.order
    }
    observe({
      updateSelectInput(session, "check.group.csv",
                        choices = csv.order)
    })
    SelectCSV <- eventReactive(input$csv.finder, {
      input$check.group.csv
    })
    observeEvent(input$csv.finder, {
      print("Parsing CSV")
      csv.path <- paste(globe.raw.FCS.dir, SelectCSV(), sep = "/")
      multi.list <- GetFilePathsfromCSV(csv.path)
      multi.list.global <<- multi.list
      updateSelectInput(session, "check.group.files", 
                        choices = c("MultiFLOWMAP"))  
    })
    ContentDiff <- eventReactive(input$csv.finder, {
      fcs.file.path <- GetMultiFilePaths(multi.list.global)
      test.globe <<- fcs.file.path
      temp.result <- GetMarkerNameParam(file.iter = fcs.file.path,
                                        order = seq(1, length(fcs.file.path)),
                                        folder.name = globe.raw.FCS.dir)
      fcs.list <- temp.result$fcs.list
      temp.list <- temp.result$temp.list
      diffs <- ComparePanels(fcs.list)$diffs
      final.new.diff <<- diffs
      diffs
    })
    ContentSame <- eventReactive(input$csv.finder, {
      fcs.file.path <- GetMultiFilePaths(multi.list.global)
      all.files <<- fcs.file.path
      temp.result <- GetMarkerNameParam(file.iter = fcs.file.path,
                                        order = seq(1, length(fcs.file.path)),
                                        folder.name = globe.raw.FCS.dir)
      fcs.list <- temp.result$fcs.list
      temp.list <- temp.result$temp.list
      same <- ComparePanels(fcs.list)$same
      final.new.same <<- same
      same
      # gives the same parameters
    })
    TableCreate <- eventReactive(input$csv.finder, {
      panel.info <<- UpdatePanel(final.new.same, final.new.diff)
      output$table <- renderRHandsontable({
        rhandsontable(panel.info) %>%
          hot_col("channels", readOnly = TRUE)
      })
      panel.info.edit <<- panel.info
    })
    observeEvent(input$csv.finder, {
      panel.info <<- UpdatePanel(final.new.same, final.new.diff)
    })
    observe({
      updateSelectInput(session, "check.group.sim", choices = ContentSame())
    })
    # updates the checkbox group to show same
    observe({
      updateSelectInput(session, "check.group.diff", choices = ContentDiff())
    })
    FileMergeDiff <- eventReactive(input$merge.button, {
      files.tbm <- input$check.group.diff
      merge.name <- input$file.merge
      new.diff <- ContentDiff() [! ContentDiff() %in% files.tbm]
      print(new.diff)
      new.diff
    })
    FileMergeSame <- eventReactive(input$merge.button, {
      new.same <- c(input$file.merge, ContentSame())
      print(new.same)
      new.same
    })
    FileMergeTable <- eventReactive(input$merge.button, {
      files.tbm <- input$check.group.diff
      new.panel.info <- panel.info.edit
      print(input$check.group.diff)
      for (i in files.tbm) {
        new.panel.info[new.panel.info$channels == i, "annotate"] <- input$file.merge
        panel.info.edit <<- new.panel.info
      }
      print("new.panel.info")
      print(new.panel.info)
      output$table <- renderRHandsontable({
        rhandsontable(new.panel.info) %>%
          hot_col("channels", readOnly = TRUE)
      })
    })
    observe({
      updateSelectInput(session, "check.group.sim", choices = FileMergeSame())
    })
    # updates the checkbox group to show same
    observe({
      updateSelectInput(session, "check.group.diff", choices = FileMergeDiff())
    })
    observe({
      FileMergeTable()
    })
    WriteFile <- eventReactive(input$start.button, {
      flowfile <- (hot_to_r(input$table))
      files <- multi.list.global
      print(multi.list.global)
      mode <- globe.inputs[["mode"]]
      save.folder <- globe.result.dir
      var.annotate <- BuildVarAnnotate(files[[1]][1], flowfile)
      var.remove <- SelectVarRemove(flowfile, var.annotate)
      clustering.var <- SelectClusteringVar(flowfile, var.annotate)
      maximum <- as.numeric(globe.inputs[["edge.max.num"]])
      minimum <- as.numeric(globe.inputs[["edge.min.num"]])
      distance.metric <- globe.inputs[["distance.metric"]]
      subsamples <- as.numeric(globe.inputs[["subsample.num"]])
      cluster.numbers <- as.numeric(globe.inputs[["cluster.num"]])
      seed.X <- as.numeric(globe.inputs[["seed.num"]])
      savePDFs <- as.logical(as.numeric(globe.inputs[["savePDFs.toggle"]]))
      which.palette <- globe.inputs[["color.palette"]]
      name.sort <- FALSE
      downsample <- as.logical(as.numeric(globe.inputs[["downsample.toggle"]]))
      
      print("var.annotate")
      print(var.annotate)
      print("clustering.var")
      print(clustering.var)
      
      # Run FLOW-MAP
      if (downsample) {
        print("Downsampling")
        target.number <- subsamples
        subsamples <- FALSE
        target.percent <- NULL
        exclude.pctile <- input$exclude.pctile
        target.pctile <- input$target.pctile
        FLOWMAP(mode = mode, files = files, var.remove = var.remove, var.annotate = var.annotate,
                clustering.var = clustering.var, cluster.numbers = cluster.numbers,
                distance.metric = distance.metric, minimum = minimum, maximum = maximum,
                save.folder = save.folder, subsamples = subsamples,
                name.sort = name.sort, downsample = downsample, seed.X = seed.X,
                savePDFs = savePDFs, which.palette = which.palette,
                exclude.pctile = exclude.pctile, target.pctile = target.pctile,
                target.number = target.number, target.percent = target.percent)
      } else {
        print("No Downsampling")
        FLOWMAP(mode = mode, files = files, var.remove = var.remove, var.annotate = var.annotate,
                clustering.var = clustering.var, cluster.numbers = cluster.numbers,
                distance.metric = distance.metric, minimum = minimum, maximum = maximum,
                save.folder = save.folder, subsamples = subsamples,
                name.sort = name.sort, downsample = downsample, seed.X = seed.X,
                savePDFs = savePDFs, which.palette = which.palette)
      }
      stopApp()
    })
    output$writefile <- renderText({
      WriteFile()
      NULL
    })
    output$vartable <- renderText({
      TableCreate()
      NULL
    })
    output$ordering <- renderText({
      NULL
    })
    output$fcsorder <- renderText({
      NULL
    })
  }
} else if (globe.inputs[["mode"]] == "one") {
  server <- function(input, output, session) {
    if (globe.inputs[["quit"]]) {
      stopApp()
    }
    options(shiny.maxRequestSize = 1000 * 1024^2)
    panel.info <- InitializePanel()
    final.new.same <- NULL
    final.new.diff <- NULL
    file.info <- FileOrder(globe.raw.FCS.dir)
    len.filenames <- file.info$len.filenames
    file.names <- file.info$file.names
    if (identical(file.names, character(0))) {
      choice <<- "No FCS Files"
    } else {
      choice <<- file.names
    }
    observe({
      updateSelectInput(session, "check.group.files",
                        choices = choice)
    })
    ContentSame <- eventReactive(input$gener.param.button, {
      one.fcs <<- input$check.group.files
      temp.result <- GetMarkerNameParam(file.iter = one.fcs,
                                        order = c(1),
                                        folder.name = globe.raw.FCS.dir)
      fcs.list <- temp.result$fcs.list
      temp.list <- temp.result$temp.list
      same <- ComparePanels(fcs.list)$same
      final.new.same <<- same
      same
      # gives the same parameters
    })
    TableCreate <- eventReactive(input$gener.param.button, {
      panel.info <<- MakePanelOneMode(final.new.same)
      output$table <- renderRHandsontable({
        rhandsontable(panel.info) %>%
          hot_col("channels", readOnly = TRUE)
      })
    })
    observeEvent(input$gener.param.button, {
      panel.info <<- MakePanelOneMode(final.new.same)
    })
    observeEvent(input$gener.param.button, {
      updateSelectInput(session, "check.group.files", choices = ContentSame())
    })
    observeEvent(input$gener.param.button, {
      updateSelectInput(session, "check.group.files", choices = choice)
    })
    WriteFile <- eventReactive(input$start.button, {
      flowfile <- (hot_to_r(input$table))
      files <- one.fcs
      print(one.fcs)
      mode <- globe.inputs[["mode"]]
      save.folder <- globe.result.dir
      var.annotate <- BuildVarAnnotate(files, flowfile)
      var.remove <- SelectVarRemove(flowfile, var.annotate)
      clustering.var <- SelectClusteringVar(flowfile, var.annotate)
      maximum <- as.numeric(globe.inputs[["edge.max.num"]])
      minimum <- as.numeric(globe.inputs[["edge.min.num"]])
      distance.metric <- globe.inputs[["distance.metric"]]
      subsamples <- as.numeric(globe.inputs[["subsample.num"]])
      cluster.numbers <- as.numeric(globe.inputs[["cluster.num"]])
      seed.X <- as.numeric(globe.inputs[["seed.num"]])
      savePDFs <- as.logical(as.numeric(globe.inputs[["savePDFs.toggle"]]))
      which.palette <- globe.inputs[["color.palette"]]
      name.sort <- FALSE
      downsample <- as.logical(as.numeric(globe.inputs[["downsample.toggle"]]))
      files <- paste(globe.raw.FCS.dir, files, sep = "/")
      
      print("var.annotate")
      print(var.annotate)
      print("clustering.var")
      print(clustering.var)
      
      # Run FLOW-MAP
      if (downsample) {
        print("Downsampling")
        target.number <- subsamples
        subsamples <- FALSE
        target.percent <- NULL
        exclude.pctile <- input$exclude.pctile
        target.pctile <- input$target.pctile
        FLOWMAP(mode = mode, files = files, var.remove = var.remove, var.annotate = var.annotate,
                clustering.var = clustering.var, cluster.numbers = cluster.numbers,
                distance.metric = distance.metric, minimum = minimum, maximum = maximum,
                save.folder = save.folder, subsamples = subsamples,
                name.sort = name.sort, downsample = downsample, seed.X = seed.X,
                savePDFs = savePDFs, which.palette = which.palette,
                exclude.pctile = exclude.pctile, target.pctile = target.pctile,
                target.number = target.number, target.percent = target.percent)
      } else {
        print("No Downsampling")
        FLOWMAP(mode = mode, files = files, var.remove = var.remove, var.annotate = var.annotate,
                clustering.var = clustering.var, cluster.numbers = cluster.numbers,
                distance.metric = distance.metric, minimum = minimum, maximum = maximum,
                save.folder = save.folder, subsamples = subsamples,
                name.sort = name.sort, downsample = downsample, seed.X = seed.X,
                savePDFs = savePDFs, which.palette = which.palette)
      }
      stopApp()
    })
    output$writefile <- renderText({
      WriteFile()
      NULL
    })
    output$vartable <- renderText({
      TableCreate()
      NULL
    })
    output$ordering <- renderText({
      NULL
    })
    output$fcsorder <- renderText({
      NULL
    })
  }
}


# build UI based on FLOW-MAP mode
if (globe.inputs[["downsample.toggle"]] == "1") {
  if (globe.inputs[["mode"]] == "single") {
    ui <- fluidPage(
      titlePanel("File Uploader"),
      fluidRow(
        column(width = 3, selectInput("check.group.files",
                                      label = h5("Uploaded Order"),
                                      choices = "Pending Upload",
                                      selected = NULL,
                                      multiple = TRUE,
                                      selectize = FALSE,
                                      size = 7),
               textInput("file.order.input", label = h5("Write the FCS File Order"),
                         placeholder = "Ex: 4, 2, 7, 5, 3, 1, 6"),
               numericInput("target.pctile", label = h5("Downsample Target Percentile"), value = 0.99),
               numericInput("exclude.pctile", label = h5("Downsample Exclude Percentile"), value = 0.01),
               actionButton("gener.param.button", "Read Panel from FCS Files"),
               textOutput("writefile"),
               textOutput("vartable"),
               textOutput("ordering"),
               textOutput("fcsorder")),
        column(width = 3, selectInput("check.group.sim",
                                      label = h5("Matching Channels in FCS Files"),
                                      choices = "Pending Upload",
                                      selected = NULL,
                                      multiple = TRUE,
                                      selectize = FALSE,
                                      size = 7),
               selectInput("check.group.diff",
                           label = h5("Different Channels in FCS Files"),
                           choices = "Pending Upload",
                           selected = NULL,
                           multiple = TRUE,
                           selectize = FALSE,
                           size = 7),
               textInput("file.merge", label = h5("Select New Merged Channel Name"), placeholder = "New Name"),
               actionButton("merge.button", "Merge Selected Channels")),
        column(width = 5,
               actionButton("start.button", "Run FLOWMAPR"),
               br(),
               rHandsontableOutput("table", width = 600))
      )
    )
  } else if (globe.inputs[["mode"]] == "multi") {
    ui <- fluidPage(
      titlePanel("File Uploader"),
      fluidRow(
        column(width = 3, selectInput("check.group.csv",
                                      label = h5("Import CSV"),
                                      choices = "Pending Upload",
                                      selected = NULL,
                                      multiple = FALSE,
                                      selectize = FALSE,
                                      size = 3),
               numericInput("target.pctile", label = h5("Downsample Target Percentile"), value = 0.99),
               numericInput("exclude.pctile", label = h5("Downsample Exclude Percentile"), value = 0.01),
               actionButton("csv.finder", "Input CSV"),
               textOutput("writefile"),
               textOutput("vartable"),
               textOutput("ordering"),
               textOutput("fcsorder")),
        column(width = 3, selectInput("check.group.sim",
                                      label = h5("Matching Channels in FCS Files"),
                                      choices = "Pending Upload",
                                      selected = NULL,
                                      multiple = TRUE,
                                      selectize = FALSE,
                                      size = 7),
               selectInput("check.group.diff",
                           label = h5("Different Channels in FCS Files"),
                           choices = "Pending Upload",
                           selected = NULL,
                           multiple = TRUE,
                           selectize = FALSE,
                           size = 7),
               textInput("file.merge", label = h5("Select New Merged Channel Name"), placeholder = "New Name"),
               actionButton("merge.button", "Merge Selected Channels")),
        column(width = 5,
               actionButton("start.button", "Run FLOWMAPR"),
               br(),
               rHandsontableOutput("table", width = 600))
      )
    )
  } else if (globe.inputs[["mode"]] == "one") {
    ui <- fluidPage(
      titlePanel("File Uploader"),
      fluidRow(
        column(width = 3, selectInput("check.group.files",
                                      label = h5("Uploaded Order"),
                                      choices = "Pending Upload",
                                      selected = NULL,
                                      multiple = TRUE,
                                      selectize = FALSE,
                                      size = 7),
               numericInput("target.pctile", label = h5("Downsample Target Percentile"), value = 0.99),
               numericInput("exclude.pctile", label = h5("Downsample Exclude Percentile"), value = 0.01),
               actionButton("gener.param.button", "Read Panel from FCS File"),
               textOutput("writefile"),
               textOutput("vartable"),
               textOutput("ordering"),
               textOutput("fcsorder")),
        column(width = 5,
               actionButton("start.button", "Run FLOWMAPR"),
               br(),
               rHandsontableOutput("table", width = 600))
      )
    )
  }
} else if (globe.inputs[["downsample.toggle"]] == "0") {
  if (globe.inputs[["mode"]] == "single") {
    ui <- fluidPage(
      titlePanel("File Uploader"),
      fluidRow(
        column(width = 3, selectInput("check.group.files",
                                      label = h5("Uploaded Order"),
                                      choices = "Pending Upload",
                                      selected = NULL,
                                      multiple = TRUE,
                                      selectize = FALSE,
                                      size = 7),
               textInput("file.order.input", label = h5("Write the FCS File Order"),
                         placeholder = "Ex: 4, 2, 7, 5, 3, 1, 6"),
               actionButton("gener.param.button", "Read Panel from FCS Files"),
               textOutput("writefile"),
               textOutput("vartable"),
               textOutput("ordering"),
               textOutput("fcsorder")),
        column(width = 3, selectInput("check.group.sim",
                                      label = h5("Matching Channels in FCS Files"),
                                      choices = "Pending Upload",
                                      selected = NULL,
                                      multiple = TRUE,
                                      selectize = FALSE,
                                      size = 7),
               selectInput("check.group.diff",
                           label = h5("Different Channels in FCS Files"),
                           choices = "Pending Upload",
                           selected = NULL,
                           multiple = TRUE,
                           selectize = FALSE,
                           size = 7),
               textInput("file.merge", label = h5("Select New Merged Channel Name"), placeholder = "New Name"),
               actionButton("merge.button", "Merge Selected Channels")),
        column(width = 5,
               actionButton("start.button", "Run FLOWMAPR"),
               br(),
               br(),
               rHandsontableOutput("table", width = 600))
      )
    )
  } else if (globe.inputs[["mode"]] == "multi") {
    ui <- fluidPage(
      titlePanel("File Uploader"),
      fluidRow(
        column(width = 3, selectInput("check.group.csv",
                                      label = h5("Import CSV"),
                                      choices = "Pending Upload",
                                      selected = NULL,
                                      multiple = FALSE,
                                      selectize = FALSE,
                                      size = 3),
               actionButton("csv.finder", "Input CSV"),
               textOutput("writefile"),
               textOutput("vartable"),
               textOutput("ordering"),
               textOutput("fcsorder")),
        column(width = 3, selectInput("check.group.sim",
                                      label = h5("Matching Channels in FCS Files"),
                                      choices = "Pending Upload",
                                      selected = NULL,
                                      multiple = TRUE,
                                      selectize = FALSE,
                                      size = 7),
               selectInput("check.group.diff",
                           label = h5("Different Channels in FCS Files"),
                           choices = "Pending Upload",
                           selected = NULL,
                           multiple = TRUE,
                           selectize = FALSE,
                           size = 7),
               textInput("file.merge", label = h5("Select New Merged Channel Name"), placeholder = "New Name"),
               actionButton("merge.button", "Merge Selected Channels")),
        column(width = 5,
               actionButton("start.button", "Run FLOWMAPR"),
               br(),
               br(),
               rHandsontableOutput("table", width = 600))
      )
    )
  } else if (globe.inputs[["mode"]] == "one") {
    ui <- fluidPage(
      titlePanel("File Uploader"),
      fluidRow(
        column(width = 3, selectInput("check.group.files",
                                      label = h5("Uploaded Order"),
                                      choices = "Pending Upload",
                                      selected = NULL,
                                      multiple = TRUE,
                                      selectize = FALSE,
                                      size = 7),
               actionButton("gener.param.button", "Read Panel from FCS File"),
               textOutput("writefile"),
               textOutput("vartable"),
               textOutput("ordering"),
               textOutput("fcsorder")),
        column(width = 8,
               actionButton("start.button", "Run FLOWMAPR"),
               br(),
               br(),
               rHandsontableOutput("table", width = 600))
      )
    )
  }
} else {
  stop("Unrecognized downsampling selection!")
}

shinyApp(ui = ui, server = server)
