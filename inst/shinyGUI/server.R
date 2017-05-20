require(flowCore)
require(shiny)
require(rhandsontable)

shinyServer(function(input, output, session) {
  print("globe.inputs")
  print(globe.inputs)
  if (globe.inputs[["quit"]]) {
    stopApp()
  }
  options(shiny.maxRequestSize = 1000 * 1024^2)
  panel.info <- data.frame(channels = c(NA), removal = c(NA), cluster = c(NA), annotate = c(NA))
  operating.system <- Sys.info()[1]
  print("operating.system")
  print(operating.system)
  final.new.same <- NULL
  final.new.diff <- NULL
  dir.now <- globe.raw.FCS.dir
  FileOrder <- function(dir.now) {
    file.names <- list.files(dir.now, pattern = "\\.fcs")
    len.filenames <- seq(1, length(file.names))
    return(list(len.filenames = len.filenames,
                file.names = file.names))
  }
  file.info <- FileOrder(dir.now)
  len.filenames <- file.info$len.filenames
  file.names <- file.info$file.names
  observe({
    updateSelectInput(session, "checkGroup_files",
                      choices = paste(len.filenames, file.names, sep = " "))
  })
  testprint <- eventReactive(input$default.button, {
    print("BEEP")
    print("BOOP")
  })
  chosen_order <- eventReactive(input$gener.param.button, {
    input$file.order.input
  })
  fcs_order <- eventReactive(input$gener.param.button, {
    order <- as.numeric(unlist(strsplit(chosen_order(), ",")))
    fcs.list <- c()
    for(i in order) {
      fcs.list <- c(fcs.list, file.names[i])
      fcs.list
    }
  })
  contentdiff <- eventReactive(input$gener.param.button, {
    # Read input Files
    # Set the names
    fcs.list <- list()
    rows <- length(FileOrder(dir.now))
    count <- 0
    order <- as.numeric(unlist(strsplit(chosen_order(), ",")))
    for(i in order) {
      setwd(dir.now)
      files <- read.FCS(file.names[i], emptyValue = FALSE)
      filed <- pData(parameters(files))[, c("name")]
      name.desc <- do.call(paste, as.data.frame(filed, stringsAsFactors = FALSE))
      fcs.list[[(count + 1)]] <- name.desc
      count <- count + 1
      # Reads FCS Files, gets name and Description, add to a list of different FCS files
    }
    if (rows > 1) {
      same <- Reduce(intersect, fcs.list)
      every <- Reduce(union, fcs.list)
      diffs <- every
      diffs <- diffs[! every %in% same]
    } else {
      diffs <- NULL
    }
    # Gets different parameters from the FCS files
    final.new.diff <<- diffs
    diffs
    # If there is 1 FCS file, then there is no difference
  })
  contentsame <- eventReactive(input$gener.param.button, {
    fcs.list <- list()
    rows <- length(FileOrder(dir.now))
    count <- 0
    order <- as.numeric(unlist(strsplit(chosen_order(), ",")))
    for(i in order)
    {
      setwd(dir.now)
      files <- read.FCS(file.names[i], emptyValue = FALSE)
      filed <- pData(parameters(files))[, c("name")]
      name.desc <- do.call(paste, as.data.frame(filed, stringsAsFactors = FALSE))
      fcs.list[[(count + 1)]] <- name.desc
      count <- count + 1
      # Does same thing as above
    }
    same <- Reduce(intersect, fcs.list)
    every <- Reduce(union, fcs.list)
    diff <- every
    diff <- diff[! every %in% same]
    final.new.same <<- same
    same
    # gives the same paramters
  })
  tablecreate <- eventReactive(input$gener.param.button, {
    if (length(final.new.diff) == 0) {
      panel.info <<- data.frame(channels = c(final.new.same, final.new.diff),
                        removal = logical(length = length(final.new.same)),
                        cluster = logical(length = length(final.new.diff) + length(final.new.same)),
                        annotate = c(final.new.same, final.new.diff), stringsAsFactors = FALSE)
    } else {
      panel.info <<- data.frame(channels = c(final.new.same, final.new.diff),
                        removal = c(logical(length = length(final.new.same)),
                                    !logical(length = length(final.new.diff))),
                        cluster = logical(length = length(final.new.diff) + length(final.new.same)),
                        annotate = c(final.new.same, final.new.diff), stringsAsFactors = FALSE)
    }
    output$table <- renderRHandsontable({
      rhandsontable(panel.info) %>%
        hot_col("channels", readOnly = TRUE)
    })
    panel.info.edit <<- panel.info
  })
  observeEvent(input$gener.param.button, {
    if (length(final.new.diff) == 0) {
      panel.info <- data.frame(channels = c(final.new.same, final.new.diff),
                       removal = logical(length = length(final.new.same)),
                       cluster = logical(length = length(final.new.diff) + length(final.new.same)),
                       annotate = c(final.new.same, final.new.diff))
    } else {
      panel.info <- data.frame(channels = c(final.new.same, final.new.diff),
                       removal = c(logical(length = length(final.new.same)),
                                   !logical(length = length(final.new.diff))),
                       cluster = logical(length = length(final.new.diff) + length(final.new.same)),
                       annotate = c(final.new.same, final.new.diff))
    }
  })
  observe({
    updateSelectInput(session, "checkGroup_sim", choices = contentsame())
  })
  # updates the checkbox group to show same
  observe({
    updateSelectInput(session, "checkGroup_diff", choices = contentdiff())
  })
  file_merge_diff <- eventReactive(input$mbutton, {
    files_tbm <- input$checkGroup_diff
    merge_name <- input$filemerge
    new_diff <- contentdiff() [! contentdiff() %in% files_tbm]
    print(new_diff)
    print("DIFF WORKS")
    new_diff
  })
  file_merge_same <- eventReactive(input$mbutton, {
    new_same <- c(input$filemerge, contentsame())
    print(new_same)
    print("SAME WORKS")
    new_same
  })
  file_merge_table <- eventReactive(input$mbutton, {
    print("TABLE STARTED")
    files_tbm <- input$checkGroup_diff
    new.panel.info <- panel.info.edit
    print(input$checkGroup_diff)
    print("got here")
    for (i in files_tbm) {
      print("round")
      new.panel.info[new.panel.info$channels == i, "annotate"] <- input$filemerge
      panel.info.edit <<- new.panel.info
    }
    print("NEW DF MADE")
    output$table <- renderRHandsontable({
      rhandsontable(new.panel.info) %>%
        hot_col("channels", readOnly = TRUE)
    })
    print("TABLE SUCCESS")
  })
  observe({
    updateSelectInput(session, "checkGroup_sim", choices = file_merge_same())
  })
  # updates the checkbox group to show same
  observe({
    updateSelectInput(session, "checkGroup_diff", choices = file_merge_diff())
  })
  observe({
    file_merge_table()
  })
  write_file <- eventReactive(input$button, {
    if (operating.system != "Windows") {
      dir <- globe.result.dir
    } else {
      dir <- globe.result.dir
    }
    setwd(dir)
    # writes the file
    file.order <- as.numeric(unlist(strsplit(input$file.order.input, split = ",")))
    flowfile <- (hot_to_r(input$table))
    # print(folder.now)
    setwd(dir.now)
    set.seed(globe.inputs[["seedNum"]])
    files <- list.files(globe.raw.FCS.dir, full.names = TRUE, pattern = "\\.fcs")[file.order]
    # NEED MULTIFLOW-MAP FIX FOR FILES
    mode <- globe.inputs[["mode"]]
    save.folder <- globe.result.dir
    var.annotate <- list()
    for (j in 1:nrow(flowfile)) {
      var.annotate[[flowfile[j, 1]]] <- flowfile[j, 4]
    }
    var.remove <- flowfile[flowfile$removal == TRUE, 1]
    clustering.var <- flowfile[flowfile$cluster == TRUE, 1]
    per <- as.numeric(globe.inputs[["edgepctNum"]])
    maximum <- as.numeric(globe.inputs[["edgeMaxNum"]])
    minimum <- as.numeric(globe.inputs[["edgeminNum"]])
    distance.metric <- globe.inputs[["distanceMetric"]]
    subsamples <- as.numeric(globe.inputs[["subsampleNum"]])
    cluster.numbers <- as.numeric(globe.inputs[["clusterNum"]])
    seed.X <- as.numeric(globe.inputs[["seedNum"]])
    savePDFs <- as.logical(as.numeric(globe.inputs[["savePDFsToggle"]]))
    which.palette <- globe.inputs[["colorpalette"]]
    for (i in 1:length(clustering.var)) {
      clustering.var[i] <- var.annotate[[clustering.var[i]]]
    }
    name.sort <- FALSE
    downsample <- as.logical(as.numeric(globe.inputs[["downsampleToggle"]]))
    # Run FLOW-MAP
    if (downsample) {
      print("Downsampling")
      target.number <- subsamples
      subsamples <- FALSE
      target.percent <- NULL
      exclude.pctile <- 0.01
      target.pctile <- 0.99
      FLOWMAP(mode = mode, files = files, var.remove = var.remove, var.annotate = var.annotate,
              clustering.var = clustering.var, cluster.numbers = cluster.numbers,
              distance.metric = distance.metric, minimum = minimum, maximum = maximum,
              per = per, save.folder = save.folder, subsamples = subsamples,
              name.sort = name.sort, downsample = downsample, seed.X = seed.X,
              savePDFs = savePDFs, which.palette = which.palette,
              exclude.pctile = exclude.pctile, target.pctile = target.pctile,
              target.number = target.number, target.percent = target.percent)
    } else {
      print("No downsampling")
      FLOWMAP(seed.X = seed.X, files = files, var.remove = var.remove, var.annotate = var.annotate,
              clustering.var = clustering.var, cluster.numbers = cluster.numbers,
              subsamples = subsamples, distance.metric = distance.metric,
              minimum = minimum, maximum = maximum, per = per,
              save.folder = save.folder, mode = mode,
              name.sort = name.sort, downsample = downsample,
              savePDFs = savePDFs, which.palette = which.palette)
    }
    stopApp()
  })
  output$TESTPRINT <- renderText({
    testprint()
  })
  output$writefile <- renderText({
    write_file()
    NULL
  })
  output$vartable <- renderText({
    tablecreate()
    NULL
  })
  output$ordering <- renderText({
    chosen_order()
    NULL
  })
  output$fcsorder <- renderText({
    fcs_order()
  })
  output$testprint <- renderText({
    testprint()
    test.print <- "TEST PRINT"
    test.print
    "TEST PRINT"
  })
})
