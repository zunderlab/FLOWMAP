require(flowCore)
require(shiny)
require(rhandsontable)
if(globe.inputs[["mode"]] == "single"){
  shinyServer(function(input, output, session) {
    print("globe.inputs")
    print(globe.inputs)
  
    print("globe.raw.FCS.dir")
    print(globe.raw.FCS.dir)
  
    print("globe.result.dir")
    print(globe.result.dir)
    
    if (globe.inputs[["quit"]]) {
      stopApp()
    }
    options(shiny.maxRequestSize = 1000 * 1024^2)
    panel.info <- data.frame(channels = c(NA), removal = c(NA), cluster = c(NA), annotate = c(NA))
    final.new.same <- NULL
    final.new.diff <- NULL
    FileOrder <- function(dir.now) {
      file.names <- list.files(dir.now, pattern = "\\.fcs")
      len.filenames <- seq(1, length(file.names))
      return(list(len.filenames = len.filenames,
                  file.names = file.names))
    }
    file.info <- FileOrder(globe.raw.FCS.dir)
    len.filenames <- file.info$len.filenames
    file.names <- file.info$file.names
    if(identical(file.names, character(0))){
      choice = "No FCS Files"
    } else {
      choice = paste(len.filenames, file.names, sep = " ")
    }
    observe({
      updateSelectInput(session, "check.group.files",
                        choices = choice)
    })
    ChosenOrder <- eventReactive(input$gener.param.button, {
        print(paste(len.filenames, sep = "", collapse = ", "))
        actual.input <- input$file.order.input
        if( actual.input == ""){
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
        fcs.list <- list()
        rows <- length(FileOrder(globe.raw.FCS.dir))
        count <- 0
        order <- as.numeric(unlist(strsplit(ChosenOrder(), ",")))
        print(order)
        file.names <- file.names[!is.na(file.names)]
        for(i in order) {
          setwd(globe.raw.FCS.dir)
          print(file.names[i])
          
          fcs.files <- read.FCS(file.names[i], emptyValue = FALSE)
          fcs.files.desc <- pData(parameters(fcs.files))[, c("name")]
          name.desc <- do.call(paste, as.data.frame(fcs.files.desc, stringsAsFactors = FALSE))
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
    ContentSame <- eventReactive(input$gener.param.button, {
        fcs.list <- list()
        rows <- length(FileOrder(globe.raw.FCS.dir))
        count <- 0
        order <- as.numeric(unlist(strsplit(ChosenOrder(), ",")))
        for(i in order) {
          setwd(globe.raw.FCS.dir)
          fcs.files <- read.FCS(file.names[i], emptyValue = FALSE)
          fcs.files.desc <- pData(parameters(fcs.files))[, c("name")]
          name.desc <- do.call(paste, as.data.frame(fcs.files.desc, stringsAsFactors = FALSE))
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
    TableCreate <- eventReactive(input$gener.param.button, {
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
      print("DIFF WORKS")
      new.diff
    })
    FileMergeSame <- eventReactive(input$merge.button, {
      new.same <- c(input$file.merge, ContentSame())
      print(new.same)
      print("SAME WORKS")
      new.same
    })
    FileMergeTable <- eventReactive(input$merge.button, {
      print("TABLE STARTED")
      files.tbm <- input$check.group.diff
      new.panel.info <- panel.info.edit
      print(input$check.group.diff)
      print("got here")
      for (i in files.tbm) {
        print("round")
        new.panel.info[new.panel.info$channels == i, "annotate"] <- input$file.merge
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
      # file.order <- as.numeric(unlist(strsplit(input$file.order.input, split = ",")))
      flowfile <- (hot_to_r(input$table))
      print("flowfile")
      print(flowfile)
      setwd(globe.raw.FCS.dir)
      set.seed(globe.inputs[["seed.num"]])
      files <- list.files(globe.raw.FCS.dir, full.names = TRUE, pattern = "\\.fcs")[file.order]
      # NEED MULTI-FLOWMAP FIX FOR FILES
      mode <- globe.inputs[["mode"]]
      save.folder <- globe.result.dir
      var.annotate <- list()
      for (j in 1:nrow(flowfile)) {
        var.annotate[[flowfile[j, 1]]] <- flowfile[j, 4]
      }
      var.remove <- flowfile[flowfile$removal == TRUE, 1]
      clustering.var <- flowfile[flowfile$cluster == TRUE, 1]
      per <- as.numeric(globe.inputs[["edge.pct.num"]])
      maximum <- as.numeric(globe.inputs[["edge.max.num"]])
      minimum <- as.numeric(globe.inputs[["edge.min.num"]])
      distance.metric <- globe.inputs[["distance.metric"]]
      subsamples <- as.numeric(globe.inputs[["subsample.num"]])
      cluster.numbers <- as.numeric(globe.inputs[["cluster.num"]])
      seed.X <- as.numeric(globe.inputs[["seed.num"]])
      savePDFs <- as.logical(as.numeric(globe.inputs[["savePDFs.toggle"]]))
      which.palette <- globe.inputs[["color.palette"]]
      for (i in 1:length(clustering.var)) {
        clustering.var[i] <- var.annotate[[clustering.var[i]]]
      }
      name.sort <- FALSE
      downsample <- as.logical(as.numeric(globe.inputs[["downsample.toggle"]]))
      
      print("output")
      print(output)
      
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
        print(files)
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
    # output$TESTPRINT <- renderText({
    #   TestPrint()
    # })
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
    output$testprint <- renderText({
      TestPrint()
      test.print <- "TEST PRINT"
      test.print
      "TEST PRINT"
    })
  })
} else if(globe.inputs[["mode"]] == "multi"){
  shinyServer(function(input, output, session) {
    print("globe.inputs")
    print(globe.inputs)
    
    print("globe.raw.FCS.dir")
    print(globe.raw.FCS.dir)
    
    print("globe.result.dir")
    print(globe.result.dir)
    
    if (globe.inputs[["quit"]]) {
      stopApp()
    }
    # multi.start <<- FALSE
    options(shiny.maxRequestSize = 1000 * 1024^2)
    panel.info <- data.frame(channels = c(NA), removal = c(NA), cluster = c(NA), annotate = c(NA))
    final.new.same <- NULL
    final.new.diff <- NULL
    FileOrder <- function(dir.now) {
      file.names <- list.files(dir.now, pattern = "\\.fcs")
      len.filenames <- seq(1, length(file.names))
      return(list(len.filenames = len.filenames,
                  file.names = file.names))
    }
    file.info <- FileOrder(globe.raw.FCS.dir)
    len.filenames <- file.info$len.filenames
    file.names <- file.info$file.names
    
    csv.order = list.files(globe.raw.FCS.dir, pattern = "\\.csv")
    if(identical(csv.order, character(0))){
      choice = "No CSV Files"
    } else {
      choice = csv.order
    }
    observe({
      updateSelectInput(session, "check.group.csv",
                        choices = csv.order)
    })
    SelectCsv <- eventReactive(input$csv.finder, {
      input$check.group.csv
    })
    observeEvent(input$csv.finder, {
      print("Parsing CSV")
      csv.path = paste(globe.raw.FCS.dir, SelectCsv(), sep = "/")
      print(csv.path)
      csv.data = read.csv(csv.path, header = FALSE)
      multi.list = list()
      temp.vec = c()
      for(i in 1:nrow(csv.data)){
        for(j in 1:length(csv.data[i,])){
          if(csv.data[i, ][j] != ""){
            temp.vec = c(temp.vec, csv.data[i, ][j])
          }
        }
        multi.list[[i]] = temp.vec
        temp.vec = c()
      }
      print(multi.list)
      multi.list.global <<- multi.list
      updateSelectInput(session, "check.group.files", 
                        choices = c("MultiFlowMap"))  
    })
    
    ContentDiff <- eventReactive(input$csv.finder, {
        count = 0
        fcs.list <- list()
        fcs.file.path <- c()
        for(i in multi.list.global){
          fcs.file.path <- c(fcs.file.path, i)
        }
        fcs.file.path <- unlist(fcs.file.path) 
        fcs.file.path <- levels(droplevels(fcs.file.path))
        print(fcs.file.path)
        test.globe <<- fcs.file.path
        for(i in fcs.file.path) {
          fcs.files <- read.FCS(i, emptyValue = FALSE)
          fcs.files.desc <- pData(parameters(fcs.files))[, c("name")]
          name.desc <- do.call(paste, as.data.frame(fcs.files.desc, stringsAsFactors = FALSE))
          fcs.list[[(count + 1)]] <- name.desc
          count <- count + 1
          # Reads FCS Files, gets name and Description, add to a list of different FCS files
        }
        print(fcs.list)
        same <- Reduce(intersect, fcs.list)
        every <- Reduce(union, fcs.list)
        diffs <- every
        diffs <- diffs[! every %in% same]
        final.new.diff <<- diffs
        diffs
      # }
    })
    ContentSame <- eventReactive(input$csv.finder, {
        fcs.list <- list()
        fcs.file.path <- c()
        for(i in multi.list.global){
          fcs.file.path <- c(fcs.file.path, i)
        }
        fcs.file.path <- unlist(fcs.file.path) 
        fcs.file.path <- levels(droplevels(fcs.file.path))
        all.files <<- fcs.file.path
        print(fcs.file.path)
        count = 0
        for(i in fcs.file.path) {
          print(i)
          fcs.files <- read.FCS(i, emptyValue = FALSE)
          fcs.files.desc <- pData(parameters(fcs.files))[, c("name")]
          name.desc <- do.call(paste, as.data.frame(fcs.files.desc, stringsAsFactors = FALSE))
          fcs.list[[(count + 1)]] <- name.desc
          count <- count + 1
          # Reads FCS Files, gets name and Description, add to a list of different FCS files
        }
        same <- Reduce(intersect, fcs.list)
        every <- Reduce(union, fcs.list)
        diff <- every
        diff <- diff[! every %in% same]
        final.new.same <<- same
        same
      # }
      # gives the same paramters
    })
    TableCreate <- eventReactive(input$csv.finder, {
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
    observeEvent(input$csv.finder, {
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
      print("DIFF WORKS")
      new.diff
    })
    FileMergeSame <- eventReactive(input$merge.button, {
      new.same <- c(input$file.merge, ContentSame())
      print(new.same)
      print("SAME WORKS")
      new.same
    })
    FileMergeTable <- eventReactive(input$merge.button, {
      print("TABLE STARTED")
      files.tbm <- input$check.group.diff
      new.panel.info <- panel.info.edit
      print(input$check.group.diff)
      print("got here")
      for (i in files.tbm) {
        print("round")
        new.panel.info[new.panel.info$channels == i, "annotate"] <- input$file.merge
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
      print("flowfile")
      print(flowfile)
      setwd(globe.raw.FCS.dir)
      set.seed(globe.inputs[["seed.num"]])
      files <- all.files
      # NEED MULTI-FLOWMAP FIX FOR FILES
      mode <- globe.inputs[["mode"]]
      save.folder <- globe.result.dir
      var.annotate <- list()
      for (j in 1:nrow(flowfile)) {
        var.annotate[[flowfile[j, 1]]] <- flowfile[j, 4]
      }
      var.remove <- flowfile[flowfile$removal == TRUE, 1]
      clustering.var <- flowfile[flowfile$cluster == TRUE, 1]
      per <- as.numeric(globe.inputs[["edge.pct.num"]])
      maximum <- as.numeric(globe.inputs[["edge.max.num"]])
      minimum <- as.numeric(globe.inputs[["edge.min.num"]])
      distance.metric <- globe.inputs[["distance.metric"]]
      subsamples <- as.numeric(globe.inputs[["subsample.num"]])
      cluster.numbers <- as.numeric(globe.inputs[["cluster.num"]])
      seed.X <- as.numeric(globe.inputs[["seed.num"]])
      savePDFs <- as.logical(as.numeric(globe.inputs[["savePDFs.toggle"]]))
      which.palette <- globe.inputs[["color.palette"]]
      for (i in 1:length(clustering.var)) {
        clustering.var[i] <- var.annotate[[clustering.var[i]]]
      }
      name.sort <- FALSE
      downsample <- as.logical(as.numeric(globe.inputs[["downsample.toggle"]]))
      
      print("output")
      print(output)
      
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
        print(files)
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
    # output$TESTPRINT <- renderText({
    #   TestPrint()
    # })
    output$writefile <- renderText({
      WriteFile()
      NULL
    })
    output$vartable <- renderText({
      TableCreate()
      NULL
    })
    output$ordering <- renderText({
      # ChosenOrder()
      NULL
    })
    output$fcsorder <- renderText({
      # GetFCSinOrder()
      NULL
    })
    output$testprint <- renderText({
      TestPrint()
      test.print <- "TEST PRINT"
      test.print
      "TEST PRINT"
    })
  })
} else if(globe.inputs[["mode"]] == "one"){
  shinyServer(function(input, output, session) {
    print("globe.inputs")
    print(globe.inputs)
    
    print("globe.raw.FCS.dir")
    print(globe.raw.FCS.dir)
    
    print("globe.result.dir")
    print(globe.result.dir)
    
    if (globe.inputs[["quit"]]) {
      stopApp()
    }
    options(shiny.maxRequestSize = 1000 * 1024^2)
    panel.info <- data.frame(channels = c(NA), removal = c(NA), cluster = c(NA), annotate = c(NA))
    final.new.same <- NULL
    final.new.diff <- NULL
    FileOrder <- function(dir.now) {
      file.names <- list.files(dir.now, pattern = "\\.fcs")
      len.filenames <- seq(1, length(file.names))
      return(list(len.filenames = len.filenames,
                  file.names = file.names))
    }
    file.info <- FileOrder(globe.raw.FCS.dir)
    len.filenames <- file.info$len.filenames
    file.names <- file.info$file.names
    if(identical(file.names, character(0))){
      choice <<- "No FCS Files"
    } else {
      choice <<- file.names
    }
    observe({
      updateSelectInput(session, "check.group.files",
                        choices = choice)
    })
    
    ContentSame <- eventReactive(input$gener.param.button, {
      fcs.list <- list()
        setwd(globe.raw.FCS.dir)
        fcs.files <- read.FCS(input$check.group.files, emptyValue = FALSE)
        fcs.files.desc <- pData(parameters(fcs.files))[, c("name")]
        name.desc <- do.call(paste, as.data.frame(fcs.files.desc, stringsAsFactors = FALSE))
        fcs.list[[(1)]] <- name.desc
        # Does same thing as above
      same <- Reduce(intersect, fcs.list)
      every <- Reduce(union, fcs.list)
      diff <- every
      diff <- diff[! every %in% same]
      final.new.same <<- same
      same
      # gives the same paramters
    })
    TableCreate <- eventReactive(input$gener.param.button, {
      panel.info <<- data.frame(channels = c(final.new.same),
                                removal = logical(length = length(final.new.same)),
                                cluster = logical(length = length(final.new.same)),
                                annotate = c(final.new.same), stringsAsFactors = FALSE)
      output$table <- renderRHandsontable({
        rhandsontable(panel.info) %>%
          hot_col("channels", readOnly = TRUE)
      })
      # panel.info.edit <<- panel.info
    })
    observeEvent(input$gener.param.button, {
        panel.info <<- data.frame(channels = c(final.new.same),
                                  removal = logical(length = length(final.new.same)),
                                  cluster = logical(length = length(final.new.same)),
                                  annotate = c(final.new.same), stringsAsFactors = FALSE)
    })
    
    observeEvent(input$gener.param.button, {
      updateSelectInput(session, "check.group.files", choices = ContentSame())
    })
    
    observeEvent(input$gener.param.button, {
      updateSelectInput(session, "check.group.files", choices = choice)
    })
    
    
    WriteFile <- eventReactive(input$start.button, {
      flowfile <- (hot_to_r(input$table))
      print("flowfile")
      print(flowfile)
      setwd(globe.raw.FCS.dir)
      set.seed(globe.inputs[["seed.num"]])
      files <- list.files(input$check.group.files)
      # NEED MULTI-FLOWMAP FIX FOR FILES
      mode <- globe.inputs[["mode"]]
      save.folder <- globe.result.dir
      var.annotate <- list()
      for (j in 1:nrow(flowfile)) {
        var.annotate[[flowfile[j, 1]]] <- flowfile[j, 4]
      }
      var.remove <- flowfile[flowfile$removal == TRUE, 1]
      clustering.var <- flowfile[flowfile$cluster == TRUE, 1]
      per <- as.numeric(globe.inputs[["edge.pct.num"]])
      maximum <- as.numeric(globe.inputs[["edge.max.num"]])
      minimum <- as.numeric(globe.inputs[["edge.min.num"]])
      distance.metric <- globe.inputs[["distance.metric"]]
      subsamples <- as.numeric(globe.inputs[["subsample.num"]])
      cluster.numbers <- as.numeric(globe.inputs[["cluster.num"]])
      seed.X <- as.numeric(globe.inputs[["seed.num"]])
      savePDFs <- as.logical(as.numeric(globe.inputs[["savePDFs.toggle"]]))
      which.palette <- globe.inputs[["color.palette"]]
      for (i in 1:length(clustering.var)) {
        clustering.var[i] <- var.annotate[[clustering.var[i]]]
      }
      name.sort <- FALSE
      downsample <- as.logical(as.numeric(globe.inputs[["downsample.toggle"]]))
      
      print("output")
      print(output)
      
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
        print(files)
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
    output$testprint <- renderText({
      TestPrint()
      test.print <- "TEST PRINT"
      test.print
      "TEST PRINT"
    })
  })
}