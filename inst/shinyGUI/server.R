require(flowCore)
require(shiny)
require(rhandsontable)

shinyServer(function(input, output, session) {
  options(shiny.maxRequestSize = 1000 * 1024^2)
  DF <- data.frame(channels = c(NA), removal = c(NA), cluster = c(NA), annotate = c(NA))
  operating_system <<- Sys.info()[1]
  print("operating_system")
  print(operating_system)
  
  # get function for FLOW-MAP
  if (operating_system == "Windows"){
    folder_now <<- paste(gsub("/", "\\\\", getwd()), "\\", sep = "")
  } else {
    folder_now <<- paste(getwd(), "/", sep = "")
  }
  # get directory for where all FLOW-MAP function files are located
  final_new_same <<- NULL
  final_new_diff <<- NULL
  # Set Global Variables
  dir_now <- globe_resdir
  fileorder <- function(dir_now) {
    file_names <- list.files(dir_now, pattern = "\\.fcs")
    name_vec <- c()
    for (i in 1:length(file_names)){
      name_vec <- c(name_vec, i)
    }
    len_filenames <- name_vec
    return(list(len_filenames = len_filenames,
                file_names = file_names))
  }
  file_info <- fileorder(dir_now)
  len_filenames <- file_info$len_filenames
  file_names <- file_info$file_names
  observe({
    updateSelectInput(session, "checkGroup_files", choices = paste(len_filenames, file_names, sep = " "))
  })
  chosen_order <- eventReactive(input$generbutton2, {
    input$fileorder
  })
  fcs_order <- eventReactive(input$generbutton2, {
    order <- as.numeric(unlist(strsplit(chosen_order(), ",")))
    fcs_list <- c()
    for(i in order) {
      fcs_list <- c(fcs_list, file_names[i])
      fcs_list
    }
  })
  contentdiff <- eventReactive(input$generbutton2, {
    # Read input Files
    # Set the names
    fcs_list <- list()
    rows <- length(fileorder(dir_now))
    count <- 0
    order <- as.numeric(unlist(strsplit(chosen_order(), ",")))
    for(i in order) {
      setwd(dir_now)
      files <- read.FCS(file_names[i], emptyValue = FALSE)
      filed <- pData(parameters(files))[, c("name")]
      name_desc <- do.call(paste, as.data.frame(filed, stringsAsFactors = FALSE))
      fcs_list[[(count + 1)]] <- name_desc
      count <- count + 1
      # Reads FCS Files, gets name and Description, add to a list of different FCS files
    }
    if (rows > 1) {
      same <- Reduce(intersect, fcs_list)
      every <- Reduce(union, fcs_list)
      diffs <- every
      diffs <- diffs[! every %in% same]
    } else {
      diffs <- NULL
    }
    # Gets different parameters from the FCS files
    final_new_diff <<- diffs
    diffs
    # If there is 1 FCS file, then there is no difference
  })
  contentsame <- eventReactive(input$generbutton2, {
    fcs_list <- list()
    rows <- length(fileorder(dir_now))
    count <- 0
    order <- as.numeric(unlist(strsplit(chosen_order(), ",")))
    for(i in order)
    {
      setwd(dir_now)
      files <- read.FCS(file_names[i], emptyValue = FALSE)
      filed <- pData(parameters(files))[, c("name")]
      name_desc <- do.call(paste, as.data.frame(filed, stringsAsFactors = FALSE))
      fcs_list[[(count + 1)]] <- name_desc
      count <- count + 1
      # Does sdame thing as above
    }
    same <- Reduce(intersect, fcs_list)
    every <- Reduce(union, fcs_list)
    diff <- every
    diff <- diff[! every %in% same]
    final_new_same <<- same
    same
    # gives the same paramters
  })
  tablecreate <- eventReactive(input$generbutton2, {
    if (length(final_new_diff) == 0) {
      DF <<- data.frame(channels = c(final_new_same, final_new_diff),
                        removal = logical(length = length(final_new_same)),
                        cluster = logical(length = length(final_new_diff) + length(final_new_same)),
                        annotate = c(final_new_same, final_new_diff), stringsAsFactors = FALSE)
    } else {
      DF <<- data.frame(channels = c(final_new_same, final_new_diff),
                        removal = c(logical(length = length(final_new_same)),
                                    !logical(length = length(final_new_diff))),
                        cluster = logical(length = length(final_new_diff) + length(final_new_same)),
                        annotate = c(final_new_same, final_new_diff), stringsAsFactors = FALSE)
    }
    output$table <- renderRHandsontable({
      rhandsontable(DF) %>%
        hot_col("channels", readOnly = TRUE)
    })
    print(DF)
    DF_edit <<- DF
  })
  observeEvent(input$generbutton2, {
    if (length(final_new_diff) == 0) {
      DF <- data.frame(channels = c(final_new_same, final_new_diff),
                       removal = logical(length = length(final_new_same)),
                       cluster = logical(length = length(final_new_diff) + length(final_new_same)),
                       annotate = c(final_new_same, final_new_diff))
    } else {
      DF <- data.frame(channels = c(final_new_same, final_new_diff),
                       removal = c(logical(length = length(final_new_same)),
                                   !logical(length = length(final_new_diff))),
                       cluster = logical(length = length(final_new_diff) + length(final_new_same)),
                       annotate = c(final_new_same, final_new_diff))
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
    new_df <- DF_edit
    print(input$checkGroup_diff)
    print("got here")
    for (i in files_tbm) {
      print("round")
      new_df[new_df$channels == i, "annotate"] <- input$filemerge
      DF_edit <<- new_df
    }
    print("NEW DF MADE")
    output$table <- renderRHandsontable({
      rhandsontable(new_df) %>%
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
    if (operating_system != "Windows") {
      dir <- globe_resdir
    } else {
      dir <- globe_resdir
    }
    setwd(dir)
    # writes the file
    flowfile <- (hot_to_r(input$table))
    print(folder_now)
    setwd(dir_now)
    set.seed(globe_input[["seedNum"]])
    files <- c()
    files <- globe_resdir
    mode <- globe_input[["multiSingle"]]
    save.folder <- globe_resdir2
    var.annotate <- list()
    for (j in 1:nrow(flowfile)) {
      var.annotate[[flowfile[j, 1]]] <- flowfile[j, 4]
    }
    var.remove <- flowfile[flowfile$removal == T, 1]
    clustering.var <- flowfile[flowfile$cluster == T, 1]
    per <- as.numeric(globe_input[["edgepctNum"]])
    maximum <- as.numeric(globe_input[["edgeMaxNum"]])
    minimum <- as.numeric(globe_input[["edgeminNum"]])
    distance.metric <- globe_input[["distanceMetric"]]
    subsamples <- as.numeric(globe_input[["subsampleNum"]])
    cluster.numbers <- as.numeric(globe_input[["clusterNum"]])
    seed.X <- as.numeric(globe_input[["seedNum"]])
    for (i in 1:length(clustering.var)) {
      clustering.var[i] <- var.annotate[[clustering.var[i]]]
    }
    name.sort <- FALSE
    downsample <- FALSE
    savePDFs <- TRUE
    which.palette <- "bluered"
    # Run FLOW-MAP
    FLOWMAPR::FLOWMAP(seed.X = seed.X, files = files, var.remove = var.remove, var.annotate = var.annotate,
                      clustering.var = clustering.var, cluster.numbers = cluster.numbers,
                      subsamples = subsamples, distance.metric = distance.metric,
                      minimum = minimum, maximum = maximum, per = per,
                      save.folder = save.folder, mode = mode,
                      name.sort = name.sort, downsample = downsample,
                      savePDFs = savePDFs, which.palette = which.palette)
    stopApp()
  })
  output$stuff <- renderText({
    write_file()
    NULL
  })
  output$stuff3 <- renderText({
    tablecreate()
    NULL
  })
  output$stuff2 <- renderText({
    chosen_order()
    NULL
  })
  output$stuff4 <- renderText({
    fcs_order()
  })
})