#SMG 5.10.18

library(shiny)
library(shinythemes)
library(shinyFiles)
library(shinydashboard)
library(shinyalert)
library(flowCore)
library(ShinyInputDir)
library(rhandsontable)

#Raw FCS Files or CSV Directory (for mode multiFLOW-MAP):
# help button
# files/directory select
#Results Directory:
# help button
# files/directory select
#SPADE Downsampling:
#help
#checkbox
#Save Graph PDFs:
#help
#checkbox
#Distance Metric
#select input
#manhattan, euclidean
#FLOW-MAP Mode
#selectInput
#bluered, jet, cb
#Subsample Number (Downsample Target Number)
#numericInput
#Number of Clusters
#Set Seed Number
#Minimun Number of Edges
#Maximum Number of Edges

#Reset, Quit, Submit buttons

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

#Initialize globe.inputs

globe.inputs <<- list()



###########################################################################

########################################################################


header <- dashboardHeader(title = "FLOWMAPR")

sidebar <- dashboardSidebar(
  sidebarMenu(
    #menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
    #menuItem("Widgets", tabName = "widgets", icon = icon("th")),
    menuItem("Settings", tabName = "params", icon = icon("cog", lib = "glyphicon")),
    menuItem("File Processing", tabName = "files", icon = icon("folder-open", lib = "glyphicon")),
    menuItem("Help", tabName = "help", icon = icon("question-sign", lib = "glyphicon")),
    actionButton("reset", "Reset"),
    actionButton("quit", "Quit")
  )#sideBarMenu
)#dashboardSidebar

body <- dashboardBody(
  tabItems(
    # First tab content
    tabItem(tabName = "params",
            fluidRow(
              uiOutput('resetable_input1')
            ),
            fluidRow(
              box( # Input: Select a file ----
                   title = "Directory Selection",
                   #tags$hr(),
                   tags$h5("Raw FCS Files or CSV Directory (for mode multiFLOW-MAP):"),
                   div(style="display: inline-block;vertical-align:top; width: 75px;",
                       shinyDirButton("dirIn", "Choose...", "Upload")),
                   div(style="display: inline-block;vertical-align:top; width: 10px;",
                       HTML("<br>")),
                   div(style="display: inline-block;vertical-align:top; width: 300px;",
                       verbatimTextOutput("dirIn", placeholder = FALSE)),
                   #actionButton("uploadFile", "Upload Files"),
                   tags$h5("Result Directory:"),
                   div(style="display: inline-block;vertical-align:top; width: 75px;",
                       shinyDirButton("dirOut", "Choose...", "Upload")),
                   div(style="display: inline-block;vertical-align:top; width: 10px;",
                       HTML("<br>")),
                   div(style="display: inline-block;vertical-align:top; width: 300px;",
                       verbatimTextOutput("dirOut", placeholder = FALSE))
              ), #box
              uiOutput('resetable_input2')
            )#fluidRow

    ),#tabItem
    tabItem(tabName = "files",
            uiOutput("ui")
    )#tabItem
  )#tabItems

)#dashboard body
ui <- dashboardPage(header, sidebar, body)


server <- function(input, output, session) {
  observeEvent(input$quit, {
    print("Exiting FLOWMAP")
    stopApp()
  })
  output$resetable_input <- renderUI({
    times <- input$reset
    print(times)
    div(id=letters[(times %% length(letters)) + 1],
    fluidRow(
      # Boxes need to be put in a row (or column)
      #column(4.5,
      tabBox( width=9, id="settings_tabset",
              tabPanel(
                "Analysis Settings",
                #Distance Metric
                #select input
                useShinyalert(),  # Set up shinyalert
                div(style="display:inline-block",actionButton("DistanceMetricHelp", label = "?")),
                div(style="display:inline-block",
                    selectInput("distMetric", "Distance Metric:", c("Choose one" = "", "manhattan", "euclidean"))),
                #FLOW-MAP Mode
                #selectInput("flowmapMode", "FLOW-MAP Mode:", c("Choose one" = "", "multi", "single", "one")),

                #Subsample Number (Downsample Target Number)
                #numericInput
                numericInput("subsampNum", "Subsample Number (Downsample Target Number):", 200),
                #Number of Clusters
                numericInput("clusterNum", "Number of Clusters:", 100),
                #Set Seed Number
                numericInput("seedNum", "Set Seed Number:", 1),
                #Minimun Number of Edges
                numericInput("minEdgeNum", "Minimum Number of Edges:", 2),
                #Maximum Number of Edges
                numericInput("maxEdgeNum", "Maximum Number of Edges:", 5)
              ),#tabPanel
              #SPADE Downsampling:
              tabPanel(
                "Downsampling",
                checkboxInput("spade", "SPADE Downsampling", FALSE),
                conditionalPanel(condition = "input.spade == TRUE",
                                 numericInput("target.pctile", label = h5("Downsample Target Percentile"), value = 0.99),
                                 numericInput("exclude.pctile", label = h5("Downsample Exclude Percentile"), value = 0.01)
                )
              ),#tabPanel
              tabPanel(
                "Graph Output",
                #Color Palette:
                selectInput("colors", "Color Palette:", c("Choose one" = "", "bluered", "jet", "cb")),
                #Save Graph PDFs:
                #checkbox
                checkboxInput("saveGraphPDFs", "Save Graph PDFs", FALSE)
              )

              #,
              #actionButton("submitParams", "Submit")
      )##tabBox
    ),#column
    #column(4.5,
    fluidRow(
      box( # Input: Select a file ----
           title = "Directory Selection",
           #tags$hr(),
           tags$h5("Raw FCS Files or CSV Directory (for mode multiFLOW-MAP):"),
           div(style="display: inline-block;vertical-align:top; width: 75px;",
               shinyDirButton("dirIn", "Choose...", "Upload")),
           div(style="display: inline-block;vertical-align:top; width: 10px;",
               HTML("<br>")),
           div(style="display: inline-block;vertical-align:top; width: 300px;",
               verbatimTextOutput("dirIn", placeholder = FALSE)),
           #actionButton("uploadFile", "Upload Files"),
           tags$h5("Result Directory:"),
           div(style="display: inline-block;vertical-align:top; width: 75px;",
               shinyDirButton("dirOut", "Choose...", "Upload")),
           div(style="display: inline-block;vertical-align:top; width: 10px;",
               HTML("<br>")),
           div(style="display: inline-block;vertical-align:top; width: 300px;",
               verbatimTextOutput("dirOut", placeholder = FALSE))
      ), #box
      box(
        title = "FLOW-MAP Mode", collapsible = FALSE,

        selectInput("flowmapMode", "Select Analysis Mode:", c("Choose one" = "", "multi", "single", "one"))
      ),#box
      box( status = "primary",
           solidHeader = TRUE,

           # div(style="display: inline-block;vertical-align:top; width: 75px;",
           #     actionButton("submitParams", "Submit")),
           div(style = "background-color:gray; text-align:center;",
               actionButton("submitParams", "Submit")),
           # div(style="display: inline-block;vertical-align:top; width: 10px;",
           #     HTML("<br>")),
           # div(style="display: inline-block;vertical-align:top; width: 300px;",
           #     verbatimTextOutput("emptyParam", placeholder = FALSE))
           #uiOutput("MissingParams")
           verbatimTextOutput("emptyParam", placeholder = FALSE)
      )#box
      #)#column
    )#fluidRow
    )#div
  })#renderUI
#########################################################HELP FUNCTIONS
  # button functions
  observeEvent(input$RawFCSDirHelp, {
    # Show a modal when the button is pressed
    shinyalert("The directory that contains the raw FCS files.", type = "info")
  })
  observeEvent(input$DistanceMetricHelp, {
    # Show a modal when the button is pressed
    shinyalert("Select the appropriate distance metric.", type = "info")
  })
  observeEvent(input$ModeHelp, {
    # Show a modal when the button is pressed
    shinyalert("Method of analysis via multiple FCS files or via one.", type = "info")
  })
  observeEvent(input$ColorHelp, {
    # Show a modal when the button is pressed
    shinyalert("Color palette to use for saved PDFs of output graph.", type = "info")
  })
  observeEvent(input$ResultDirHelp, {
    # Show a modal when the button is pressed
    shinyalert("The directory where FLOW-MAP results will be saved.", type = "info")
  })
  #######
  observeEvent(input$DownsampleHelp, {
    # Show a modal when the button is pressed
    shinyalert("Turn on/off SPADE downsampling.", type = "info")
  })
  observeEvent(input$SavePDFsHelp, {
    # Show a modal when the button is pressed
    shinyalert("Turn on/off saving PDFs of output graph.", type = "info")
  })
  observeEvent(input$SubsampleNumHelp, {
    # Show a modal when the button is pressed
    shinyalert("The number of cells to sample from each FCS file. \
               If SPADE downsampling is selected, this field \
               specifies the target number.", type = "info")
  })
  observeEvent(input$ClusterNumHelp, {
    # Show a modal when the button is pressed
    shinyalert("The number of clusters to generate from each FCS file.", type = "info")
  })
  observeEvent(input$SeedNumHelp, {
    # Show a modal when the button is pressed
    shinyalert("The seed number for reproducible analysis.", type = "info")
  })

  observeEvent(input$EdgeNumMinHelp, {
    # Show a modal when the button is pressed
    shinyalert("The lower bound for number of edges during \
               density-dependent edge drawing steps.", type = "info")
  })
  observeEvent(input$EdgeNumMaxHelp, {
    # Show a modal when the button is pressed
    shinyalert("The upper bound for number of edges during \
               density-dependent edge drawing steps.", type = "info")
  })

######################################################### File Input tab server functionality
  #Allow for files as large as FCS file
  options(shiny.maxRequestSize=300*1024^2)
  # currdir <- getwd()
  # output$dirIn <- renderText({
  #   currdir
  # })
  # output$dirOut <- renderText({
  #   currdir
  # })

  #Access directory containing files to upload
  #need to use isolate so I can freeze each of them after they're selected
  roots = getVolumes()
  print(roots)
  shinyDirChoose(input, 'dirIn', roots = roots)
  global <- reactiveValues(datapath = getwd())
  dirIn <- reactive(input$dirIn)
  output$dirIn <- renderText({
    parseDirPath(roots, dirIn())
  })
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$dirIn
               },
               handlerExpr = {
                 home <- normalizePath("~")
                 global$datapath <-
                   file.path(home, paste(unlist(dirIn()$path[-1]), collapse = .Platform$file.sep))
               })
  #Access directory specifying desired results file location
  roots = getVolumes()
  shinyDirChoose(input, 'dirOut', roots = roots)
  global <- reactiveValues(datapath = getwd())
  dirOut <- reactive(input$dirOut)
  output$dirOut <- renderText({
    parseDirPath(roots, dirOut())
  })
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$dirOut
               },
               handlerExpr = {
                 home <- normalizePath("~")
                 global$datapath <-
                   file.path(home, paste(unlist(dirOut()$path[-1]), collapse = .Platform$file.sep))
                 #print(global$datapath)
                 #print(parseDirPath(roots, dirOut()))
               })
  #collect parameters to pass to flowmap algorithm
  params <- reactiveValues(inputs=list())
  globe.toggle <<- 0

  observeEvent(input$submitParams, {
    emptyParams <- ""
    params$inputs <- list()
    params$inputs[["mode"]] <- input$flowmapMode
    params$inputs[["color.palette"]] <- input$colors
    params$inputs[["downsample.toggle"]] <- input$spade
    params$inputs[["savePDFs.toggle"]] <- input$saveGraphPDFs
    params$inputs[["subsample.num"]] <- input$subsampNum
    params$inputs[["distance.metric"]] <- input$distMetric
    params$inputs[["cluster.num"]] <- input$clusterNum
    params$inputs[["seed.num"]] <- input$seedNum
    # inputs[["edge.pct.num"]] <- tclvalue(edge.pct.num)
    params$inputs[["edge.max.num"]] <- input$maxEdgeNum
    params$inputs[["edge.min.num"]] <- input$minEdgeNum

    for (item in names(params$inputs)){
      #print(item)
      #print(names(params$inputs))
      if (params$inputs[[item]] == ""){
        emptyParams <- paste(emptyParams, "\n", item)
      }
    }
    #print(emptyParams)
    if (emptyParams != "") {
      # output$MissingParams = renderUI({
      #     renderText(paste(" Missing parameters! \n",
      #               "Please make a selection for: \n", emptyParams))
      # })
      output$emptyParam <- renderText(paste(" Missing parameters! \n",
                                            "Please make a selection for: \n",emptyParams))
    } else {
      output$emptyParam <- renderText("Successful parameter and directory selection!")
      #inputs[["quit"]] <- quit.var
      globe.toggle <<- 1
      globe.inputs <<- params$inputs
      globe.raw.FCS.dir <<- parseDirPath(roots, dirIn())
      #print(globe.raw.FCS.dir)
      globe.result.dir <<- parseDirPath(roots, dirOut())
    }

  })##observeEvent


    #print(globe.inputs)
    #print(globe.raw.FCS.dir)
######################################################### Build file process tab
    output$ui <- renderUI({
      # build UI based on FLOW-MAP mode
      #print("globe.toggle")

      if (length(globe.inputs) != 0) {
        if (globe.inputs[["downsample.toggle"]] == TRUE) {
          if (globe.inputs[["mode"]] == "single") {
            ui <- fluidPage(
              titlePanel("File Uploader"),
              fluidRow(
                #column(width = 3,
                box( title="Upload and Order Files",
                     actionButton("loadFiles", "Load Files"),
                     selectInput("check.group.files",
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
                     textOutput("fcsorder")
                ),
                # box(
                #   numericInput("target.pctile", label = h5("Downsample Target Percentile"), value = 0.99),
                #   numericInput("exclude.pctile", label = h5("Downsample Exclude Percentile"), value = 0.01)
                # ),
                #),#column

                #column(width = 3,
                box(
                  selectInput("check.group.sim",
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
                  actionButton("merge.button", "Merge Selected Channels")
                ),
                #),#box, column
                #column(width = 5,
                box(
                  rHandsontableOutput("table", width = 600),
                  br(),
                  actionButton("start.button", "Run FLOWMAPR")
                )

              )#fluidRow
            )#fluidPage
            #)#
          } else if (globe.inputs[["mode"]] == "multi") {
            ui <- fluidPage(
              titlePanel("File Processing"),
              fluidRow(
                column(width = 3,
                       box(
                         selectInput("check.group.csv",
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
                         textOutput("fcsorder"))
                ),#box
                box(
                  numericInput("target.pctile", label = h5("Downsample Target Percentile"), value = 0.99),
                  numericInput("exclude.pctile", label = h5("Downsample Exclude Percentile"), value = 0.01)
                ),#box


                column(width = 3,
                       box(
                         selectInput("check.group.sim",
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
                         actionButton("merge.button", "Merge Selected Channels"))
                ),#box
                column(width = 5,
                       box(
                         rHandsontableOutput("table", width = 600)),
                       br(),
                       actionButton("start.button", "Run FLOWMAPR")
                )

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
                       rHandsontableOutput("table", width = 600)),
                br(),
                actionButton("start.button", "Run FLOWMAPR")
              )
            )
          }
        } else if (globe.inputs[["downsample.toggle"]] == FALSE) {
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
                       rHandsontableOutput("table", width = 600)),
                br(),
                actionButton("start.button", "Run FLOWMAPR")
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
                       rHandsontableOutput("table", width = 600)),
                br(),
                actionButton("start.button", "Run FLOWMAPR")
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
                       rHandsontableOutput("table", width = 600)),
                br(),
                actionButton("start.button", "Run FLOWMAPR")
              )
            )
          }
        } else {
          stop("Unrecognized downsampling selection!")
        }
      }#else if globe.toggle == 1
      else  {
        fluidRow(
          column(8, align="center",
                 icon("info-sign", lib = "glyphicon"),
                 textOutput("infoText")
          )#column
        )#fluidRow
      }#else
    })#renderUI

  #print info text regarding upload files page
  output$infoText <- renderText("FLOW-MAP Mode selection is required to load this page")
#######################################################################File process tab server functionality
  # build server based on FLOW-MAP mode
  if (length(globe.inputs) != 0) {
    print("1")
    if (globe.inputs[["mode"]] == "single") {
      print("single")
      print("observing...")
      observeEvent(input$loadFiles, {
        options(shiny.maxRequestSize = 1000 * 1024^2)
        panel.info <- InitializePanel()
        final.new.same <- NULL
        final.new.diff <- NULL
        print(globe.raw.FCS.dir)
        file.info <- FileOrder(globe.raw.FCS.dir)
        print(file.info)
        len.filenames <- file.info$len.filenames
        file.names <- file.info$file.names
        if (identical(file.names, character(0))) {
          choice <- "No FCS Files"
        } else {
          choice <- paste(len.filenames, file.names, sep = " ")
        }
        updateSelectInput(session, input$check.group.files,
                          choices = choice)
      })#observeEvent "loadFiles" button
      # observe({
      #   updateSelectInput(session, "check.group.files",
      #                     choices = choice)
      # })
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

    } else if (globe.inputs[["mode"]] == "multi") {

      if (globe.inputs[["quit"]]) {
        stopApp()
      }
      options(shiny.maxRequestSize = 1000 * 1024^2)
      panel.info <- InitializePanel()
      final.new.same <- NULL
      final.new.diff <- NULL
      file.info <- FileOrder(globe.raw.FCS.dir)
      print(file.info)
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

    } else if (globe.inputs[["mode"]] == "one") {

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
    }#server stuff for depending on mode
  }#if (globe.toggle == 1)
}



shinyApp(ui = ui, server = server)
