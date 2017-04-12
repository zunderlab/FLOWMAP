
shinyUI(
  fluidPage(
    titlePanel("File Uploader"),
    br(),
    br(),
    fluidRow(
      column(width = 3,
             # fileInput('file1',label = h3("Choose File to Upload"), multiple = TRUE, accept = ".fcs"
             # ),
             # actionButton("folderbutton", "Choose folder of FCS Files"),
             # br(),
             # br(),
             actionButton("generbutton", "Generate Order List"),
             selectInput("checkGroup_files",
                         label = h3("Uploaded Order"),
                         choices = "Pending Upload",
                         selected = NULL,
                         multiple = TRUE,
                         selectize = FALSE,
                         size = 10
             ),
             textInput("fileorder", label = h3("Write the FCS File Order, if none, leave blank"), placeholder = "Ex: 4,2,7,5,3,1,6"),
             # actionButton("orderbutton", "Select Order"),
             actionButton("generbutton2", "Generate Parameters"),
             textOutput(
               "stuff"
             ),
             textOutput(
               "stuff2"
             ),
             textOutput(
               "stuff3"
             ),
             textOutput(
               "stuff4"
             )
             
             
             
      ),
      column(width = 3,
             
             br(),
             br(),
             selectInput("checkGroup_sim",
                         label = h3("Similar Fields"),
                         choices = "Pending Upload",
                         selected = NULL,
                         multiple = TRUE,
                         selectize = FALSE,
                         size = 10
             ),
             selectInput("checkGroup_diff",
                         label = h3("Different Fields"),
                         choices = "Pending upload",
                         selected = NULL,
                         multiple = TRUE,
                         selectize = FALSE,
                         size = 10
             )
      ),
      # column(width = 2,
      #        h4('If Needed, select fields to be merged From "Different Fields"; select 1.'
      #        ),
      #        selectInput("checkGroup_merge",
      #                    label = h3("Choose New Merge Name"),
      #                    choices = "Pending Merger Selection",
      #                    selected = NULL,
      #                    multiple = TRUE,
      #                    selectize = FALSE,
      #                    size = 3
      #        ),
      #        actionButton("button2", "Select Merger Name"),
      #        br(),
      #        br(),
      #        h4("Choose a file from the 'Similar Fields' area, and type a new name."),
      #        textInput("rename", label = h3("Rename Parameter"), value = "", placeholder = "Choose something to rename"),
      #        actionButton("button_rename", "Rename Parameter"),
      #        br(),
      #        selectInput("checkGroup_rename",
      #                    label = h3("Similar Fields with New Names"),
      #                    choices = "Pending upload",
      #                    selected = NULL,
      #                    multiple = TRUE,
      #                    selectize = FALSE,
      #                    size = 5
      #        )
      # ),
      
      
      
      
      # column(width = 2,
      #        radioButtons("radio_dist", label = h3("Select Distance Metric"),
      #                     choices = list("manhattan" = "manhattan", "euclidean" = "euclidean"), 
      #                     selected = "manhattan"),
      #        br(),
      #        h4('Enter FLOWMAP graph settings (if you are unsure what to use, try default settings)'),
      #        numericInput("slider_pct", label = h3("Percent of edges (ordered by smallest distance) to use to calculate density"), value = 1),
      #        sliderInput("slider_maxmin", label = h3("Minimum / Maximum number of edges to allot based on density"), min = 0, 
      #                    max = 100, value = c(3, 4)),
      #        br()
      # ),
      # column(width = 2,
      #        h4("Select number of cells to subsample from each FCS file and number of clusters for each FCS file (recommended ratio of 1:2, subsample to cluster)"),
      #        numericInput("sample_num", label = h3("Subsample Number"), value = 600),
      #        numericInput("cluster_num", label = h3("Cluster Number"), value = 300),
      #        br(),
      #        numericInput("seed_num", label = h3("Select the seed for reproducible analysis (change if you want a different result from a previous run)"), value = 1),
      #        radioButtons("radio_singmulti", label = h3("Do you want to run single or multi FLOWMAP"),
      #                     choices = list("multi-FLOWMAP" = 1, "single-FLOWMAP" = 2),
      #                     selected = "2")
      # ),
      
      
      
      
      column(width = 3,
             rHandsontableOutput("table", width = 600)
      ),
      column(width = 3,
             # textInput("file_name", label = h3("Export File"), value = "", placeholder = "Do Not Include Extension"),
             actionButton("button", "Write File")
      )
    )
  ))