require(shiny)
require(rhandsontable)
if(globe.inputs[["mode"]] == "single"){
  shinyUI(
    fluidPage(
      titlePanel("File Uploader"),
      fluidRow(
        column(width = 3,
               selectInput("check.group.files",
                           label = h5("Uploaded Order"),
                           choices = "Pending Upload",
                           selected = NULL,
                           multiple = TRUE,
                           selectize = FALSE,
                           size = 7
               ),
               textInput("file.order.input", label = h5("Write the FCS File Order"),
                         placeholder = "Ex: 4, 2, 7, 5, 3, 1, 6"),
               actionButton("gener.param.button", "Generate Parameters"),
               textOutput(
                 "TESTPRINT"
               ),
               textOutput(
                 "writefile"
               ),
               textOutput(
                 "vartable"
               ),
               textOutput(
                 "ordering"
               ),
               textOutput(
                 "fcsorder"
               )
        ),
        column(width = 3,
               selectInput("check.group.sim",
                           label = h5("Similar Fields"),
                           choices = "Pending Upload",
                           selected = NULL,
                           multiple = TRUE,
                           selectize = FALSE,
                           size = 7
               ),
               selectInput("check.group.diff",
                           label = h5("Different Fields"),
                           choices = "Pending upload",
                           selected = NULL,
                           multiple = TRUE,
                           selectize = FALSE,
                           size = 7
               ),
               textInput("file.merge", label = h5("Select New Merge Name"), placeholder = "New Name"),
               actionButton("merge.button", "Merge Selected Diff")
        ),
        column(width = 5,
               actionButton("start.button", "Run FLOWMAPR"),
               br(),
               br(),
               rHandsontableOutput("table", width = 600)
               # br(),
               # br(),
               # actionButton("down.gener.param", "Run Downsample")
        )
      )
    ))
} else if(globe.inputs[["mode"]] == "multi"){
  shinyUI(
    fluidPage(
      titlePanel("File Uploader"),
      fluidRow(
        column(width = 3,
               selectInput("check.group.csv",
                           label = h5("Import CSV"),
                           choices = "Pending Upload",
                           selected = NULL,
                           multiple = FALSE,
                           selectize = FALSE,
                           size = 3
               ),
               actionButton("csv.finder", "Input CSV"),
               br(),
               br(),
               selectInput("check.group.files",
                           label = h5("Uploaded Order"),
                           choices = "Pending Upload",
                           selected = NULL,
                           multiple = TRUE,
                           selectize = FALSE,
                           size = 7
               ),
               textInput("file.order.input", label = h5("Write the FCS File Order"),
                         placeholder = "Ex: 4, 2, 7, 5, 3, 1, 6"),
               actionButton("gener.param.button", "Generate Parameters"),
               # actionButton("default.button", "Use Default Order"),
               textOutput(
                 "TESTPRINT"
               ),
               textOutput(
                 "writefile"
               ),
               textOutput(
                 "vartable"
               ),
               textOutput(
                 "ordering"
               ),
               textOutput(
                 "fcsorder"
               )
        ),
        column(width = 3,
               selectInput("check.group.sim",
                           label = h5("Similar Fields"),
                           choices = "Pending Upload",
                           selected = NULL,
                           multiple = TRUE,
                           selectize = FALSE,
                           size = 7
               ),
               selectInput("check.group.diff",
                           label = h5("Different Fields"),
                           choices = "Pending upload",
                           selected = NULL,
                           multiple = TRUE,
                           selectize = FALSE,
                           size = 7
               ),
               textInput("file.merge", label = h5("Select New Merge Name"), placeholder = "New Name"),
               actionButton("merge.button", "Merge Selected Diff")
        ),
        column(width = 5,
               actionButton("start.button", "Run FLOWMAPR"),
               br(),
               br(),
               rHandsontableOutput("table", width = 600),
               br(),
               br(),
               actionButton("down.gener.param", "Run Downsample")
        )
      )
    ))
} else if(globe.inputs[["mode"]] == "one"){
  shinyUI(
    fluidPage(
      titlePanel("File Uploader"),
      fluidRow(
        column(width = 3,
               selectInput("check.group.csv",
                           label = h5("Import CSV"),
                           choices = "Pending Upload",
                           selected = NULL,
                           multiple = FALSE,
                           selectize = FALSE,
                           size = 3
               ),
               actionButton("csv.finder", "Input CSV"),
               br(),
               br(),
               selectInput("check.group.files",
                           label = h5("Uploaded Order"),
                           choices = "Pending Upload",
                           selected = NULL,
                           multiple = TRUE,
                           selectize = FALSE,
                           size = 7
               ),
               textInput("file.order.input", label = h5("Write the FCS File Order"),
                         placeholder = "Ex: 4, 2, 7, 5, 3, 1, 6"),
               actionButton("gener.param.button", "Generate Parameters"),
               # actionButton("default.button", "Use Default Order"),
               textOutput(
                 "TESTPRINT"
               ),
               textOutput(
                 "writefile"
               ),
               textOutput(
                 "vartable"
               ),
               textOutput(
                 "ordering"
               ),
               textOutput(
                 "fcsorder"
               )
        ),
        column(width = 3,
               selectInput("check.group.sim",
                           label = h5("Similar Fields"),
                           choices = "Pending Upload",
                           selected = NULL,
                           multiple = TRUE,
                           selectize = FALSE,
                           size = 7
               ),
               selectInput("check.group.diff",
                           label = h5("Different Fields"),
                           choices = "Pending upload",
                           selected = NULL,
                           multiple = TRUE,
                           selectize = FALSE,
                           size = 7
               ),
               textInput("file.merge", label = h5("Select New Merge Name"), placeholder = "New Name"),
               actionButton("merge.button", "Merge Selected Diff")
        ),
        column(width = 5,
               actionButton("start.button", "Run FLOWMAPR"),
               br(),
               br(),
               rHandsontableOutput("table", width = 600),
               br(),
               br(),
               actionButton("down.gener.param", "Run Downsample")
        )
      )
    ))
}
# exclude.pctile <- 0.01
# target.pctile <- 0.99