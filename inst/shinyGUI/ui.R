require(shiny)
require(rhandsontable)

shinyUI(
  fluidPage(
    titlePanel("File Uploader"),
    fluidRow(
      column(width = 3,
             selectInput("checkGroup_files",
                         label = h5("Uploaded Order"),
                         choices = "Pending Upload",
                         selected = NULL,
                         multiple = TRUE,
                         selectize = FALSE,
                         size = 10
             ),
             textInput("fileorder", label = h5("Write the FCS File Order"),
                       placeholder = "Ex: 4, 2, 7, 5, 3, 1, 6"),
             actionButton("generbutton", "Generate Parameters"),
             actionButton("defaultbutton", "Use Default Order"),
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
             selectInput("checkGroup_sim",
                         label = h5("Similar Fields"),
                         choices = "Pending Upload",
                         selected = NULL,
                         multiple = TRUE,
                         selectize = FALSE,
                         size = 7
             ),
             selectInput("checkGroup_diff",
                         label = h5("Different Fields"),
                         choices = "Pending upload",
                         selected = NULL,
                         multiple = TRUE,
                         selectize = FALSE,
                         size = 7
             ),
             textInput("filemerge", label = h5("Select New Merge Name"), placeholder = "New Name"),
             actionButton("mbutton", "Merge Selected Diff")
      ),
      column(width = 5,
             actionButton("button", "Run FLOWMAPR"),
             br(),
             br(),
             rHandsontableOutput("table", width = 600)
      )
    )
  ))

# exclude.pctile <- 0.01
# target.pctile <- 0.99