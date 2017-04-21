require(shiny)
require(rhandsontable)

shinyUI(
  fluidPage(
    titlePanel("File Uploader"),
    br(),
    br(),
    fluidRow(
      column(width = 3,
             actionButton("generbutton", "Generate Order List"),
             selectInput("checkGroup_files",
                         label = h5("Uploaded Order"),
                         choices = "Pending Upload",
                         selected = NULL,
                         multiple = TRUE,
                         selectize = FALSE,
                         size = 10
             ),
             textInput("fileorder", label = h5("Write the FCS File Order"), placeholder = "Ex: 4,2,7,5,3,1,6"),
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
                         label = h5("Similar Fields"),
                         choices = "Pending Upload",
                         selected = NULL,
                         multiple = TRUE,
                         selectize = FALSE,
                         size = 10
             ),
             selectInput("checkGroup_diff",
                         label = h5("Different Fields"),
                         choices = "Pending upload",
                         selected = NULL,
                         multiple = TRUE,
                         selectize = FALSE,
                         size = 10
             ),
             textInput("filemerge", label = h5("Select New Merge Name"), placeholder = "New Name"),
             actionButton("mbutton", "Merge Selected Diff")
      ),
      column(width = 5,
             actionButton("button", "Run FLOWMAP"),
             br(),
             rHandsontableOutput("table", width = 600)
             
      )
    )
  ))