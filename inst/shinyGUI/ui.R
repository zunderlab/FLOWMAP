
shinyUI(
  fluidPage(
    titlePanel("File Uploader"),
    br(),
    br(),
    fluidRow(
      column(width = 3,
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
      column(width = 3,
             rHandsontableOutput("table", width = 600)
      ),
      column(width = 3,
             actionButton("button", "Write File")
      )
    )
  ))
