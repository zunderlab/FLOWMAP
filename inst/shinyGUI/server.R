require(flowCore)
require(shiny)

shinyServer(function(input, output, session) {
  # FLOW
  options(shiny.maxRequestSize = 1000 * 1024^2)
  # try(rm(count), silent = TRUE)
  # try(rm(final_new_diff), silent = TRUE)
  # try(rm(final_new_same), silent = TRUE)
  # try(rm(folder_now), silent = TRUE)
  # try(rm(FCS_filepath), silent = TRUE)
  # try(rm(order_sim), silent = TRUE)
  # try(rm(order_sim_list), silent = TRUE)
  # try(rm(order_while), silent = TRUE)
  # try(rm(operating_system), silent = TRUE)
  # try(rm(dir_now), silent = TRUE)
  # try(rm(len_filenames), silent = TRUE)
  # # DELETE ALL GLOBAL
  DF = data.frame(channels = c(NA), removal = c(NA), cluster = c(NA), annotate = c(NA))
  # order_while <<- 0
  operating_system <<- Sys.info()[1]
  print(operating_system)
  
  #get function for FLOW_MAP
  if (operating_system == "Windows"){
    folder_now <<- paste(gsub("/", "\\\\", getwd()), "\\", sep = "")
  } else{
    folder_now <<- paste(getwd(), "/", sep = "")
  }
  #get directory for where all flowmap function files are located
  # count <<- 0
  final_new_same <<- NULL
  final_new_diff <<- NULL
  
  #Set GLobal Variables
  dir_now = globe_resdir
  fileorder = eventReactive(
    input$generbutton, {
      file_names <<- list.files(dir_now, pattern = ".fcs")
      name_vec = c()
      for (i in 1:length(file_names)){
        name_vec = c(name_vec, i)
      }
      len_filenames <<- name_vec
      file_names
    })
  len_filen = eventReactive(
    input$generbutton, {
      file_names <<- list.files(dir_now, pattern = ".fcs")
      name_vec = c()
      for (i in 1:length(file_names)){
        name_vec = c(name_vec, i)
      }
      len_filenames <<- name_vec
      len_filenames
    })
  observe({
    updateSelectInput(session, "checkGroup_files", choices = paste(len_filen(), fileorder(), sep = " "))
  })
  chosen_order = eventReactive(input$generbutton2, {
    input$fileorder
  })
  fcs_order = eventReactive(input$generbutton2, {
    order = as.numeric(unlist(strsplit(chosen_order(), ",")))
    fcs_list = c()
    for(i in order){
      fcs_list = c(fcs_list, fileorder()[i])
      fcs_list
    }
  })
  contentdiff = eventReactive(input$generbutton2, {
    # Read input Files
    # Set the names
    print("HELLo")
    fcs_list = list()
    rows = length(len_filen())
    count = 0
    order = as.numeric(unlist(strsplit(chosen_order(), ",")))
    for(i in order)
    {
      setwd(dir_now)
      files = read.FCS(fileorder()[i], emptyValue = FALSE)
      filed = pData(parameters(files))[,c("name")]
      name_desc = do.call(paste, as.data.frame(filed, stringsAsFactors=FALSE))
      fcs_list[[(count + 1)]] = name_desc
      count = count + 1
      # Reads FCS Files, gets name and Description, add to a list of different FCS files
    }
    if(rows > 1){
      same = Reduce(intersect, fcs_list)
      every = Reduce(union, fcs_list)
      diff = every
      diff = diff[! every %in% same]
    }
    # Gets different parameters from the FCS files
    else{
      diff = NULL
    }
    final_new_diff <<- diff
    print(diff)
    diff
    #If there is 1 FCS file, then there is no difference
  })
  contentsame = eventReactive(input$generbutton2, {
    print("HELLO")
    fcs_list = list()
    rows = length(len_filen())
    count = 0
    order = as.numeric(unlist(strsplit(chosen_order(), ",")))
    for(i in order)
    {
      setwd(dir_now)
      files = read.FCS(fileorder()[i], emptyValue = FALSE)
      filed = pData(parameters(files))[,c("name")]
      name_desc = do.call(paste, as.data.frame(filed, stringsAsFactors=FALSE))
      fcs_list[[(count + 1)]] = name_desc
      count = count + 1
      #Does sdame thing as above
    }
    same = Reduce(intersect, fcs_list)
    every = Reduce(union, fcs_list)
    diff = every
    diff = diff[! every %in% same]
    final_new_same <<- same
    print(same)
    same
    #gives the same paramters
  })
  
  tablecreate = eventReactive(input$generbutton2, {
    if(length(final_new_diff) == 0){
      DF = data.frame(channels = c(final_new_same, final_new_diff), removal = logical(length = length(final_new_same)), cluster = logical(length = length(final_new_diff)+ length(final_new_same)), annotate = c(final_new_same, final_new_diff), stringsAsFactors = FALSE)
    }else{
      DF = data.frame(channels = c(final_new_same, final_new_diff), removal = c(logical(length = length(final_new_same)), !logical(length = length(final_new_diff))), cluster = logical(length = length(final_new_diff)+ length(final_new_same)), annotate = c(final_new_same, final_new_diff), stringsAsFactors = FALSE)
    }
    output$table = renderRHandsontable({
      rhandsontable(DF) %>%
        hot_col("channels", readOnly = TRUE)
    })
    print(DF)
  })
  
  observeEvent(input$generbutton2,{
    if(length(final_new_diff) == 0){
      DF = data.frame(channels = c(final_new_same, final_new_diff), removal = logical(length = length(final_new_same)), cluster = logical(length = length(final_new_diff)+ length(final_new_same)), annotate = c(final_new_same, final_new_diff))
    }else{
      DF = data.frame(channels = c(final_new_same, final_new_diff), removal = c(logical(length = length(final_new_same)), !logical(length = length(final_new_diff))), cluster = logical(length = length(final_new_diff)+ length(final_new_same)), annotate = c(final_new_same, final_new_diff))
    }
  })
  observe({
    updateSelectInput(session, "checkGroup_sim", choices = contentsame())
  })
  # updates the checkbox group to show same
  observe({
    updateSelectInput(session, "checkGroup_diff", choices = contentdiff())
  })
  write_file = eventReactive(input$button, {
    if(operating_system != "Windows"){
      dir = globe_resdir
    } else{
      dir = globe_resdir
    }
    setwd(dir)
    
    # Chooses export directory
    # COME BACK TO THIS
    # files_name = input$file_name
    # same_stuff = input$checkGroup_sim
    # diff_stuff = input$checkGroup_diff
    # if(files_name == ""){
    #   files_name = Sys.Date()
    #   files_name = paste("file", files_name, sep = "_")
    # }
    #creates file name
    # file_name_txt = paste(files_name, "txt", sep = ".")
    # sink(file_name_txt)
    # for(x in 1:length(same_stuff)){
    #   cat(same_stuff[x], "\n")
    # }
    # for(y in 1:length(diff_stuff)){
    #   cat(diff_stuff[y], "\n")
    # }
    #created the file to have same then different parameters
    # sink()
    # vec_name = paste(files_name, "vector", sep = "_")
    # vec_name_txt = paste(vec_name, "txt", sep = ".")
    # final = c(same_stuff, diff_stuff)
    # sink(vec_name_txt)
    # cat(final)
    # sink()
    
    # writes the file
    print(unlist(input$table))
    inputhello <<- input$table
    write.csv(unlist(input$table), paste0(gsub("/", "//", globe_resdir), "//", globe_proj, ".csv"))
    print(folder_now)
    set.seed(input$seed_num)
    print("works")
    setwd(dir_now)
    
    # Run FLOW-MAP
    # suppressWarnings(FLOWMAP(files = files, var.remove = var.remove, var.annotate = var.annotate,
    #                          clustering.var = clustering.var, cluster.numbers = cluster.numbers,
    #                          subsamples = subsamples, distance.metric = distance.metric,
    #                          minimum = minimum, maximum = maximum, per = per, save.folder = save.folder,
    #                          mode = mode, shuffle = TRUE, name.sort = FALSE, downsample = FALSE))
    
  })
  output$stuff = renderText({
    write_file()
    NULL
  })
  output$stuff3 = renderText({
    # foldern()
    tablecreate()
    NULL
  })
  output$stuff2 = renderText({
    chosen_order()
    NULL
  })
  output$stuff4 = renderText({
    fcs_order()
  })
})

#############################################################################################################
########## I STILL NEED TO FIGURE OUT HOW TO CHANGE FCS FILE, AND COPY OLD ONE AND CHANGE THE COPY ##########
#############################################################################################################
