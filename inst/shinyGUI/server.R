
shinyServer(function(input, output, session) {
  #FLOW
  options(shiny.maxRequestSize = 1000*1024^2)
  try(rm(count), silent = TRUE)
  try(rm(final_new_diff), silent = TRUE)
  try(rm(final_new_same), silent = TRUE)
  try(rm(folder_now), silent = TRUE)
  try(rm(FCS_filepath), silent = TRUE)
  try(rm(order_sim), silent = TRUE)
  try(rm(order_sim_list), silent = TRUE)
  try(rm(order_while), silent = TRUE)
  try(rm(operating_system), silent = TRUE)
  try(rm(dir_now), silent = TRUE)
  try(rm(len_filenames), silent = TRUE)
  #DELETE ALL GLOBAL
  DF = data.frame(channels = c(NA), removal = c(NA), cluster = c(NA), annotate = c(NA))
  # order_while <<- 0
  operating_system <<- Sys.info()[1]
  print(operating_system)
  
  #get function for FLOW_MAP
  
  if (operating_system == "Windows"){
    folder_now <<- paste(gsub("/", "\\\\", getwd()), "\\", sep = "")
  } else{
    paste(getwd(), "/", sep = "")
  }
  #get directory for where all flowmap function files are located
  
  # count <<- 0
  final_new_same <<- NULL
  final_new_diff <<- NULL
  
  #Set GLobal Variables
  
  # foldern = eventReactive(input$folderbutton, {
  #   if(operating_system != "Windows"){
  #     library(tcltk)
  #     dir_now <<- tclvalue(tkchooseDirectory())
  #   } else{
  #     dir_now <<- choose.dir(getwd(), "Choose a suitable folder")
  #   }
  # })
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
    #     inFile = input$file1
    #     if (is.null(inFile))
    #       return(NULL)
    # #Read input Files
    #     directory = getwd()
    #     names = inFile[ ,"name"]
    #Set the names
    library(flowCore)
    print("HELLo")
    fcs_list = list()
    rows = length(len_filen())
    count = 0
    order = as.numeric(unlist(strsplit(chosen_order(), ",")))
    for(i in order)
    {
      # file_num = as.character(count)
      # row_atm = inFile[i, ]
      # 
      # if (operating_system != "Windows"){
      # file_path = gsub("//", "/", row_atm[ ,"datapath"])
      # split_path = unlist(strsplit(file_path, "/"))
      # file_path = paste(head(split_path, n = (length(split_path) - 1)), collapse = "/")
      # 
      # } else{
      # split_res = strsplit(row_atm[ ,"datapath"], "/")                    #Comment out in Mac
      # file_path = paste(split_res[[1]][1], split_res[[1]][2], sep = "\\") #Comment out in Mac
      # }
      
      # setwd(file_path)
      #Sets path to temp loaction of FCS in Shiny
      # FCS_filepath <<- file_path
      # print(file_path)
      setwd(dir_now)
      files = read.FCS(fileorder()[i], emptyValue = FALSE)
      filed = pData(parameters(files))[,c("name")]
      name_desc = do.call(paste, as.data.frame(filed, stringsAsFactors=FALSE))
      fcs_list[[(count + 1)]] = name_desc
      count = count + 1
      #Reads FCS Files, gets name and Description, add to a list of different FCS files
    }
    if(rows > 1){
      same = Reduce(intersect, fcs_list)
      every = Reduce(union, fcs_list)
      diff = every
      diff = diff[! every %in% same]
    }
    #Gets different parameters from the FCS files
    else{
      diff = NULL
    }
    final_new_diff <<- diff
    print(diff)
    diff
    #If there is 1 FCS file, then there is no difference
  })
  contentsame = eventReactive(input$generbutton2, {
    # inFile = input$file1
    # if (is.null(inFile))
    #   return(NULL)
    # print(inFile)
    # directory = getwd()
    # names = inFile[ ,"name"]
    library(flowCore)
    print("HELLO")
    fcs_list = list()
    rows = length(len_filen())
    count = 0
    order = as.numeric(unlist(strsplit(chosen_order(), ",")))
    for(i in order)
    {
      
      # file_num = as.character(count)
      # row_atm = inFile[i, ]
      # if (operating_system != "Windows"){
      # #uncomment in Mac
      # file_path = gsub("//", "/", row_atm[ ,"datapath"])
      # split_path = unlist(strsplit(file_path, "/"))
      # file_path = paste(head(split_path, n = (length(split_path) - 1)), collapse = "/")
      # } else{
      # split_res = strsplit(row_atm[ ,"datapath"], "/")                    #Comment out in Mac
      # file_path = paste(split_res[[1]][1], split_res[[1]][2], sep = "\\") #Comment out in Mac
      # }
      # setwd(file_path)
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
    library(rhandsontable)
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
  # updates the checkbox group to show different
  # contentmerge = reactive({
  #   merge = input$checkGroup_diff
  #   merge
  # })
  #selects what is to be merged
  # observe({
  #   updateSelectInput(session, "checkGroup_merge", choices = contentmerge())
  # })
  #shows files seleteced to be merged
  write_file = eventReactive(input$button, {
    if(operating_system != "Windows"){
      library(tcltk)
      dir = globe_resdir
    } else{
      dir = globe_resdir
    }
    setwd(dir)
    #Chooses export directory
    
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
    #writes the file
    print(unlist(input$table))
    inputhello <<- input$table
    write.csv(unlist(input$table), paste0(gsub("/", "//", globe_resdir), "//", globe_proj, ".csv"))
    rm(final_new_same)
    rm(count)
    rm(final_new_diff)
    #removes most global variables
    print(folder_now)
    # set.seed(input$seed_num)
    print("works")
    setwd(dir_now)
    # flowmaster_run(prefolder = folder_now, single.folder = FCS_filepath, save.folder = dir, folder = FCS_filepath,
    #                file.format = "*.fcs", var.annotate = list("marker1" = "marker1", "marker2" = "marker2"),
    #                var.remove = c(), per = input$slider_pct, minimum = input$slider_maxmin[1], maximum = input$slider_maxmin[2],
    #                distance.metric = input$radio_dist, subsample = input$sample_num, cluster.number = input$cluster_num,
    #                seed.X = input$seed_num, clustering.var = c("marker1", "marker2"))
    # setup(folder_now)
    # print(input$checkGroup_sim)
    # SingleFLOWMAP(folder = dir_now, file.format = "*.fcs", var.remove = c(),
    #               var.annotate = list("marker1" = "marker1", "marker2" = "marker2"), clustering.var = input$checkGroup_sim,
    #               cluster.number = input$cluster_num, subsample = input$sample_num,
    #               distance.metric = input$radio_dist, minimum = input$slider_maxmin[1], maximum = input$slider_maxmin[2],
    #               per = input$slider_pct, save.folder = dir, shuffle = TRUE)
    #Run FLOW_Map
    
    # order_while <<- 1
    
    try(rm(operating_system), silent = TRUE)
    try(rm(dir_now), silent = TRUE)
    try(rm(len_filenames), silent = TRUE)
    
    #remove final global variables
  })
  # new_diff = eventReactive(input$button2,{
  #   diff_stuff = input$checkGroup_diff
  #   diff_stuff
  # })
  # #saves selected choices for different paramters
  # new_same = eventReactive(input$button2,{
  #   merge_stuff = input$checkGroup_merge
  #   merge_stuff
  # })
  # #saves selected parameters for merged stuff
  # final_same = eventReactive(input$button2, {
  #   length_same = length(final_new_same)
  #   final_new_same <<- append(final_new_same, new_same(), after = length_same + count)
  #   final_new_same
  # })
  # 
  # final_diff = eventReactive(input$button2, {
  #   final_new_diff <<- final_new_diff [! final_new_diff %in% new_diff()]
  #   final_new_diff
  # })
  # final_merge = eventReactive(input$button2, {
  #   if(identical(final_new_diff, character(0))){
  #     NULL
  #   }
  #   else{
  #     contentmerge()
  #   }
  # })
  # observe({
  #   updateSelectInput(session, "checkGroup_sim", choices = final_same())
  # })
  # observe({
  #   updateSelectInput(session, "checkGroup_diff", choices = final_diff())
  # })
  # observe({
  #   updateSelectInput(session, "checkGroup_merge", choices = final_merge())
  # })
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




