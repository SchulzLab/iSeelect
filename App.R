#load packages
library("org.Mm.eg.db")
library(BiocManager)
options(repos = BiocManager::repositories())
library(SingleCellExperiment)
library(SummarizedExperiment)
library(Biobase)
library(shiny)
library(irlba)
library(Rtsne)
library(shinyhelper)
library(magrittr)
library(BiocSingular)     
library(scater)
library(scran)
library(shinyalert)

#ui of iSEElect
ui <- fluidPage(
  #pop up boxes
  useShinyalert(),
  #title of the app
  titlePanel(
    fluidRow(
      column(3,div("iSEElect")%>% helper(type="inline", title = "Welcome to iSEElect!", colour = "red", size = "s",
                                         content = c("This App analysis your Expressiondata by converting your data into an RDS - file and running iSEE!",
                                                     "Have fun!"))))),
  # first phase of iSEElect, data input
  h4("Select your data for iSEE analysis"),
  column(10, wellPanel(
    fluidRow(
      #expression data upload
      column(3, fileInput("target_upload", h4("Expression Data"),
                          accept = c("text/csv","text/comma-separated-values",".csv"))
             %>% helper(type="inline", title = "Expression Data", colour = "red",
                        content = c("Please upload your data as CSV - file.",
                                    "You have to upload your data with the following <b>Layout</b>:",
                                    "<b>First row:</b>  Gene names or ID etc. (no duplicates)",
                                    "<b>Following rows:</b>  Time-steps (day_cell-type_replicate, d1_FB_4)",
                                    "<em>Note that you sort your daypoints by cell-types.</em>")),
             #separator
             radioButtons("separator","Separator: ",choices = c(";",",",":"), 
                          selected=";")
             %>% helper(type="inline", title = "Seperator CSV-file", colour = "red",
                        content = c("Please choose the seperator of your CSV - file.",
                                    "<b>,</b> or <b>;</b> or <b>:</b> is possible."))
      ),
      #phenotypic data input
      column(3, fileInput("row_data", h4("Phenotypic Data"),
                          accept = c("text/csv","text/comma-separated-values",".csv"))
             %>% helper(type="inline", title = "Phenotypic Data", colour = "red",
                        content = c("Please upload your <b>Phenotypic data</b> as CSV - file.",
                                    "You have to upload your data with the following <b>Layout</b>:",
                                    "<b>First row:</b>  Samples. (no duplicates)",
                                    "<b>Following rows:</b>  all other data (like dates, days, sex, etc.)",
                                    "<em>Note that length of row.names(Pheno) == length of column.names(Data) has to be TRUE.</em>")),
             #annotation
             radioButtons("annotation", "Annotation: ", choices = c("Human", "Mouse", "C.elegans", "Drosophila", "Other"))
      )),
    fluidRow(
      #get activated if both input files are selected
      conditionalPanel("output.fileUploaded && output.phenoUploaded",
                       #shows inputted data
                       column(3, actionButton(inputId = "show_table", label = h4("Show Data"))),
                       #starts second phase, condition selection
                       column(3, actionButton(inputId = "compute", label = h4("Start")
                       )%>% helper(type="inline", title = "Start Selection", colour = "red",
                                   content = c("<b>Start:</b> Start your Selection if you´ve uploaded the correct data.",
                                               "<b>show data table:</b> Controll your uploaded file by pressing show data table")
                       ))))
    
  ),
  #starts by clicking "Start", scond phase of iSEElect -> condition selection
  conditionalPanel("input.compute%2==1",
                   wellPanel(
                     fluidRow(
                       #selection of cell-type
                       column(3, selectInput(inputId = "cell_type", label = h4("Celltype/s"), choices = c("HC","FB","EC","CM"), selected = NULL, multiple = T)
                              %>% helper(type="inline", title = "Cell - type", colour = "red",
                                         content = c("Choose the Cell-types you want to have compared.", "Cannot be left empty!",
                                                     "<em>Dynamically adapted through your input!</em>")) 
                       ),
                       #selection of replicates
                       column(4, checkboxGroupInput(inputId = "replicate", label = h4("Replicate/s"), choices = c(1,2,3,4), inline = TRUE)
                              %>% helper(type="inline", title = "Replicates", colour = "red",
                                         content = c("Choose the replicates you want to have compared.","If you choose more than one Replicate, you can calculate the Median",
                                                     "Cannot be left empty!",
                                                     "<em>Dynamically adapted through your input!</em>")) 
                       )
                     ),
                     textOutput("text"),
                     fluidRow(
                       #selection of day points
                       column(3, selectInput(inputId = "day", label = h4("Timepoint/s"), choices = c(0,1,3,7,14,28), width = 150, multiple = T, selected = NULL)
                              %>% helper(type="inline", title = "Daypoints", colour = "red",
                                         content = c("Choose the Daypoints you want to have compared.", "Cannot be left empty",
                                                     "<em>Dynamically adapted through your input!</em>"))
                       )
                     ),
                     #select to use Median, of all selected replicates, possible if more than 1 replicate is selected
                     conditionalPanel("input.replicate.length > 1",
                                      htmlOutput("text2"),
                                      fluidRow(
                                        #median of all replicates
                                        column(3, checkboxInput(inputId = "median_plot", label = "use median", FALSE)
                                               %>% helper(type="inline", title = "Use Median", colour = "red",
                                                          content = c("<b>Selected:<b/> Use Median", "<b>Not selected:</b> Don´t use median",
                                                                      "Computes the median of your selected replicates.",
                                                                      "<em>You can change that at any time.</em>"))
                                        ),
                                        #median of all replicates + all single replicates
                                        column(4, checkboxInput(inputId = "median_rep", label = "median + replicates", FALSE)
                                               %>% helper(type="inline", title = "Use Median + Replicates", colour = "red",
                                                          content = c("<b>Selected:<b/> Use Median + single Replicates", "<b>Not selected:</b> Just Replicates",
                                                                      "Computes the median of your selected replicates. Displays Median + selected Replictaes.",
                                                                      "<em>You can change that at any time.</em>"))
                                        )
                                      )
                     ),
                     # buttons get activated when everything is selected (filled)
                     conditionalPanel("input.replicate.length > 0 && input.day.length > 0 && input.cell_type.length > 0",
                                      fluidRow(
                                        #show seleted data as data frame
                                        column(3, actionButton(inputId = "show_usertable", label = h4("Show Table"))
                                               %>% helper(type="inline", title = "iSEE", colour = "red",
                                                          content = c("That button shows your selection as data frame.", 
                                                                      "You will analyze THAT data by starting iSEE")) 
                                        ),
                                        #creates an SingleCellExperiment-object and downloads it as a RDS file
                                        column(4, downloadButton(outputId = "iSEE", label = h4("Create RDS"))
                                               %>% helper(type="inline", title = "iSEE", colour = "red",
                                                          content = c("That button creates the RDS - file.", 
                                                                      "You can analyze your data by starting iSEE",
                                                                      "Repress the Button to reproduce the RDS - file.")) 
                                        ),
                                        #pop up window to start iSEE manually on ones machine with created RDS
                                        column(5,
                                             actionButton(inputId = "GO", label = h4("start ISEE()"))
                                             %>% helper(type="inline", title = "iSEE", colour = "red",
                                                        content = c("That button starts the ineteractive shiny app iSEE.", 
                                                                    "You can analyze your data there.",
                                                                    "Restart iSEElect to reselect all the data.")) 
                                      )))
                    )
  ),
  
  #if one file input empty print "Please select a file or use default data."
  conditionalPanel("!(output.fileUploaded) && !(output.phenoUploaded)",
                   textOutput("please")),
  # if both input files are filled print "Please press start to start the UI."
  conditionalPanel("output.fileUploaded && output.phenoUpladed && input.compute%2==0",
                   textOutput("startit")), 
  # if everything is selected in phase 2 of iSEElect print "Your data is ready for your analysis."
  conditionalPanel("input.cell_type != '' && input.day != '' && input.replicate.length > 0", 
                   textOutput("result")),
  
  conditionalPanel("input.compute%2==1 && output.fileUploaded && output.phenoUpladed",
                   textOutput("prepare_data")),
  #prints 
  textOutput("data"),
  
  #prints "You are using the median for your analysis!" if checkbox median is selected
  textOutput("use_median"),
  
  #prints table
  column(dataTableOutput("table"), width = 12),
  
  #prints "You have started the UI succesfully!"
  textOutput("SCE")
  ))



server <- function(input, output, session){
  
  #increase maximum upload size
  shiny.maxRequestSize=35*1024^2
  
  #create help boxes (question marks)
  observe_helpers()
  
  #text outputs
  output$please <- renderText({"Please select a file or use default data."})
  output$startit <- renderText({"Please press start to start the UI."})
  output$result <- renderText({
    paste("Your data is ready for your analysis.")
  })
  output$text <- renderText({"Please choose days to compare your data!"})
  output$text2 <- renderUI({HTML(paste("You have picked more than one replicate.",
                                       "Decide if you want to use the median or not.", sep = "<br/>"))})
  #if median is selected
  output$use_median <- renderText({if(input$median_plot){"You are using the median for your analysis!"
  }else{""}
  })
  
  #function to create selectable conditions by reading column description, clicking the button compute is required
  observeEvent(input$compute, {
    observe({
      #get all column descriptions split up with "_"
      selections <- get.selection()
      #every condition as single vector
      reps <- as.vector(selections[[3]])
      cell <- as.vector(selections[[2]])
      days <- as.vector(selections[[1]])
      #update checkbox replicates with conditions
      updateCheckboxGroupInput(session, "replicate",
                               choices = reps,
                               inline=TRUE)
      #update selectinput cell-type with conditions
      updateSelectInput(session, "cell_type",
                        choices = cell)
      #update selectinput day points with conditions
      updateSelectInput(session, "day",
                        choices = days)
    })
    #print output that UI has started succesfully
    output$SCE <- renderText({
      paste("You have started the UI succesfully!")
    })
  })
  
  #function to get all conditions (cell-type, replicates, day points)
  get.selection <- reactive({
    #function to read csv expression data
    og <- get.data()
    #create empty vectors
    input_replicates <- numeric(0)
    input_days <- numeric(0)
    input_cell <- character(0)
    
    #iterate over all colnames of expression data
    for (i in colnames(og)){
      #split data at _
      splitted <- strsplit(i, "_")
      #append days (first condition) to day vector
      input_days <- append(input_days, as.numeric(substring(splitted[[1]][1], 2)))
      #append cell-type (second condition) to cell-type vector
      input_cell <- append(input_cell, splitted[[1]][2])
      #append replicates (third condition) to replicate vector
      input_replicates <- append(input_replicates, as.numeric(splitted[[1]][3]))
    }
    
    #delete duplicates and sort by size
    input_days <- input_days[!duplicated(input_days)]
    input_days <- sort(input_days)
    input_replicates <- input_replicates[!duplicated(input_replicates)]
    input_replicates <- sort(input_replicates)
    input_cell <- input_cell[!duplicated(input_cell)]
    
    #save as list and retur list
    new_list <- list(input_days, input_cell, input_replicates)
    return(new_list)
  })
  
  #function to read the inputted expression data
  get.data <- reactive({
    #inputted file
    inFile <- input$target_upload
    #if input is empty
    if (is.null(inFile)){
      df <- NULL
    #read the csv with the selected separator, decimal separator comma and row names 
    }else{
      df <- read.csv(inFile$datapath, header = TRUE,sep = input$separator, fill = TRUE, dec = ",", stringsAsFactors = FALSE, row.names = 1, as.is = TRUE)
    }
    return(df)
  })
  
  #function to read the inputted phenotypic data
  get.pheno <- reactive({
    #inputted file
    inFile <- input$row_data
    #if input is empty
    if (is.null(inFile)){
      pheno <- NULL
    #read the csv with the selected separator and row names 
    }else{
      pheno <- read.csv(inFile$datapath, header = TRUE,sep = input$separator, fill = TRUE, stringsAsFactors = FALSE, row.names = 1)
    }
    return(pheno)
  })
  
  #if button "show table" is clicked, display the inputted data as data frame
  observeEvent(input$show_table, {
    output$table<- renderDataTable({
      #get expression data
      df <- get.data()
      df
    }, options = list(scrollX = TRUE)) 
  })
  
  #if button "table" is clicked, the through specific conditions updated table gets displayed
  observeEvent(input$show_usertable, {
    if (input$median_plot){
      output$table<- renderDataTable({
        # function to calculate the median table, if median is used
        df <- get.median()
        df
      }, options = list(scrollX = TRUE))
    }
    else{
      output$table<- renderDataTable({
        #selected conditions w/o median
        df <- get.userdata()
        df
      }, options = list(scrollX = TRUE)) 
    }
  })
  
  #checks if gene expression file upload is empty or not
  output$fileUploaded <- reactive({
    return(!is.null(input$target_upload))
  })
  outputOptions(output, 'fileUploaded', suspendWhenHidden = FALSE)
  
  #checks if phenotypic file upload is empty
  output$phenoUploaded <- reactive({
    return(!is.null(input$row_data))
  })
  outputOptions(output, 'phenoUploaded', suspendWhenHidden = FALSE)
  
  #function to create the gene expression data of selected conditions
  get.userdata <- reactive({
    #read original csv
    og <- get.data()
    #create empty vector
    string_day <- character(0)
    #for every input create the sample name and append to created vector
    for (i in input$day){
      string_day_bef <- paste(sep="", "d", i, "_", rep(input$cell_type, each = length(input$replicate)), "_", input$replicate)
      string_day <- append(string_day, string_day_bef)
    }
    #creates vector
    true_dayta <- startsWith(colnames(og),string_day)
    #iterates over all created samples
    for (v in string_day){
      #creates vector and sets true if this sample is in original data frame else false
      true <- startsWith(colnames(og),v)
      #iterates from 0 to length of the vector
      for (n in (1:length(true))){
        #if the sample is true and sets True in created vector true_dayta
        if (true[n] == TRUE){
          true_dayta[n] <- TRUE
        }
      }
    }
    #creates empty vector
    user_colnames <- character(0)
    #iterates from 0 to the length of vector true_dayta
    for (i in (1:length(true_dayta))){
      #if the entry is true, append the sample to the created empty vector user_colnames
      if (true_dayta[i] == TRUE){
        user_colnames <- append(user_colnames, colnames(og)[i])
      }
    }
    #get all data from the original data frame, but just for the selected samples
    user_data_frame <- og[, user_colnames]
    return(user_data_frame)
  })
  
  #function to create the phenotypic data of selected conditions
  get.userphenodata <- reactive({
    #same procedure as get.data() but with rownames instead
    #get the original phenotypic data
    og_pheno <- get.pheno()
    string_day_ph <- character(0)
    for (i in input$day){
      string_day_ph_bef <- paste(sep="", "d", i, "_", rep(input$cell_type, each = length(input$replicate)), "_", input$replicate)
      string_day_ph <- append(string_day_ph, string_day_ph_bef)
    }
    
    true_dayta_ph <- startsWith(rownames(og_pheno),string_day_ph)
    for (s in string_day_ph){
      true_ph <- startsWith(rownames(og_pheno),s)
      for (l in (1:length(true_ph))){
        if (true_ph[l] == TRUE){
          true_dayta_ph[l] <- TRUE
        }
      }
    }
    user_rownames <- character(0)
    for (i in (1:length(true_dayta_ph))){
      if (true_dayta_ph[i] == TRUE){
        user_rownames <- append(user_rownames, rownames(og_pheno)[i])
      }
    }
    user_pheno_frame <- og_pheno[user_rownames,]
    return(user_pheno_frame)
  })
  
  #function to calculate the median for all samples and return the data frame
  get.median <- reactive({
    #get original data
    og_median <- get.data()
    #create empty vectors
    allmedians <- numeric(0)
    med_string_true <- character(0)
    med_rownames_before <- character(0)
    #create all strings but w/o replicate definition
    med_string <- paste(sep = "","d", input$day, "_", rep(input$cell_type, each = length(input$day)))
    #iterate over all selected day points
    for (i in input$day){
      #create samples for all day poitns
      med_rownames <- paste(sep="", "d", i, "_", rep(input$cell_type, each = length(input$replicate)), "_", input$replicate)
      for (z in med_rownames){
        if (z %in% colnames(og_median)){
          #if sample in original data append to vector
          med_rownames_before <- append(med_rownames_before, z)
        }
      }
    }
    #create table with selected samples
    med_table_before <- og_median[,med_rownames_before]
    
    
    for (i in med_string){
      #create vector with true for all med_string in the produced table (starts with)
      true_dayta1_med <- startsWith(colnames(med_table_before),i)
      #if the entry is true append to med_string_true
      if (TRUE %in% true_dayta1_med){
        med_string_true <- append(med_string_true, i)
      }
    }
    
    #iterate over all samples in med_string_true
    for (v in med_string_true){
      #create empty vector
      user_mednames <- character(0)
      #create vector with true for all med_string_true in the produced table (starts with)
      true_med <- startsWith(colnames(med_table_before),v)
      
      #iterate over all entries in true_med
      for (i in (1:length(true_med))){
        #if true append the sample to user_mednames
        if (true_med[i] == TRUE){
          user_mednames <- append(user_mednames, colnames(med_table_before)[i])
        }
      }
      #create the new table, with user_mednames
      user_median_table <- med_table_before[,user_mednames]
      
      #calculate the median 
      #if length of user med names 0 (just 1 replicate) calculate median for all samples and stor as vector
      if(length(user_mednames) >1){
        for(l in (1:length(med_table_before[,1]))){
          allmedians <- append(allmedians, (sum(user_median_table[l,])/length(user_median_table[l,])))
        }
      }
      else{
        #else split all replicates and save their values in a vector
        for (x in (1:length(med_table_before[,1]))){
          allmedians <- append(allmedians, user_median_table[x])
        }
      }
    }
    #create an empty matrix suitable for the median values
    output <- matrix(ncol = length(med_string_true), nrow = length(allmedians)/length(med_string_true))
    
    #fill the matrix with median values
    for (i in (1:length(allmedians))){
      output[i] <- allmedians[i]
    }
    
    #set rownames and colnames
    colnames(output) = med_string_true
    rownames(output) = rownames(user_median_table)
    #convert to data.frame
    df_output <- as.data.frame(output)
    
    return(df_output)
  })
  
  #function to get the phenotypic median data
  get.pheno_median <- reactive({
    #get selected pheno data
    og_pheno <- get.userphenodata()
    #get the calculated median expression data 
    og_med_data <- get.median()
    #get colnames of the median expression data
    names_col_pheno <- colnames(og_med_data)
    #create empty vector
    new_rownames <- character(0)
    
    #iterates over all colnames in the expression median data
    for (i in names_col_pheno){
      #create vector with true and fals (true if sample of mexpr. median data is in rownames of the user pheno data (starts with))
      right_lines <- startsWith(rownames(og_pheno), i)
      #iterate over created vector
      for (s in (1:length(right_lines))){
        #if true append sample and break for next replicate
        if (right_lines[s] == TRUE){
          new_rownames <- append(new_rownames, rownames(og_pheno)[s])
          break
        }
      }
    }
    #create new phenotypic median table
    new_phen_df <- og_pheno[new_rownames,]
    new_rownames_med <- character(0)
    #create and set rownames for that table
    for (i in new_rownames){
      new_rownames_med <- append(new_rownames_med, substr(i,1,5))
    }
    rownames(new_phen_df) <- new_rownames_med
    
    return(new_phen_df)
    
  })
  
  #function to calculate the input object for iSEE (SCE-object)
  get.download <- reactive({
    #get userdata 
     if (input$median_plot == FALSE){
       csv <- get.userdata()
       p_Data <- get.userphenodata()
     }
    #get median data, if selected
     else{
       csv <- get.median()
       p_Data <- get.pheno_median()
     }
  
     #create Matrix for Expressiondata
     exprs <- as.matrix(csv)
     #create Metadata
     meta_data <- data.frame(labelDescription =
                               c("Mouse Gender", "Sample Date", "tissue of Sample", "Day"),
                             row.names = colnames(p_Data))
  
     #create Phenodata
     pheno_data <- new("AnnotatedDataFrame", data = p_Data, varMetadata=meta_data)
  
     #create ExpressionSet ADD: User Input for annotation
     set <- ExpressionSet(assayData = exprs, phenoData = pheno_data, annotation = "org.Mm.eg.db")
  
     #create SE for iSEE
     se1 <- as(set, "SummarizedExperiment")
     #create SCE for iSEE (reduced dimension)
     sce1 <- as(se1, "SingleCellExperiment")
  
     #PCA
     sce1 <- scater::runPCA(sce1, exprs_values = "exprs")
  
     #TSNE
     irlba_out <- irlba(assay(sce1, "exprs"), nv=1)
     tsne_out <- Rtsne(irlba_out$v, pca = FALSE, perplexity = 0.33, verbose = TRUE)
     reducedDim(sce1, "TSNE") <- tsne_out$Y
  
     return(sce1)
   })
  
  #download the just created object as RDS
   output$iSEE <- downloadHandler(
     #file name with date
     filename = function(){
       paste("iSEE_object-", Sys.Date(), ".rds", sep="")
     },
     content = function(file){
       saveRDS(get.download(), file = file)
     }
   )
  
  #advices to start iSEE
  observeEvent(input$GO, {
    shinyalert("iSEE", "You have to start iSEE manually with created RSV - file in R. iSEE is not published.", type = "warning")
  })
}

#start the shiny app
shinyApp(ui = ui, server = server)
