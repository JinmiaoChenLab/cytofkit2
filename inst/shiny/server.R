library(shiny)
# library(shinyFiles)
# source("plot_functions.R")
# library(Seurat)
library(data.table)
source('./ui.R')
library(flowCore)
#library(cytofkit)
library(shinyalert)
source('./global.R')
library(shinyWidgets)

# library(ezTools)
# library(monocle)
# library(jsTree)

options(bitmapType='cairo')
# maximum size 10000MB
options(shiny.maxRequestSize = 1024*1024*100*100) 
# plan('multicore')

shinyServer = function(input, output, session)
{
  v = reactiveValues(
    # fcsFile = NULL
    markers = NULL
    , selected_markers = NULL
    , data = NULL
    , sampleInfo = NULL
    , sample_selected_index = NULL
    , sample_choices = NULL
    , sample_init = FALSE
    , sample_ready = FALSE
    , sample_update = FALSE
  )
  
  inputs = reactiveValues(
    fcsFiles = NULL
    , markers = NULL
    , l_w = 0.5
    , l_t = 500000
    , l_m = 4.5
    , l_a = 0
  )
  
  observeEvent(input$merge_method, {
    # browser()
    if(input$merge_method == "ceil" || input$merge_method == "fixed"){
      shinyjs::enable("fix_number")
    } else {
      shinyjs::disable('fix_number')
    }
  })
  
  observeEvent(input$rawfcs, {
    # browser()
    
    cur_dir = paste0(dirname(input$rawfcs$datapath)[1], "/")
    file.rename(input$rawfcs$datapath, paste0(cur_dir, input$rawfcs$name))
    inputs$fcsFiles = paste0(cur_dir, input$rawfcs$name)
    # input$rawfcs$datapath = paste0(cur_dir, input$rawfcs$name)
    
    fcs <- suppressWarnings(read.FCS(inputs$fcsFiles[1]))
    pd <- fcs@parameters@data
    markers <- paste(pd$name, "<", pd$desc, ">", sep = "")
    inputs$markers = markers
    updateSelectInput(session, "markers", choices = markers)
    # updatePickerInput(session, 'pickermarker', choices = markers)
    if (length(markers) == 0) {
      shinyalert(title = "Error", text = "No markers found in the FCS file!", type = "error")
      # stop()
    }
  })
  
  observeEvent(input$submit, {
    # browser()
    ### need to check data
    withProgress(
      {
        # browser()
        cytofkit(fcsFiles = inputs[["fcsFiles"]],
                 markers = input$markers,
                 projectName = input$project_name,
                 mergeMethod = input$merge_method,
                 fixedNum = input$fix_number,
                 transformMethod = input$transform_method,
                 dimReductionMethod = "tsne",
                 clusterMethods = input$cluster_method,
                 visualizationMethods = tolower(input$dr_method),
                 progressionMethod = input$progressionMethods,
                 Rphenograph_k = input$rphenograph_k,
                 FlowSOM_k = input$flowsom_k,
                 seed = input$seed,
                 clusterSampleSize = 500,
                 resultDir = tempdir(),
                 saveResults = TRUE,
                 saveObject = TRUE,
                 l_w = as.numeric(inputs[["l_w"]]), 
                 l_t = as.numeric(inputs[["l_t"]]), 
                 l_m = as.numeric(inputs[["l_m"]]), 
                 l_a = as.numeric(inputs[["l_a"]]))
      }
      
      , value = 0.5, message = "Analysing")
    if(!is.null(input$sample_anno)){
      res_fn = paste0(input$project_name, ".RData")
      load(res_fn)
      meta_data = read.table(input$sample_anno$datapath, sep = "\t", header = TRUE, row.names = 1)
      analysis_results$meta_data = meta_data
      analysis_results$sampel_name_ori = analysis_results$sampleNames
      # analysis_results$sampleInfo_ori = analysis_results$sampleInfo
      save(analysis_results, file = res_fn)
    }
  })
  
  output$download_analysis_res = downloadHandler(
    filename = function() {
      paste0(input$project_name, '.zip')
    },
    content = function(file) {
      # browser()
      cur_dir = './'
      files = dir(path = cur_dir, pattern = input$project_name)
      zip::zip(file, files)
      #stopApp(returnValue = invisible())
    }
  )
  
  ##------------------Reactive Values and Reactive Objects-------------------
  
  #if?
  c <- reactiveValues(clusterCol = list())
  p <- reactiveValues(progressionCluster = NULL)
  
  #parseQueryString to use .RData path as analysis results
  output$queryText <- renderText({
    query <- parseQueryString(session$clientData$url_search)
    if(!length(query) == 0){
      load(query[["fcspath"]])
      if(exists("analysis_results")){
        if(!is.null(analysis_results)) {
          v$data <- analysis_results
          v$sampleInfo <- data.frame(cellID = row.names(analysis_results$expressionData),
                                     cellSample = factor(sub("_[0-9]*$", "", row.names(analysis_results$expressionData))),
                                     stringsAsFactors = FALSE)
          p$progressionCluster <- names(analysis_results$clusterRes)[1]
          paste0("Loaded: ", query[["fcspath"]])
        }
      }
    }else{
      return(NULL)
    }
  })
  
  ## Scatter plot methods
  visualizationMethods <- reactive({
    if(is.null(v$data) || is.null(v$data$visualizationMethods)){
      return(NULL)
    }else{
      return(v$data$visualizationMethods)
    }
  })
  
  ## Scatter plot functions
  visualizationFunctions <- reactive({
    if(is.null(v$data) || is.null(v$data$clusterRes)){
      return(NULL)
    }else{
      return(c(names(v$data$clusterRes), 
               "Sample",
               "Density",
               "None"))
    }
  })
  
  ## cluster methods
  clusterMethods <- reactive({
    if(is.null(v$data))
      return(NULL)
    cMethods <- names(v$data$clusterRes)
    return(cMethods)
  })
  
  ## progression labs
  progressionLabs <- reactive({
    if(is.null(v$data))
      return(NULL)
    if(is.null(v$data$progressionRes))
      return(NULL)
    progressionLabs <- colnames(v$data$progressionRes[[3]])
    return(progressionLabs)
  })
  
  
  ##--------------------------------Side Panel-------------------------------
  
  ## Load cytofkit RData object
  observeEvent(input$goButton, {
    cytofkitObj <- input$cytofkitObj
    if (is.null(cytofkitObj)){
      v$data <- NULL
    }else{
      # browser()
      message(cytofkitObj$datapath)
      load(cytofkitObj$datapath)
      v$data <- analysis_results      
      
      if(is.null(v$data$projectName)){
        v$data$projectName <- "cytofkit_shinyAPP_output"
      }
      
      if(!is.null(v$data$progressionRes)){
        ## default the first cluster results are used for progression analysis
        p$progressionCluster <- names(v$data$clusterRes)[1]
      }
      
      
      # Need modification later
      # currently doesn't update sampleInfo with v$data$sampleInfo
      v$sampleInfo <- data.frame(cellID = row.names(v$data$expressionData),
                                 cellSample = factor(sub("_[0-9]*$", "", row.names(v$data$expressionData))),
                                 stringsAsFactors = FALSE)
      v$data$sampleInfo <- v$sampleInfo
      v$sample_ready = TRUE
      v$ori_sampleInfo = v$sampleInfo
      v$ori_data = v$data
      updateSelectInput(session, "sample_info_selection", choices = colnames(v$data$meta_data))
    }
  })
  
  ## For user, set roots option to your server directory 
  roots <- c(a=a)
  shinyFileChoose(input, 'serverObj', session = session, roots = roots, filetypes = "RData")
  
  observeEvent(input$serverObj, {
    inServer <- parseFilePaths(roots= roots, input$serverObj)
    print(inServer$datapath)
    load(as.character(inServer$datapath))
    v$data <- analysis_results
    if(is.null(v$data$projectName)){
      v$data$projectName <- "cytofkit_shinyAPP_output"
    }
    if(!is.null(v$data$progressionRes)){
      ## default the first cluster results are used for progression analysis
      p$progressionCluster <- names(v$data$clusterRes)[1]
    }
    # Need modification later
    # currently doesn't update sampleInfo with v$data$sampleInfo
    v$sampleInfo <- data.frame(cellID = row.names(v$data$expressionData),
                               cellSample = factor(sub("_[0-9]*$", "", row.names(v$data$expressionData))),
                               stringsAsFactors = FALSE)
    v$data$sampleInfo <- v$sampleInfo
  })
  
  output$rdata_desc <- renderText({
    if(is.null(v$data)){
      paste0("No .RData loaded yet")
    }else if(!length(session$clientData$url_search) == 0){
      return(NULL) 
    }else{
      paste0("Loaded: ", v$data$resultDir, "/", v$data$projectName, ".RData")
    }
  })
  
  observeEvent(input$reset, {
    analysis_results <- NULL
    session$reload()
    print("Reset done")
  })
  
  output$selectAll <- renderUI({
    if(is.null(v$data) || is.null(v$sampleInfo)){
      return(NULL)
    }else{
      checkboxInput('selectDeselectAll', label = "Select/Deselect All", value = TRUE)
    }   
  })
    
  #### only initial UI once.
  observeEvent(v$sample_ready, {
    output$sampleSelect <- renderUI({
      if(is.null(isolate(v$data)) || is.null(isolate(v$sampleInfo))){
        return(NULL)
      }else{
        # browser()
        sampleNames <- isolate(unique(as.character(v$sampleInfo$cellSample)))
        v$sample_choices = sampleNames
        v$sample_selected_index = 1:length(sampleNames)
        v$sample_init = TRUE
        checkboxGroupInput(inputId = 'samples', label = NULL,
                           choices = sampleNames, selected = sampleNames)
      }
    })
    if(v$sample_init){
      # browser()
      update_samples()
    }
  })
  
  
  observeEvent(unique(as.character(v$sampleInfo$cellSample)), {
    # browser()
    if(v$sample_init){
      # browser()
      sampleNames = unique(as.character(v$sampleInfo$cellSample))
      ### get selected item index
      # selected_index = match(input$samples, v$sample_choices)
      # v$sample_selected_index = selected_index
      v$sample_selected_index = 1:length(sampleNames)
      v$sample_choices = sampleNames
      updateCheckboxGroupInput(session, 'samples', choices = v$sample_choices)
      # browser()
      v$sample_update = !(v$sample_update)
      # updateCheckboxGroupInput(session, 'samples', selected = v$sample_choices[selected_index])
    }
  })
  
  observeEvent(v$sample_update, {
    # browser()
    if (length(v$sample_selected_index) > length(v$sample_choices)){
      v$sample_selected_index = 1:length(v$sample_choices)
    }
    updateCheckboxGroupInput(session, 'samples', selected = v$sample_choices[v$sample_selected_index])
  })
  
  update_samples = function(){
    # v$sampleInfo$cellSample = v$data$meta_data[match(v$data$sampleInfo_ori$cellSample, rownames(v$data$meta_data))
                                               # , input$sample_info_selection]
    # browser()
    allSamp <- input$selectDeselectAll
    sampleNames <- isolate(unique(as.character(v$sampleInfo$cellSample)))
    if(allSamp == TRUE){
      updateCheckboxGroupInput(session, 'samples', selected = sampleNames)
    }else{
      updateCheckboxGroupInput(session, 'samples', selected = character(0))
    }
    # browser()
  }
  observeEvent(input$selectDeselectAll, {
    # browser()
    update_samples()
  })
  
  observeEvent(input$sample_info_selection, {

    if(!is.null(v$data) && !is.null(input$sample_info_selection)){
      # browser()
      
      temp_name1 = sapply(1:length(v$data$sampleNames), function(i){
        v$data$sampleNames[[i]][1]
      })
      temp_name2 = sapply(1:length(v$data$sampleNames), function(i){
        v$data$sampleNames[[i]][2]
      })
      new_sample_names = as.character(v$data$meta_data[temp_name1[order(temp_name2)], input$sample_info_selection])
      reset_sample()
      rename_sample(new_sample_names)
      # v$sample_update = !v$sample_update
      # browser()
      # v$data$sampleNames
      # v$data$sampleNames = v$data$mea_data[unlist(v$data$sampel_name_ori), input$sample_info_selection]
      # v$data$sampleInfo$cellSample = v$data$meta_data[match(v$data$sampleInfo_ori$cellSample, rownames(v$data$meta_data))
                                                      # , input$sample_info_selection]
      # update_samples()
      # analysis_results = v$data
      # analysis_results$sampleInfo_ori = analysis_results$sampleInfo
      # save(analysis_results, file = "test2.Rdata")
    }
  })
  
  observe({
    if(!is.null(v$data) && !is.null(v$sampleInfo) && !is.null(input$samples)){
      x <- input$samples
      sampleNames <- unique(as.character(v$sampleInfo$cellSample))
      if(length(x) == 0){
        x <- character(0)
        updateCheckboxInput(session, 'selectDeselectAll', value = FALSE)
      }
      if(length(x) == length(sampleNames)){
        updateCheckboxInput(session, 'selectDeselectAll', value = TRUE)
      }
    }
  })
  
  output$summaryText1 <- renderText({
    if(is.null(v$data))
      return(NULL)
    paste0("-- ", nrow(v$data[[1]]), " cells x ", ncol(v$data[[1]]), " markers")
  })
  
  output$summaryText2 <- renderText({
    if(is.null(v$data))
      return(NULL)
    paste0("-- ", paste(names(v$data$clusterRes), collapse = " | "))
  })
  
  output$summaryText3 <- renderText({
    if(is.null(v$data))
      return(NULL)
    paste0("-- ", paste(v$data$visualizationMethods, collapse =  " | "))
  })
  
  output$summaryText4 <- renderText({
    if(is.null(v$data))
      return(NULL)
    paste0("-- ", ifelse(is.null(v$data$progressionRes), "NULL", 
                         sub("_[0-9]*$", "", colnames(v$data$progressionRes$progressionData)[1])))
  })
  
  output$summaryText5 <- renderText({
    if(is.null(v$data))
      return(NULL)
    paste0("-- ", paste(v$data$dimRedMarkers, collapse =  " | "))
  })
  
  # ## Save and parse cytofkit RData object
  # observeEvent(input$saveButton, {
  #   if (!is.null(v$data)){
  #     withProgress(message='Saving Results ', value=0, {
  #       ## check results saving path
  #       if(is.null(v$data$resultDir) || !dir.exists(v$data$resultDir)){
  #         v$data$resultDir <- path.expand("~")  ## default save to home if not specified
  #       }
  #       saveToFCS <- input$saveFCS
  #       if(is.null(v$data$rawFCSdir)){
  #         saveToFCS <- FALSE
  #         warning("Path for original FCS files is not provided, 
  #                 data cannnot be saved to new copies of FCS files.")
  #       }else if(!dir.exists(v$data$rawFCSdir)){
  #         saveToFCS <- FALSE
  #         warning(paste0("Path for original FCS files doesn't exist, 
  #                        data cannnot be saved to new copies of FCS files.", 
  #                        "Please check path: ", v$data$rawFCSdir))
  #       }
  #       
  #       ## NOTE: if samples are regrouped, then new FCS file cannot be saved
  #       incProgress(1/2, message = paste0("To ", v$data$resultDir))
  #       v$data$sampleInfo <- v$sampleInfo
  #       analysis_results <<- v$data
  #       cytof_writeResults(analysis_results,
  #                          saveToRData = input$saveRData,
  #                          saveToFCS = saveToFCS,
  #                          saveToFiles = input$saveCsv)
  #       incProgress(1/2)
  #       ## open the results directory
  #       opendir(v$data$resultDir)
  #     })
  #     }
  # })
  
  output$saveButton = downloadHandler(
    filename = function() {
      paste0(input$project_name, '_result.zip')
    },
    content = function(file) {
      # browser()
      withProgress(message='Saving Results ', value=0, {

        # ## check results saving path
        # if(is.null(v$data$resultDir) || !dir.exists(v$data$resultDir)){
        #   v$data$resultDir <- path.expand("~")  ## default save to home if not specified
        # }
        # saveToFCS <- input$saveFCS
        # if(is.null(v$data$rawFCSdir)){
        #   saveToFCS <- FALSE
        #   warning("Path for original FCS files is not provided, 
        #           data cannnot be saved to new copies of FCS files.")
        # }else if(!dir.exists(v$data$rawFCSdir)){
        #   saveToFCS <- FALSE
        #   warning(paste0("Path for original FCS files doesn't exist, 
        #                  data cannnot be saved to new copies of FCS files.", 
        #                  "Please check path: ", v$data$rawFCSdir))
        # }

        # browser()
        res_folder = paste0(input$project_name, "_results")
        if(!dir.exists(res_folder)){
          dir.create(res_folder, recursive = TRUE)
        }
        ## NOTE: if samples are regrouped, then new FCS file cannot be saved
        incProgress(1/2)
        v$data$sampleInfo <- v$sampleInfo
        analysis_results <<- v$data
        cytof_writeResults(analysis_results,
                           saveToRData = input$saveRData,
                           saveToFCS = FALSE,
                           saveToFiles = input$saveCsv
                           , resultDir = res_folder)
        incProgress(1/2)
        # ## open the results directory
        # opendir(v$data$resultDir)
        files = dir(path = res_folder, pattern = ".*", full.names = TRUE)
        zip::zip(file, files)
        #stopApp(returnValue = invisible())
      })
      
    }
  )
  
  output$reportButton = downloadHandler(
    filename = function() {
      paste0(input$project_name, '_report.html')
    },
    content = function(file) {
      # browser()
      withProgress(message='Generating report ', value=0, {
        
        # browser()
        # library(ezTools)
        # library(cytofkit2)
        
        #### suppress warning
        options(warn=-1)
        analysis_results = v$data
        
        rownames(analysis_results$sampleInfo) = analysis_results$sampleInfo$cellID
        analysis_results$clusterRes$Rphenograph = as.data.frame(analysis_results$clusterRes$Rphenograph)
        colnames(analysis_results$clusterRes$Rphenograph) = "Rphenograph Cluster"
        #### set output markdown file
        output_file_name = './cytofkit_report.RMD'
        
        create_script = function(title = "", script = "", values = NULL, fig_width = NULL, fig_height = NULL
                                 , chunk_option = NULL){
          # browser()
          param_names = as.list(match.call()$script)
          param_names = param_names[-1]
          r_option = "```{r, echo = FALSE"
          if (!is.null(chunk_option)) {
            r_option = paste0(r_option, ", ", chunk_option)
          }
          if (!is.null(fig_width)) {
            r_option = paste0(r_option, ", fig.width = ", fig_width)
          }
          if (!is.null(fig_height)) {
            r_option = paste0(r_option, ", fig.height = ", fig_height)
          }
          r_option = paste0(r_option, "}")
          res = paste(title, r_option, paste(param_names, collapse = "\n"), "```", sep = "\n")
          if (!is.null(values)) {
            for (i in 1:length(values)) {
              res = gsub(values[i], paste0("\"", eval_string(values[i], envir = parent.frame()), "\""), res, fixed = TRUE)
            }
          }
          res
        }
        
        substitute_script = function(script_file = "", script_id = "", scripts = "", output_file = ""){
          # browser()
          all_script = read_string(script_file)
          script_id = paste0("####* ", script_id, " *####")
          #### substitue script
          script = paste(scripts, collapse = "\n")
          script_res = gsub(script_id, script, all_script, fixed = TRUE)
          write_string(script_res, output_name = output_file)
        }
        
        create_dr_scripts = function(analysis_results){
          # browser()
          dr_names = names(analysis_results$dimReducedRes)
          scripts = lapply(1:length(dr_names), function(i){
            if (dr_names[i] == "tsne") {
              return("")
              # temp_name = 'tSNE'
            } else if (dr_names[i] == "umap") {
              temp_name = 'UMAP'
            } else {
              temp_name = dr_names[i]
            }
            
            analysis_results$temp1 <<- str_replace_all(analysis_results$dimRedMarkers, "<", "&lt;")
            analysis_results$temp2 <<- str_replace_all(analysis_results$temp1, ">", "&gt;")
            script1 = create_script(paste0("## ", temp_name, " plot color by sample\n", "Based on the markers: "
                                        , paste0(analysis_results$temp2, collapse = ', ')), {
                                          plot_scatter(analysis_results$dimReducedRes[[dr_names[i]]]
                                                       , analysis_results$sampleInfo[, "cellSample", drop = FALSE]) + coord_fixed()
                                        }, values = c("dr_names[i]"))
            script2 = create_script(paste0("## ", temp_name, " plot color by sample (splitted version)"), {
              p = plot_split_scatter(analysis_results$dimReducedRes[[dr_names[i]]]
                                     , analysis_results$sampleInfo[, "cellSample", drop = FALSE], ncol = 2, show_legend = FALSE)
              p[[1]]
            }, values = c("dr_names[i]"), fig_width = 8, fig_height = 20)
            forplot3 <<- as.data.frame(cbind(analysis_results$dimReducedRes[[dr_names[i]]], analysis_results$expressionData[,analysis_results$dimRedMarkers, drop = FALSE]))
            script3 = create_script(paste0("## ", temp_name, " plot color by marker expression"), {
              for (x in 1:length(analysis_results$dimRedMarkers)){
                print(ggplot(forplot3, aes(x = umap_1, y = umap_2)) + geom_point(aes(color = get(colnames(forplot3)[x+2]))) + 
                  scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07")) +
                  theme(legend.position = "right") +
                  coord_fixed() +
                  labs(colour = colnames(forplot3)[x+2]))
              }                               
            }, values = c("dr_names[i]"))
            return(c(script1, script2, script3))
          })
          unlist(scripts)
        }
        
        create_cluster_scripts = function(analysis_results){
          # browser()
          dr_names = names(analysis_results$dimReducedRes)
          scripts = lapply(1:length(dr_names), function(i){
            if (dr_names[i] == "tsne") {
              return("")
              # temp_name = 'tSNE'
            } else if (dr_names[i] == "umap") {
              temp_name = 'UMAP'
            } else {
              temp_name = dr_names[i]
            }
            script1 = create_script(paste0("## ", temp_name, " plot color by cluster"), {
              plot_scatter(analysis_results$dimReducedRes[[dr_names[i]]]
                           , analysis_results$clusterRes$Rphenograph) + coord_fixed()
            }, values = c("dr_names[i]"))
            # fig_height = ceiling(length(unique(analysis_results$clusterRes$Rphenograph[, 1]))/4)*2.5
            # script2 = create_script(paste0("## ", temp_name, " plot color by cluster (splitted version)"), {
            #   p = plot_split_scatter(analysis_results$dimReducedRes[[dr_names[i]]]
            #                          , analysis_results$clusterRes$Rphenograph, ncol = 4, show_legend = F) 
            #   p[[1]]
            # }, values = c("dr_names[i]"), fig_height = fig_height, fig_width = 9)
            # return(c(script1, script2))
            return(c(script1))
          })
          script2 = create_script(paste0("## ", " Percentage heatmap"), {
            
            temp = as.data.frame(analysis_results$clusterRes$Rphenograph)
            rownames(analysis_results$sampleInfo) = analysis_results$sampleInfo[, 1]
            temp = ezcbind(temp, analysis_results$sampleInfo[, "cellSample", drop = FALSE])
            freq = as.data.frame.matrix(table(temp))
            freq = freq/rowSums(freq)
            pheatmap(freq, silent = FALSE)
            # print(p)
          })
          unlist(scripts)
          c(scripts, script2)
        }
        
        create_markers_script = function(analysis_results){
          dr_names = names(analysis_results$dimReducedRes)
          scripts = lapply(1:1, function(i){
            if (dr_names[i] == "tsne") {
              return("")
              # temp_name = 'tSNE'
            } else if (dr_names[i] == "umap") {
              temp_name = 'UMAP'
            } else {
              temp_name = dr_names[i]
            }
            
            # selected_markers = which(!(grepl("NA", colnames(analysis_results$expressionData))))
            # fig_height = ceiling(length(unique(selected_markers))/6)*2.5
            script2 = create_script(paste0("## ", temp_name, " plot color by cluster (splitted version)"), {
              selected_markers = which(!(grepl("NA", colnames(analysis_results$expressionData))))
              marker_list = ez_chunk(selected_markers, ceiling(length(selected_markers)/3/7))
              lapply(1:length(marker_list), function(k){
                all_plots = lapply(1:length(marker_list[[k]]), function(j){
                  plot_scatter(analysis_results$dimReducedRes[[dr_names[i]]]
                               , analysis_results$expressionData[, marker_list[[k]][j], drop = FALSE]
                               , colors = c("#BEBEBE",brewer.pal(9,"Reds")), color_as_factor = FALSE) +
                    coord_fixed() 
                })
                p = plot_grid(plotlist = all_plots, ncol = 3)
                p
              })
            }, values = c("dr_names[i]"), fig_height = 15, fig_width = 9
            , chunk_option = "results='hide', fig.keep='all', message=FALSE")
            return(c(script2))
          })
          unlist(scripts)
        }
        
        create_express_heatmap = function(analysis_results){
          res = create_script(paste0("## Expression heatmap"), {
            selected_markers = which(!(grepl("NA", colnames(analysis_results$expressionData))))
            expression = as.data.frame(analysis_results$expressionData[, selected_markers, drop = FALSE])
            expression$clusters = analysis_results$clusterRes$Rphenograph[, 1]
            cluster_expression = fast_aggr(expression, ncol(expression))
            dt = seurat_sacle_data(cluster_expression)
            pheatmap(dt)
          }, chunk_option = "results='hide', fig.keep='all', message=FALSE")
          res
          
          
        }
        
        create_expression_histogram = function(analysis_results){
          script2 = create_script(paste0("## Expression histograms"), {
            dt = as.data.frame(analysis_results$expressionData)
            selected_markers = which(!(grepl("NA", colnames(analysis_results$expressionData))))
            selected_markers = colnames(analysis_results$expressionData)[selected_markers]
            marker_list = ez_chunk(selected_markers, ceiling(length(selected_markers)/3/6))
            invisible(lapply(1:length(marker_list), function(j){
              print(stackDenistyPlot(dt, marker_list[[j]], stackFactor = analysis_results$clusterRes$Rphenograph[, 1]))
            }))
          }, fig_width = 9, fig_height = 15, chunk_option = "results='hide', fig.keep='all', message=FALSE")
          return(c(script2))
        }
        
        create_abstract_scripts = function(analysis_results){
          paste0("The project \"", analysis_results$projectName, "\" has ")
          sample_num = length(unique(analysis_results$sampleInfo$cellSample))
          cell_num = nrow(analysis_results$sampleInfo)
          cluster_num = length(unique(analysis_results$clusterRes$Rphenograph$`Rphenograph Cluster`))
          res = paste0("There are ", sample_num, " sample", ifelse(sample_num > 1, "s", "")
                    , ", ", cluster_num, " cluster", ifelse(cluster_num > 1, "s", "")
                    , ", ", cell_num, " cell", ifelse(cluster_num > 1, "s", "")
                    , " in the project \"", analysis_results$projectName, "\".\n\n")
          res = paste0(res, "  The dimensionality reduction, clustering and markers expression anaylysis were conducted on the project data.")
          res
        }
        
        
        dr_scripts = create_dr_scripts(analysis_results)
        substitute_script('./pdf_report_template.Rmd', script_id = "DR analysis"
                          , scripts = dr_scripts
                          , output_file = output_file_name)
        
        abstract_scripts = create_abstract_scripts(analysis_results)
        substitute_script(output_file_name, script_id = "Abstract"
                          , scripts = abstract_scripts
                          , output_file = output_file_name)
        # render(output_file_name)

        cluster_scripts = create_cluster_scripts(analysis_results)
        substitute_script(output_file_name, script_id = "Cluster analysis"
                          , scripts = cluster_scripts
                          , output_file = output_file_name)

        markers_scripts = create_markers_script(analysis_results)
        substitute_script(output_file_name, script_id = "Expression on DR"
                          , scripts = markers_scripts
                          , output_file = output_file_name)

        expression_heatmap_scrip = create_express_heatmap(analysis_results)
        substitute_script(output_file_name, script_id = "Expression heatmap"
                          , scripts = expression_heatmap_scrip
                          , output_file = output_file_name)

        expression_histogram_script = create_expression_histogram(analysis_results)
        substitute_script(output_file_name, script_id = "Expression histogram"
                          , scripts = expression_histogram_script
                          , output_file = output_file_name)
        
        render(output_file_name, output_file = file)
        # output_pdf = basename(output_file_name)
        # output_pdf = paste0(get_file_name(output_pdf, with_ext = F), ".pdf")
        # system(paste0("cp ", output_file_name, " ", file))
        
      })
      
    }
  )
  
  # 
  # observeEvent(input$OpenDir, {
  #   pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
  #   if(dir.exists(pdfDir)){
  #     opendir(pdfDir)
  #   }else{
  #     stop("PDF not created yet!")
  #   }
  # })
  
  output$logo <- renderImage({
    return(list(
      src = "vignettes/logo.png",
      contentType = "image/png",
      alt = "Singapore Immunology Network"
    ))
  }, deleteFile = FALSE)
  
  ##------------------------------Cluster Panel------------------------------
  
  ##-----cluster plot-----
  output$C_PlotMethod <- renderUI({
    if(is.null(v$data) || is.null(visualizationMethods())){
      return(NULL)
    }else{
      selectInput('c_PlotMethod', 'Visualization Method:', choices = visualizationMethods(), 
                  selected = visualizationMethods()[1], width = "100%")
    }   
  })
  
  output$C_PlotFunction <- renderUI({
    if(is.null(v$data) || is.null(visualizationFunctions())){
      return(NULL)
    }else{
      selectInput('c_PlotFunction', 'Cluster By:', choices = visualizationFunctions(), 
                  selected = visualizationFunctions()[1], width = "100%")
    }   
  })
  
  output$C_markerSelect <- renderUI({
    if(is.null(v$data)){
      return(NULL)
    }else{
      markerNames <- colnames(v$data$expressionData)
      markerNames <- markerNames[order(markerNames)]
      checkboxGroupInput('c_markerSelect', strong('Select Markers:'),
                         markerNames, selected = markerNames, inline = TRUE)
    }   
  })
  
  output$C_clusterSelect <- renderUI({
    if(is.null(v$data) || is.null(v$data$clusterRes) || is.null(input$c_PlotFunction))
      return(NULL)
    if(input$c_PlotFunction %in% c("Sample", "Density","None")){
      return(NULL)
    }else{
      clusterMethod <- input$c_PlotFunction
      clusterIDs <- sort(unique(v$data$clusterRes[[clusterMethod]]))
      selectizeInput('c_clusterSelect', 'Clusters Filter:', 
                     choices = clusterIDs, selected = clusterIDs, 
                     multiple = TRUE, width = "100%")
      # checkboxGroupInput('p_clusterSelect', strong('Select Clusters:'), 
      #                    clusterIDs, selected = clusterIDs, inline = TRUE)
    }   
  })
  
  ## Complex dependencies here: --> (depends on)
  ## C_ScatterPlotInput --> c_PlotMethod + c_clusterSelect 
  ## c_clusterSelect --> c_PlotMethod
  ## carefull checkings are applied to solve concurrency conflicts
  C_ScatterPlotInput <- function(){
    if(is.null(v$data) || is.null(input$c_PlotMethod) || 
       is.null(input$c_PlotFunction) || is.null(input$c_clusterSelect)){
      return(NULL)
    }else if(!all(input$c_clusterSelect %in% v$data$clusterRes[[input$c_PlotFunction]]) &&
             !(input$c_PlotFunction %in% c("Sample", "Density","None"))){
      return(NULL)
    }else{
    # browser()
      
      withProgress(message="Generating Cluster Scatter Plot", value=0, {
        if(input$c_PlotFunction %in% c("Sample", "Density", "None")){
          clusterSelect <- NULL
          clusterColor <- NULL
        }else{
          clusterSelect <- input$c_clusterSelect
          clusterMethod <- input$c_PlotFunction
          if(!is.null(c$clusterCol[[clusterMethod]])){
            clusterColor <- c$clusterCol[[clusterMethod]]
          }else{
            cluster_num <- length(unique(v$data$clusterRes[[clusterMethod]]))
            clusterColor <- rainbow(cluster_num)
          }
        }
        temp_sample_names = lapply(1:length(v$data$sampleNames), function(i){
          v$data$sampleNames[[i]][length(v$data$sampleNames[[i]])]
        })
        if(!all(temp_sample_names %in% input$samples)){
          return(NULL)
        }
        # browser()
        gp <- scatterPlot(obj = v$data,
                          plotMethod = input$c_PlotMethod,
                          plotFunction = input$c_PlotFunction,
                          pointSize = input$C_PointSize,
                          addLabel = input$C_addLabel,
                          labelSize = input$C_LabelSize,
                          sampleLabel = FALSE,
                          FlowSOM_k = input$C_FlowSOM_k, 
                          selectCluster = clusterSelect,
                          selectSamples = input$samples, 
                          facetPlot = input$C_facetPlot,
                          labelRepel = input$C_labelRepel,
                          removeOutlier = TRUE,
                          clusterColor = clusterColor)
        incProgress(1/2)
        plot(gp)
        incProgress(1/2)
      })
    }
  }
  
  output$C_ScatterPlot <- renderPlot({
    C_ScatterPlotInput()
  }, height = 900, width = 950)
  
  output$PDFClusterPlot = downloadHandler(
    filename = function() {
      filename1 <- paste0(input$project_name, "_shinyAPP_Clusterplot_", Sys.Date(), ".pdf")
      i = 0
      while(file.exists(filename1)){
        filename1 <- paste0(input$project_name, "_shinyAPP_Clusterplot_",
                            Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
        i = i + 1;
      }
      filename1
    },
    content = function(file) {
      if(!is.null(v$data)){
        # browser()
        withProgress(message="Downloading Clusterplot PDF files...", value=0, {
          print(getwd())

          pdf(file, 
              width=as.integer(input$tab_w), 
              height=as.integer(input$tab_h))
          C_ScatterPlotInput()
          dev.off()
        })
      }
      
    }
  )
  
  observeEvent(input$PDFClusterPlot, {

  })
  
  
  ##----- change cluster colour -----
  output$C_colourCluster <- renderUI({
    if(is.null(v$data) || is.null(v$data$clusterRes)){
      return(NULL)
    }else{
      clusterMethods <- c(names(v$data$clusterRes)) 
      #clusterMethods <- clusterMethods[!grepl("Subset", clusterMethods)]
      selectInput('c_colourCluster', 'Choose Cluster to Change the Colour :', 
                  choices = clusterMethods, 
                  selected = clusterMethods[1], width = "50%")
    }   
  })
  
  ## currently use 100 as a limit for cluster numbers 
  ## --- TODO: use reactiveValues to automatically retrive cluster numbers --- ## 
  lapply(1:100, function(i) {
    output[[paste0('Cluster_', i, "_col")]] <- renderUI({
      if(is.null(v$data) || is.null(v$data$clusterRes) || is.null(input$c_colourCluster)){
        return(NULL)
      }
      
      clusters <- v$data$clusterRes[[input$c_colourCluster]]
      clusterLabel <- levels(as.factor(clusters))
      if(is.null(c$clusterCol[[input$c_colourCluster]])){
        clusterColor <- rainbow(length(unique(clusters)))
      }else{
        clusterColor <- c$clusterCol[[input$c_colourCluster]]
      }
      
      if (i <= length(clusterLabel)){
        x <- clusterLabel[i]
        colourpicker::colourInput(inputId=paste0('cluster_', i, '_col'), 
                    label=paste0('Cluster ', x," Colour :"), 
                    value = clusterColor[i], showColour = "both", 
                    palette = "square")
      }
    })
  })
  
  ## update cluster color
  observeEvent(input$C_updateClusterColor, {
    if(!is.null(v$data) && !is.null(input$c_colourCluster)){
      clusterMethod <- input$c_colourCluster
      clusterVec<- v$data$clusterRes[[clusterMethod]]
      clusters <- levels(as.factor(clusterVec))
      clusterCols <- NULL
      for (i in 1:length(clusters)){
        clusteri <- clusters[i]
        iCol <- input[[paste0('cluster_', i, '_col')]]
        clusterCols <- c(clusterCols, iCol)
      }
      
      ## update new cluster colours
      c$clusterCol[[clusterMethod]] <- clusterCols
      
      ## jump to C_tab1
      updateTabsetPanel(session, "C_clusterTabs", selected = "C_tab1")
    }
  })
  
  ## revert default cluster colors
  observeEvent(input$C_revertClusterColor, {
    if(!is.null(v$data) && !is.null(input$c_colourCluster)){
      clusterMethod <- input$c_colourCluster
      c$clusterCol[[clusterMethod]] <- NULL
      
      ## jump to C_tab1
      updateTabsetPanel(session, "C_clusterTabs", selected = "C_tab1")
    }
  })
  
  
  ## ------annotate clusters-----
  output$C_labelCluster <- renderUI({
    if(is.null(v$data) || is.null(v$data$clusterRes)){
      return(NULL)
    }else{
      clusterMethods <- c(names(v$data$clusterRes)) 
      #clusterMethods <- clusterMethods[!grepl("Subset", clusterMethods)]
      selectInput('c_labelCluster', 'Choose Cluster Results to Annotate:', 
                  choices = clusterMethods, 
                  selected = clusterMethods[1], width = "50%")
    }   
  })
  
  output$C_labelCluster_name <- renderUI({
    if(is.null(v$data) || is.null(v$data$clusterRes) || is.null(input$c_labelCluster)){
      return(NULL)
    }else{
      textInput("c_labelCluster_name", label = "Type In Your Name for Annotated Cluster", 
                value = paste0("Annotated_", input$c_labelCluster), width = "50%")
    }
  })
  
  
  ## currently use 100 as a limit for cluster numbers 
  ## --- TODO: use reactiveValues to automatically retrive cluster numbers --- ## 
  lapply(1:100, function(i) {
    output[[paste0('Cluster', i)]] <- renderUI({
      if(is.null(v$data) || is.null(v$data$clusterRes) || is.null(input$c_labelCluster)){
        return(NULL)
      }
      
      # create new item in RData object
      clusters <- sort(unique(v$data$clusterRes[[input$c_labelCluster]]))
      if (i <= length(clusters)){
        x <- clusters[i]
        textInput(paste0('cluster', i), paste0('Cluster ', x," :"), 
                  value = "", width = "30%", placeholder = "Type in the cell type")
      }
    })
  })
  
  ## update cluster labels
  observeEvent(input$updatelabel, {
    if(!is.null(v$data) && !is.null(input$c_labelCluster) && !is.null(input$c_labelCluster_name)){
      obj <- v$data
      clusterMethod <- input$c_labelCluster
      clusterVec<- obj$clusterRes[[clusterMethod]]
      clusterLabels <- clusterVec
      clusters <- sort(unique(clusterVec))
      
      for (i in 1:length(clusters)){
        clusteri <- clusters[i]
        ilabel <- input[[paste0('cluster', i)]]
        # if(ilabel == ""){
        #   clusterLabels[clusterLabels==clusteri] <- "Unknown"
        # }else{
        #   clusterLabels[clusterLabels==clusteri] <- ilabel
        # }
        if(ilabel != ""){
          clusterLabels[clusterLabels==clusteri] <- ilabel
        }
      }
      
      ## update new cluster results
      labelName <- input$c_labelCluster_name
      obj$clusterRes[[labelName]] <- clusterLabels
      
      ## update the project name
      obj$projectName <- paste0(obj$projectName, "_annotated_", clusterMethod)
      
      v$data <- obj
      
      ## jump to C_tab1
      updateTabsetPanel(session, "C_clusterTabs", selected = "C_tab1")
    }
  })
  
  
  
  ##-----RUN flowSOM-----
  ## result object which will be updated by C_runFlowSOM
  observeEvent(input$C_runFlowSOM, {
    if(!is.null(v$data) && !is.null(input$c_markerSelect)){
      obj <- v$data
      withProgress(message=paste0('Running FlowSOM using k=', input$C_FlowSOM_k), value=0, {
        FlowSOM_cluster <- cytof_cluster(xdata = obj$expressionData[ ,input$c_markerSelect],
                                         method = "FlowSOM",
                                         FlowSOM_k = input$C_FlowSOM_k)
        incProgress(1/2)
        ## update FlowSOM cluster results
        obj$clusterRes[["FlowSOM"]] <- FlowSOM_cluster
        ## update the project name
        obj$projectName <- paste0(obj$projectName, "_add_FlowSOM")
        v$data <- obj
        incProgress(1/2)
      })
      
      ## jump to C_tab1
      updateTabsetPanel(session, "C_clusterTabs", selected = "C_tab1")
    }
  })
  
  
  ##------------------------------Marker Panel-------------------------------
  
  ##-----heat map plot-----
  output$M_plotCluster <- renderUI({
    if(is.null(v$data) || is.null(clusterMethods())){
      return(NULL)
    }else{
      selectInput('m_plotCluster', 'Cluster Method:', choices = clusterMethods(), 
                  selected = clusterMethods()[1], width = "100%")
    }   
  })
  
  output$M_heatmapmarkerSelect <- renderUI({
    if(is.null(v$data)){
      return(NULL)
    }else{
      sorted_markerNames <- colnames(v$data$expressionData)
      markerNames <- sorted_markerNames[order(sorted_markerNames)]
      initNum <- ifelse(length(markerNames) >=4, 4, 1)
      selectizeInput('m_heatmapmarkerSelect', 'Select Markers:', 
                     choices = markerNames, selected = markerNames[1:initNum], 
                     multiple = TRUE, width = "100%")
    }   
  })
  
  observeEvent(input$M_heatmapSelectAll, {
    raw_markers <- colnames(v$data$expressionData)
    markers <- raw_markers[order(raw_markers)]
    updateSelectizeInput(session, "m_heatmapmarkerSelect", selected = markers)
  })
  
  M_heatmapPlotInput <- reactive({
    if(is.null(v$data) || is.null(input$m_plotCluster) || is.null(input$m_heatmapmarkerSelect))
      return(NULL)
    heatMap(data = v$data, 
            clusterMethod = input$m_plotCluster, 
            type = input$M_plotMethod, 
            dendrogram = input$M_heatmap_dendrogram,
            colPalette = input$M_heatmap_colorPalette,
            selectSamples = input$samples,
            selectMarkers = input$m_heatmapmarkerSelect,
            cex_row_label= input$M_rowLabelSize, 
            cex_col_label= input$M_colLabelSize, 
            scaleMethod = input$M_scaleMethod)
  })
  
  output$M_heatmapPlot <- renderPlot({
    M_heatmapPlotInput()
  }, height = 900, width = 950)
  
  output$PDFHeatmap = downloadHandler(
    filename = function() {
      filename1 <- paste0(input$project_name, "_shinyAPP_Marker_Heatmap_", Sys.Date(), ".pdf")
      i = 0
      while(file.exists(filename1)){
        filename1 <- paste0(input$project_name,
                            "_shinyAPP_Marker_Heatmap_",
                            Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
        i = i + 1;
      }
      filename1
    },
    content = function(file) {
      if(!is.null(v$data)){
        withProgress(message="Downloading Marker Heatmap PDF files...", value=0, {
          pdf(file,
              width=as.integer(input$tab_w),
              height=as.integer(input$tab_h))
          heatMap(data = v$data, 
                  clusterMethod = input$m_plotCluster, 
                  type = input$M_plotMethod, 
                  dendrogram = input$M_heatmap_dendrogram,
                  colPalette = input$M_heatmap_colorPalette,
                  selectSamples = input$samples,
                  selectMarkers = input$m_heatmapmarkerSelect,
                  cex_row_label= input$M_rowLabelSize, 
                  cex_col_label= input$M_colLabelSize, 
                  scaleMethod = input$M_scaleMethod)
          dev.off()
        })
      }
    }
  )
  

  
  session$onSessionEnded(function(){
    # file.remove("cytofkit_shinyAPP_marker_heatmap.pdf")
  })
  
  ##-----level plot-----
  output$M_PlotMethod <- renderUI({
    if(is.null(v$data) || is.null(visualizationMethods())){
      return(NULL)
    }else{
      selectInput('m_PlotMethod', 'Visualization Method:', choices = visualizationMethods(), 
                  selected = visualizationMethods()[1], width = "100%")
    }   
  })
  
  output$M_PlotMarker <- renderUI({
    if(is.null(v$data)){
      return(NULL)
    }else{
      sorted_markers <- colnames(v$data$expressionData)
      sorted_markers <- sorted_markers[order(sorted_markers)]
      #markers <- c(sorted_markers, "All Markers", "All Markers(scaled)")
      selectizeInput('m_PlotMarker', 'Plot Marker:', choices = sorted_markers, 
                     selected = sorted_markers[1], multiple = TRUE, width = "100%")
    }   
  })
  
  observeEvent(input$M_chooseAllMarker, {
    raw_markers <- colnames(v$data$expressionData)
    markers <- raw_markers[order(raw_markers)]
    updateSelectizeInput(session, "m_PlotMarker", selected = markers)
  })
  
  M_markerExpressionPlotInput <- function(){
    if(is.null(v$data) || is.null(input$m_PlotMethod) || is.null(isolate(input$m_PlotMarker))){
      return(NULL)
    }else{
      withProgress(message="Generating Marker Expression Plot", value=0, {
        gp <- scatterPlot(obj = v$data,
                          plotMethod = input$m_PlotMethod,
                          plotFunction = isolate(input$m_PlotMarker),
                          pointSize = input$M_PointSize,
                          alpha = input$M_Alpha,
                          addLabel = FALSE,
                          labelSize = input$S_LabelSize,
                          sampleLabel = FALSE,
                          FlowSOM_k = input$C_FlowSOM_k, 
                          selectSamples = input$samples, 
                          facetPlot = FALSE,
                          colorPalette = input$M_colorPalette,
                          labelRepel = FALSE,
                          removeOutlier = TRUE,
                          globalScale = ifelse(input$M_ScaleOptions == "Global", TRUE, FALSE),
                          centerScale = ifelse(input$M_scaledData == "Centered", TRUE, FALSE))
        incProgress(1/2)
        plot(gp)
        incProgress(1/2)
      })
    }
  }
  
  output$M_markerExpressionPlot <- renderPlot({
    M_markerExpressionPlotInput()
  }, height = 900, width = 950)
  
  observeEvent({
    input$M_updateExPlot
    input$m_PlotMethod
    input$M_PointSize
    input$S_LabelSize
    input$M_colorPalette
    input$M_ScaleOptions
    input$M_scaledData
  }, {
    output$M_markerExpressionPlot <- renderPlot({
      M_markerExpressionPlotInput()
    }, height = 900, width = 950)
  })
  
  
  output$PDFExpPlot = downloadHandler(
    filename = function() {
      paste0(input$project_name, '_expression_plots.zip')
      filename1 <- paste0(input$project_name, "_shinyAPP_Marker_Expression_Plot_", Sys.Date(), ".pdf")
      i = 0
      while(file.exists(filename1)){
        filename1 <- paste0(input$project_name,
                            "_shinyAPP_Marker_Expression_Plot_",
                            Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
        i = i + 1;
      }
      filename1
    },
    content = function(file) {
      if(!is.null(v$data)){
        withProgress(message="Downloading Marker Expression Plot PDF files...", value=0, {
          print(getwd())

          pdf(file, 
              width=as.integer(input$tab_w), 
              height=as.integer(input$tab_h))
          M_markerExpressionPlotInput()
          dev.off()
        })
      }
    }
  )
  
  ##-----histogram plot-----
  output$M_stackFactor <- renderUI({
    if(is.null(v$data)){
      return(NULL)
    }else{
      stackFactorChoice <- c(names(v$data$clusterRes), "sample") 
      selectInput('m_stackFactor', 'Stack Factor:', choices = stackFactorChoice, 
                  selected = stackFactorChoice[1], width = "100%")
    }   
  })
  
  output$M_markerSelect <- renderUI({
    if(is.null(v$data)){
      return(NULL)
    }else{
      sorted_markerNames <- colnames(v$data$expressionData)
      markerNames <- sorted_markerNames[order(sorted_markerNames)]
      initNum <- ifelse(length(markerNames) >=4, 4, 1)
      selectizeInput('m_markerSelect', 'Select Markers:', 
                     choices = markerNames, selected = markerNames[1:initNum], 
                     multiple = TRUE, width = "100%")
    }   
  })
  
  observeEvent(input$M_histSelectAll, {
    raw_markers <- colnames(v$data$expressionData)
    markers <- raw_markers[order(raw_markers)]
    updateSelectizeInput(session, "m_markerSelect", selected = markers)
  })
  
  M_stackDensityPlotInput <- function(){
    m_markerSelect <- isolate(input$m_markerSelect)
    if(is.null(v$data) || is.null(input$m_stackFactor) || is.null(m_markerSelect)){
      return(NULL)
    }else{
      withProgress(message="Generating Stack Density Plot", value=0, {
        data <- data.frame(v$data$expressionData, check.names = FALSE)
        samples <- as.character(v$sampleInfo$cellSample)
        mySamples <- samples %in% input$samples
        sfactors <- data.frame(do.call(cbind, v$data$clusterRes), 
                               sample = samples, 
                               stringsAsFactors = FALSE, 
                               check.names = FALSE)
        data <- data[mySamples, ,drop=FALSE]
        stackFactor <- sfactors[mySamples, input$m_stackFactor]
        
        if(input$m_stackFactor == "sample"){
          stackFactorColours <- NULL
        }else{
          clusterMethod <- input$m_stackFactor
          clusterVec <- v$data$clusterRes[[clusterMethod]]
          cluster_num <- length(unique(clusterVec))
          selectColors <- match(levels(as.factor(stackFactor)), levels(as.factor(clusterVec)))
          if(!is.null(c$clusterCol[[clusterMethod]])){
            stackFactorColours <- c$clusterCol[[clusterMethod]][selectColors]
          }else{
            stackFactorColours <- rainbow(cluster_num)[selectColors]
          }
        }
        
        incProgress(1/3)
        gp <- stackDenistyPlot(data = data, 
                               densityCols=m_markerSelect, 
                               stackFactor = stackFactor,
                               kernel = "gaussian",
                               bw = "nrd0", 
                               adjust = 1,
                               stackRotation = 0, 
                               stackSeperation = "auto",
                               x_text_size = input$M_xlab_size, 
                               strip_text_size = input$M_markerTextSize,
                               legend_text_size = input$M_legendTextSize, 
                               legendRow = input$M_legendRow,
                               legend_title = input$m_stackFactor,
                               stackFactorColours = stackFactorColours)
        incProgress(1/3)
        plot(gp)
        incProgress(1/3)
      })
    }
  }
  
  observeEvent(input$M_updateDensityPlot, {
    output$M_stackDensityPlot <- renderPlot({
      M_stackDensityPlotInput()
    }, height = 900, width = 950)
  })
  
  output$PDFHistogram = downloadHandler(
    filename = function() {
      # browser()
      filename1 <- paste0(input$project_name, "_shinyAPP_Stack_Density_Plot_", Sys.Date(), ".pdf")
      i = 0
      while(file.exists(filename1)){
        filename1 <- paste0(input$project_name,
                            "_shinyAPP_Stack_Density_Plot_",
                            Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
        i = i + 1;
      }
      filename1
    },
    content = function(file) {
      if(!is.null(v$data)){
        withProgress(message="Downloading Stack Density Plot PDF files...", value=0, {
          print(getwd())
          pdf(file, 
              width=as.integer(input$tab_w), 
              height=as.integer(input$tab_h))
          M_stackDensityPlotInput()
          dev.off()
        })
      }
    }
  )
  
  observeEvent(input$PDFHistogram, {

  })
  
  
  ##----- update marker names -----
  
  ## currently use 100 as a limit for marker number
  ## --- TODO: use reactiveValues to automatically retrive marker numbers --- ## 
  lapply(1:100, function(i) {
    output[[paste0('Marker_', i, "_name")]] <- renderUI({
      if(is.null(v$data)){
        return(NULL)
      }
      sorted_markerNames <- colnames(v$data$expressionData)
      markerNames <- sorted_markerNames[order(sorted_markerNames)]
      
      if (i <= length(markerNames)){
        markeri <- markerNames[i]
        textInput(inputId = paste0('marker_', i, "_name"), 
                  label = markeri, value = markeri, width = "30%", 
                  placeholder = "Type in your new name for this marker")
      }
    })
  })
  
  
  ## update cluster labels
  observeEvent(input$C_updateMarkerNames, {
    if(!is.null(v$data)){
      sorted_markerNames <- colnames(v$data$expressionData)
      markerNames <- sorted_markerNames[order(sorted_markerNames)]
      newMarkerNames <- NULL
      for (i in 1:length(markerNames)){
        iName <- input[[paste0('marker_', i, '_name')]]
        newMarkerNames <- c(newMarkerNames, iName)
        
      }
      ## update new cluster colours
      mark_pos = which(colnames(v$data$expressionData) %in% markerNames)
      v$data$expressionData[, mark_pos] = v$data$expressionData[, markerNames]
      colnames(v$data$expressionData)[mark_pos] <- newMarkerNames
      ## jump to C_tab1
      updateTabsetPanel(session, "M_markerTabs", selected = "M_tab1")
    }
  })
  
  
  ##------------------------------Sample Panel-------------------------------
  
  ##-----cell percentage heatmap-----
  output$S_plotCluster <- renderUI({
    if(is.null(v$data) || is.null(clusterMethods())){
      return(NULL)
    }else{
      selectInput('s_plotCluster', 'Cluster Method:', choices = clusterMethods(), 
                  selected = clusterMethods()[1], width = "100%")
    }   
  })
  
  S_heatmapPlotInput <- reactive({
    if(is.null(v$data) || is.null(clusterMethods()) || is.null(input$s_plotCluster))
      return(NULL)
    
    heatMap(data = v$data, 
            clusterMethod = input$s_plotCluster, 
            type = input$S_plotMethod, 
            dendrogram = input$S_heatmap_dendrogram,
            colPalette = input$S_heatmap_colorPalette,
            selectSamples = input$samples,
            cex_row_label= input$S_rowLabelSize, 
            cex_col_label= input$S_colLabelSize, 
            scaleMethod = input$S_scaleMethod)
    
  })
  
  output$S_heatmapPlot <- renderPlot({
    S_heatmapPlotInput()
  }, height = 900, width = 950)
  
  output$PDFSamHeat = downloadHandler(
    filename = function() {
      filename1 <- paste0(input$project_name, "_shinyAPP_Sample_Heatmap_", Sys.Date(), ".pdf")
      i = 0
      while(file.exists(filename1)){
        filename1 <- paste0(input$project_name,
                            "_shinyAPP_Sample_Heatmap_",
                            Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
        i = i + 1;
      }
      filename1
    },
    content = function(file) {
      if(!is.null(v$data)){
        withProgress(message="Downloading Sample Heatmap PDF files...", value=0, {
          pdf(file,
              width=as.integer(input$tab_w),
              height=as.integer(input$tab_h))
          heatMap(data = v$data, 
                  clusterMethod = input$s_plotCluster, 
                  type = input$S_plotMethod, 
                  dendrogram = input$S_heatmap_dendrogram,
                  colPalette = input$S_heatmap_colorPalette,
                  selectSamples = input$samples,
                  cex_row_label= input$S_rowLabelSize, 
                  cex_col_label= input$S_colLabelSize, 
                  scaleMethod = input$S_scaleMethod)
          dev.off()
        })
      }
    }
  )
  
  observeEvent(input$PDFSamHeat, {
    
  })
  
  
  ##-----cell percentage line chart-----
  output$S_clusterMethod2 <- renderUI({
    if(is.null(v$data) || is.null(clusterMethods())){
      return(NULL)
    }else{
      selectInput('s_clusterMethod2', 'Cluster Method:', choices = clusterMethods(), 
                  selected = clusterMethods()[1], width = "100%")
    }   
  })
  
  output$S_clusterFilter <- renderUI({
    if(is.null(v$data) || is.null(clusterMethods()) || is.null(input$s_clusterMethod2)){
      return(NULL)
    }else{
      clusterIDs <- sort(unique(v$data$clusterRes[[input$s_clusterMethod2]]))
      selectizeInput('s_clusterFilter', 'Filter Clusters:', 
                     choices = clusterIDs, selected = clusterIDs, 
                     multiple = TRUE, width = "100%")
    }   
  })
  
  S_rateChangePlotInput <- function(){
    if(is.null(v$data) || is.null(clusterMethods()) || is.null(input$s_clusterMethod2) || is.null(input$s_clusterFilter))
      return(NULL)
    withProgress(message="Generating Rate Change Plot", value=0, {
      ## percentage stat
      data <- data.frame(sample = v$sampleInfo$cellSample,
                         cluster = as.factor(v$data$clusterRes[[input$s_clusterMethod2]]),
                         counts = 1)
      statData1 <- aggregate(counts ~ ., data = data, sum)
      statData2 <- aggregate(counts ~ sample, data = data, sum)
      statData <- merge(statData1, statData2, by="sample", suffixes = c("InAll","InSample"))
      statData$percentageInSample <- statData$countsInAll/statData$countsInSample
      incProgress(1/3)
      ## filter clusters
      usedClusters <- input$s_clusterFilter
      clusterCheck <- as.character(statData$cluster) %in% usedClusters
      statData <- statData[clusterCheck, ,drop=FALSE]
      incProgress(1/3)
      gp <- ggplot(data = statData, aes_string(x="sample", 
                                               y="percentageInSample", 
                                               color = "cluster",
                                               group = "cluster")) + 
        geom_point(size = 2) + geom_line(size = 1.5) + 
        xlab("Cell Group") + ylab("Percentage of Cells in Group") + theme_bw() + 
        theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
      incProgress(1/3)
      plot(gp)
    })
  }
  
  output$S_rateChangePlot <- renderPlot({
    S_rateChangePlotInput()
  }, height = 500, width = 950)
  
  output$PDFrateChange = downloadHandler(
    filename = function() {
      filename1 <- paste0(input$project_name, "_shinyAPP_Rate_Change_Plot_", Sys.Date(), ".pdf")
      i = 0
      while(file.exists(filename1)){
        filename1 <- paste0(input$project_name,
                            "_shinyAPP_Rate_Change_Plot_",
                            Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
        i = i + 1;
      }
      filename1
    },
    content = function(file) {
      if(!is.null(v$data)){
        withProgress(message="Downloading Rate Change Plot PDF files...", value=0, {
          pdf(file, 
              width=as.integer(input$tab_w), 
              height=as.integer(input$tab_h))
          S_rateChangePlotInput()
          dev.off()
        })
      }
    }
  )
  
  observeEvent(input$PDFrateChange, {
    
  })
  
  
  # output$S_clusterTable <- renderTable({
  #     if(is.null(v$data) || is.null(clusterMethods()) || is.null(input$s_clusterMethod2)){
  #         return(NULL)
  #     }else{
  #         data <- data.frame(sample = v$sampleInfo$cellSample,
  #                            cluster = as.factor(v$data$clusterRes[[input$s_clusterMethod2]]),
  #                            counts = 1)
  #         
  #         statData1 <- aggregate(counts ~ ., data = data, sum)
  #         statData2 <- aggregate(counts ~ sample, data = data, sum)
  #         statData <- merge(statData1, statData2, by="sample", suffixes = c("InAll","InSample"))
  #         if(is.numeric(statData$cluster)) statData$cluster <- as.integer(statData$cluster)
  #         statData$counts <- as.integer(statData$countsInAll)
  #         statData$percentageInAll <- round(statData$countsInAll/nrow(data), 4)
  #         statData$percentageInSample <- round(statData$countsInAll/statData$countsInSample, 2)
  #         statData[, c("sample", "cluster", "counts", "percentageInSample", "percentageInAll")]
  #     }   
  # }) 
  
  
  ##-----Regroup samples-----
  output$S_groupSamples <- renderUI({
    if(is.null(v$data) || is.null(v$data$clusterRes)){
      return(NULL)
    }else{
      clusterMethods <- c(names(v$data$clusterRes)) 
      #clusterMethods <- clusterMethods[!grepl("Subset", clusterMethods)]
      selectInput('c_labelCluster', 'Choose Cluster Results to Annotate:', 
                  choices = clusterMethods, 
                  selected = clusterMethods[1], width = "30%")
    }   
  })
  
  ## currently use 100 as a limit for sample numbers 
  ## --- TODO: use reactiveValues to automatically retrive sample numbers --- ## 
  lapply(1:100, function(i) {
    output[[paste0('S_sample', i)]] <- renderUI({
      if(is.null(v$data) || is.null(v$sampleInfo)){
        return(NULL)
      }
      
      uniqueSampleNames <- sort(unique(v$sampleInfo$cellSample))
      if (i <= length(uniqueSampleNames)){
        x <- uniqueSampleNames[i]
        textInput(paste0('Sample', i), paste0(x," :"), 
                  value = "", width = "40%", 
                  placeholder = "Type in the group name for this sample")
      }
    })
  })
  
  
  rename_sample = function(new_sample_name){
    # browser()
    v$sampleInfo$originalCellSample <- v$sampleInfo$cellSample
    uniqueSampleNames <- sort(unique(v$sampleInfo$originalCellSample))

    temp_names = sapply(1:length(uniqueSampleNames), function(i){
      sample_name_length = length(v$data$sampleNames[[i]])
      v$data$sampleNames[[i]][sample_name_length]
    })
    v$data$sampleNames = v$data$sampleNames[order(temp_names)]

    sampleGroupNames <- NULL
    for(i in 1:length(uniqueSampleNames)){
      sample_name_length = length(v$data$sampleNames[[i]])
      if (new_sample_name[i] != "") {
        sampleGroupNames <- c(sampleGroupNames, new_sample_name[i])
        v$data$sampleNames[[i]] <- c(v$data$sampleNames[[i]][1], new_sample_name[i])
      } else {
        sampleGroupNames <- c(sampleGroupNames, v$data$sampleNames[[i]][sample_name_length])
        # v$data$sampleNames[[i]] <- c(v$data$sampleNames[[i]][sample_name_length], v$data$sampleNames[[i]][sample_name_length])
      }
    }

    groupNameLevels <- strsplit(input$sampleGroupLevels, ";", fixed = TRUE)[[1]]

    if(groupNameLevels != "" && all(sampleGroupNames != "")
       && length(groupNameLevels) == length(unique(sampleGroupNames))
       && all(as.character(groupNameLevels) %in% sampleGroupNames)){
      sampleMatchID <- match(v$sampleInfo$originalCellSample, uniqueSampleNames)
      v$sampleInfo$cellSample <- factor(sampleGroupNames[sampleMatchID],
                                        levels = groupNameLevels)
    }else{
      sampleGroupNames[sampleGroupNames == ""] <- uniqueSampleNames[sampleGroupNames == ""]
      sampleMatchID <- match(v$sampleInfo$originalCellSample, uniqueSampleNames)
      v$sampleInfo$cellSample <- factor(sampleGroupNames[sampleMatchID])
    }

    cellID_number <- do.call(base::c, regmatches(v$sampleInfo$cellID,
                                                 gregexpr("_[0-9]*$", v$sampleInfo$cellID, perl=TRUE)))

    ## update reactive object v$sampleInfo
    ## newCellID = "sampleGroup" + "_cellID" + "globalID" to avoid duplicates
    v$sampleInfo$newCellID <- paste0(as.character(v$sampleInfo$cellSample),
                                     "_",
                                     1:length(cellID_number))


    ## update reactive object v$data
    expressionData <- v$data$expressionData
    row.names(expressionData) <- v$sampleInfo$newCellID
    v$data$expressionData <- expressionData

    ## update the project name
    v$data$projectName <- paste0(v$data$projectName, "_grouped_samples")

    ## update v$data$progressionRes
    if(!is.null(v$data$progressionRes)){
      sampleExpressData <- v$data$progressionRes$sampleData
      row.names(sampleExpressData) <- v$sampleInfo$newCellID[match(row.names(sampleExpressData),
                                                                   v$sampleInfo$cellID)]
      v$data$progressionRes$sampleData <- sampleExpressData
    }

    ## jump to S_tab1
    # updateTabsetPanel(session, "S_sampleTabs", selected = "S_tab1")
  }

  
  reset_sample = function(){
    # browser()
    v$sampleInfo <- v$ori_sampleInfo
    v$data = v$ori_data
    v$sample_selected_index = 1:length(unique(as.character(v$sampleInfo$cellSample)))
    
  }
  
  ## update sample groups
  observeEvent(input$updateSampleGroups, {
    if(!is.null(v$data) && !is.null(v$sampleInfo)){
      # browser()
      v$sampleInfo$originalCellSample <- v$sampleInfo$cellSample
      uniqueSampleNames <- sort(unique(v$sampleInfo$originalCellSample))
      new_sample_names <- NULL
      for(i in 1:length(uniqueSampleNames)){
        new_sample_names <- c(new_sample_names, input[[paste0("Sample", i)]])
      }
      rename_sample(new_sample_names)
      
      ## jump to S_tab1
      updateTabsetPanel(session, "S_sampleTabs", selected = "S_tab1")
    }
  })
  
  ## revert old sample names
  observeEvent(input$revertSampleNames, {
    if(!is.null(v$data) && !is.null(v$sampleInfo)){
      if(!is.null(v$sampleInfo$originalCellSample)){
        v$sampleInfo$cellSample <- v$sampleInfo$originalCellSample
        v$sampleInfo$originalCellSample <- NULL
        
        ## update reactive object v$data
        expressionData <- v$data$expressionData
        row.names(expressionData) <- v$sampleInfo$cellID
        v$data$expressionData <- expressionData
        
        ## update the project name
        v$data$projectName <- sub("_grouped_samples", "", v$data$projectName)
        
        ## update reactive object v$sampleInfo
        if(!is.null(v$data$progressionRes)){
          sampleExpressData <- v$data$progressionRes$sampleData
          row.names(sampleExpressData) <- v$sampleInfo$cellID[match(row.names(sampleExpressData),
                                                                    v$sampleInfo$newCellID)]
          v$data$progressionRes$sampleData <- sampleExpressData
        }
      }
      ## jump to S_tab1
      updateTabsetPanel(session, "S_sampleTabs", selected = "S_tab1")
    }
  })
  
  
  
  ##---------------------------Progression Panel------------------------------
  
  ##-----subset relationship plot-----
  
  output$P_xlab <- renderUI({
    if(is.null(v$data) || is.null(progressionLabs())){
      return(NULL)
    }else{
      selectInput('p_xlab', 'Plot X:', choices = progressionLabs(), 
                  selected = progressionLabs()[1], width = "100%")
    }   
  })
  
  output$P_ylab <- renderUI({
    if(is.null(v$data) || is.null(progressionLabs())){
      return(NULL)
    }else{
      selectInput('p_ylab', 'Plot Y:', choices = progressionLabs(), 
                  selected = progressionLabs()[2], width = "100%")
    }   
  })
  
  P_scatterPlotInput <- function(){
    if(is.null(v$data) || is.null(v$data$progressionRes) || is.null(input$p_xlab) || is.null(input$p_ylab)){
      return(NULL)
    }else{
      withProgress(message="Generating Progression Scatter Plot", value=0, {
        obj <- v$data$progressionRes
        data <- data.frame(obj$progressionData, 
                           cluster = obj$sampleCluster,
                           sample = sub("_[0-9]*$", "", row.names(obj$sampleData)))
        incProgress(1/3)
        data <- data[data$sample %in% input$samples, ,drop=FALSE]
        
        clusterMethod <- p$progressionCluster
        clusterVec <- v$data$clusterRes[[clusterMethod]]
        cluster_num <- length(unique(clusterVec))
        selectColors <- match(levels(as.factor(data$cluster)), levels(as.factor(clusterVec)))
        
        if(!is.null(c$clusterCol[[clusterMethod]])){
          clusterColor <- c$clusterCol[[clusterMethod]][selectColors]
        }else{
          clusterColor <- rainbow(cluster_num)[selectColors]
        }
        
        gp <- cytof_clusterPlot(data = data, 
                                xlab = input$p_xlab, 
                                ylab = input$p_ylab, 
                                cluster = "cluster", 
                                sample = "sample",
                                title = "Subset Relationship", 
                                type = ifelse(input$P_facetPlot, 2, 1),
                                point_size = input$P_PointSize, 
                                addLabel = input$P_addLabel, 
                                labelSize = input$P_LabelSize, 
                                sampleLabel = FALSE, 
                                labelRepel = input$P_labelRepel,
                                fixCoord = FALSE,
                                clusterColor = clusterColor)
        incProgress(1/3)
        plot(gp)
        incProgress(1/3)
      })
    }
  }
  
  output$P_scatterPlot <- renderPlot({
    P_scatterPlotInput()
  }, height = 900, width = 950)
  
  output$PDFScatter = downloadHandler(
    filename = function() {
      filename1 <- paste0(input$project_name, "_shinyAPP_Scatterplot_", Sys.Date(), ".pdf")
      i = 0
      while(file.exists(filename1)){
        filename1 <- paste0(input$project_name,
                            "_shinyAPP_Scatterplot_",
                            Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
        i = i + 1;
      }
      filename1
    },
    content = function(file) {
      if(!is.null(v$data)){
        withProgress(message="Downloading Progression Scatterplot PDF files...", value=0, {
          print(getwd())
          pdf(file, 
              width=as.integer(input$tab_w), 
              height=as.integer(input$tab_h))
          P_scatterPlotInput()
          dev.off()
        })
      }
    }
  )
  
  
  
  ##-----marker expression profile-----
  
  output$P_orderBy <- renderUI({
    if(is.null(v$data) || is.null(progressionLabs())){
      return(NULL)
    }else{
      selectInput('p_orderBy', 'Cell Order By:', choices = progressionLabs(), 
                  selected = progressionLabs()[1], width = "100%")
    }   
  })
  
  output$P_markerSelect <- renderUI({
    if(is.null(v$data) || is.null(v$data$progressionRes)){
      return(NULL)
    }else{
      sorted_markerNames <- colnames(v$data$progressionRes$sampleData)  
      markerNames <- sorted_markerNames[order(sorted_markerNames)]
      initNum <- ifelse(length(markerNames) >=4, 4, 1)
      selectizeInput('p_markerSelect', 'Select Markers:', 
                     choices = markerNames, selected = markerNames[1:initNum], 
                     multiple = TRUE, width = "100%")
      # checkboxGroupInput('p_markerSelect', strong('Select Markers:'), 
      #                    markerNames, selected = markerNames, inline = TRUE)
    }   
  })
  
  output$P_clusterSelect <- renderUI({
    if(is.null(v$data) || is.null(v$data$progressionRes)){
      return(NULL)
    }else{
      clusterIDs <- sort(unique(v$data$progressionRes$sampleCluster))
      selectizeInput('p_clusterSelect', 'Select Clusters:', 
                     choices = clusterIDs, selected = clusterIDs, 
                     multiple = TRUE, width = "100%")
      # checkboxGroupInput('p_clusterSelect', strong('Select Clusters:'), 
      #                    clusterIDs, selected = clusterIDs, inline = TRUE)
    }   
  })
  
  P_markerPlotInput <- function(){
    p_markerSelect <- isolate(input$p_markerSelect)
    p_clusterSelect <- isolate(input$p_clusterSelect)
    if(is.null(v$data) || is.null(v$data$progressionRes) || is.null(p_markerSelect) || is.null(p_clusterSelect) || is.null(input$p_orderBy))
      return(NULL)
    
    withProgress(message="Generating Marker Expression Profile", value=0, {
      data <- data.frame(v$data$progressionRes$sampleData,
                         cluster = v$data$progressionRes$sampleCluster, 
                         v$data$progressionRes$progressionData,
                         check.names = FALSE)
      
      sampleNames <- sub("_[0-9]*$", "", row.names(v$data$progressionRes$sampleData))
      data <- data[sampleNames %in% input$samples, ,drop=FALSE]
      incProgress(1/3)
      if(input$P_combineTrends){
        pp <- cytof_expressionTrends(data, 
                                     markers = p_markerSelect, 
                                     clusters = p_clusterSelect, 
                                     orderCol = input$p_orderBy, 
                                     clusterCol = "cluster", 
                                     reverseOrder = input$P_reverseOrder,
                                     addClusterLabel = input$P_addLabel2,
                                     clusterLabelSize = input$P_LabelSize2,
                                     segmentSize = 0.5,
                                     min_expr = NULL) 
      }else{
        pp <- cytof_progressionPlot(data, 
                                    markers = p_markerSelect, 
                                    clusters = p_clusterSelect, 
                                    orderCol = input$p_orderBy, 
                                    clusterCol = "cluster", 
                                    reverseOrder = input$P_reverseOrder,
                                    addClusterLabel = input$P_addLabel2,
                                    clusterLabelSize = input$P_LabelSize2,
                                    segmentSize = 0.5,
                                    min_expr = NULL) 
      }
      incProgress(1/3)
      plot(pp)
      incProgress(1/3)
    })
  }
  
  observeEvent(input$P_updateRegressionPlot, {
    output$P_markerPlot <- renderPlot({
      P_markerPlotInput()
    }, height = 900, width = 950)  
  })
  
  output$PDFmarkerPlot = downloadHandler(
    filename = function() {
      filename1 <- paste0(input$project_name, "_shinyAPP_Marker_Plot_", Sys.Date(), ".pdf")
      i = 0
      while(file.exists(filename1)){
        filename1 <- paste0(project_name,
                            "_shinyAPP_Marker_Plot_",
                            Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
        i = i + 1;
      }
      filename1
    },
    content = function(file) {
      if(!is.null(v$data)){
        withProgress(message="Downloading Marker Plot PDF files...", value=0, {
          print(getwd())
          pdf(file, 
              width=as.integer(input$tab_w), 
              height=as.integer(input$tab_h))
          P_markerPlotInput()
          dev.off()
        })
      }
    }
  )
  
  
  
  ##-----Run Diffusionmap-----
  
  output$P_clusterMethod <- renderUI({
    if(is.null(v$data) || is.null(clusterMethods())){
      return(NULL)
    }else{
      selectInput('p_clusterMethod', 'Cluster Method:', choices = clusterMethods(), 
                  selected = clusterMethods()[1], width = "100%")
    }   
  })
  
  output$P_clusterTable <- renderTable({
    if(is.null(v$data) || is.null(clusterMethods())){
      return(NULL)
    }else{
      clusterTable <- t(as.matrix(table(v$data$clusterRes[[input$p_clusterMethod]])))
      out <- as.data.frame(clusterTable, row.names = "Cell Counts")
      colnames(out) <- paste("Cluster", colnames(out))
      out
    }   
  })
  
  output$P_clusterFilter <- renderUI({
    if(is.null(v$data) || is.null(clusterMethods())){
      return(NULL)
    }else{
      obj <- v$data
      clusterIDs <- sort(unique(obj$clusterRes[[input$p_clusterMethod]]))
      selectizeInput('p_clusterFilter', 'Filter Clusters:', 
                     choices = clusterIDs, selected = clusterIDs, 
                     multiple = TRUE, width = "100%")
    }   
  })
  
  
  ## result object which will be updated by P_runDiffusionmap
  observeEvent(input$P_runDiffusionmap, {
    
    if(!is.null(v$data)){
      obj <- v$data
      usedClusters <- input$p_clusterFilter
      clusterCheck <- obj$clusterRes[[input$p_clusterMethod]] %in% usedClusters
      mdata <- obj$expressionData[clusterCheck, ]
      mcluster <- obj$clusterRes[[input$p_clusterMethod]][clusterCheck]
      withProgress(message="Running Diffusionmap", value=0, {
        diffmapRes <- cytof_progression(data = mdata, 
                                        cluster = mcluster, 
                                        method = "diffusionmap", 
                                        distMethod = input$P_distMethod,
                                        out_dim = input$P_outDim,
                                        clusterSampleMethod = input$P_sampleMethod,
                                        clusterSampleSize = input$P_clusterSampleSize)
        incProgress(1/2)
        ## update progressionRes results
        obj$progressionRes <- diffmapRes
        
        ## update the project name
        obj$projectName <- paste0(obj$projectName, "_added_diffusionmap")
        
        v$data <- obj
        incProgress(1/2)
      })
      p$progressionCluster <- input$p_clusterMethod
      ## jump to P_tab1
      updateTabsetPanel(session, "P_progressionTabs", selected = "P_tab1")
    }
  })
}

# app = shinyApp(ui = shiny_UI, server = shinyServer)
# runApp(app, host = '0.0.0.0',
#        port = 5678)


































