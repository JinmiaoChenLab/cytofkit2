library(shiny)
# library(shinyFiles)
# source("plot_functions.R")
# library(Seurat)
library(data.table)
source('./ui.R')
library(flowCore)
library(shinyalert)
source('./global.R')
library(shinyWidgets)
# library(ezTools)
# library(monocle)
# library(jsTree)

# maximum size 10000MB
options(shiny.maxRequestSize = 1024*1024*100*100) 
# plan('multicore')

gatePlot <- function (outputId) {
  HTML(paste("<div id=\"", outputId, "\" class=\"shiny-gateplot\"><canvas id=\"gatePlotCanvas\"></canvas></div>", sep=""))
}

render_beadremoval_ui <- function(working.directory, ...) {renderUI({
  fluidPage(
    fluidRow(
      column(12,
             selectizeInput("beadremovalui_beads_type", "Select beads type", multiple = FALSE, width = "100%",
                            choices = c("Fluidigm Beads (140,151,153,165,175)", "Beta Beads (139,141,159,169,175)")),
             selectizeInput("beadremovalui_selected_fcs", "Select FCS file",
                            choices = c("", list.files(file.path(working.directory, "normed"), pattern = "*.fcs$", ignore.case = T)), multiple = FALSE, width = "100%"),
             numericInput("beadremovalui_cutoff", "Cutoff for bead removal", value = 0, min = 0, max = 20),
             actionButton("beadremovalui_remove_beads", "Remove beads (current file)"),
             actionButton("beadremovalui_remove_beads_all_files", "Remove beads (all files)")
             
      )
    )
  )
  
  
})}


render_normalizer_ui <- function(working.directory, ...){renderUI({
  #Remove this fluidpage?
  fluidPage(
    fluidRow(
      column(12,
             selectizeInput("normalizerui_beads_type", "Select beads type", multiple = FALSE, width = "100%",
                            choices = c("Fluidigm Beads (140,151,153,165,175)", "Beta Beads (139,141,159,169,175)")),
             selectizeInput("normalizerui_selected_fcs", "Select FCS file",
                            choices = c("", list.files(working.directory, pattern = "*.fcs$", ignore.case = T)), multiple = FALSE, width = "100%"),
             fluidRow(
               column(6,
                      selectizeInput("normalizerui_baseline", "Select baseline for normalization", multiple = FALSE, width = "100%",
                                     choices = c("Current files", "Existing folder of beads files"))
               ),
               column(6,
                      conditionalPanel(
                        condition <- "input.normalizerui_baseline != 'Current files'",
                        p("Selected baseline beads folder"),
                        verbatimTextOutput("normalizerui_dialog1")
                      )
               )
             ),
             p("You have gated beads for the following files (Only these files will be normalized):"),
             verbatimTextOutput("normalizerui_dialog2"),
             actionButton("normalizerui_visualize_beads", "Visualize beads"),
             actionButton("normalizerui_apply_gates_all_files", "Apply current gates to all files"),
             actionButton("normalizerui_normalize_files", "Normalize")
      )
    )
  )
})}

generate_normalizerui_plot_outputs <- function(n) {renderUI({
  lapply(1:n, function(i) {
    column(2,
           gatePlot(paste("normalizerui_gateplot", i, sep = ""))
    )
  })
  
})}


generate_beadremovalui_plot_outputs <- function(n) {renderUI({
  lapply(1:n, function(i) {
    column(4,
           plotOutput(paste("beadremovalui_plot", i, sep = ""))
    )
  })
})}



shinyServer = function(input, output, session)
{
  
  #### debarcode raw page ####
  # output$debarcoderUI <- render_debarcoder_ui(input, output, session)
  
  
  #debarcoderUI functions
  
  debarcoderui.reactive.values <- reactiveValues(bc.key.fname = file.path(system.file(package = "premessa"), "Fluidigm_20plex_barcode_key.csv"),
                                                 fcs.fname = NULL)
  debarcoderui.reactive.values <- reactiveValues(bc.key.fname = "")
  
  
  debarcoderui_get_bc_key <- reactive({
    bc.key = NULL
    if(is.null(debarcoderui.reactive.values$bc.key.fname) || debarcoderui.reactive.values$bc.key.fname == ""){
      if(is.null(bc.key)){
        if(!is.null(input$pd_used)){
          bcs = gsub(".*?(\\d+).*", "\\1", input$pd_used)
          pops = as.data.frame(t(combn(bcs, input$pd_num)))
          bc.key = as.data.frame(matrix(0, nrow = nrow(pops), ncol = length(input$pd_used)))
          colnames(bc.key) = bcs
          for (i in 1:nrow(bc.key)) {
            bc.key[i, as.character(unlist(pops[i, , drop = T]))] = 1
            rownames(bc.key)[i] = paste0(as.character(unlist(pops[i, , drop = T])), collapse = "_")
          }
          # rownames(bc.key) = paste0('Pop', 1:nrow(bc.key))
        }
      }
    } else {
      bc.key = premessa:::read_barcode_key(debarcoderui.reactive.values$bc.key.fname)
    }
      return(bc.key)
  })
  
  debarcoderui_get_fcs <- reactive({
    ret <- NULL
    if(!is.null(debarcoderui.reactive.values$fcs.fname) && debarcoderui.reactive.values$fcs.fname != "")
      ret <- flowCore::read.FCS(debarcoderui.reactive.values$fcs.fname)
    return(ret)
  })
  
  debarcoderui_get_exprs <- reactive({
    fcs <- debarcoderui_get_fcs()
    
    ret <- NULL
    if(!is.null(fcs)) {
      m <- flowCore::exprs(fcs)
      ret <- asinh(m / 10)
    }
    
    
    return(ret)
  })
  
  debarcoderui_get_mahalanobis_distance <- reactive({
    m <- debarcoderui_get_exprs()
    bc.res <- debarcoderui_get_bc_results()
    if(!is.null(m) && !is.null(bc.res))
      return(premessa:::get_mahalanobis_distance(m, bc.res,
                                                 input$debarcoderui_separation_threshold))
  })
  
  
  debarcoderui_get_bc_results <- reactive({
    # browser()
    bc.key <- debarcoderui_get_bc_key()

    m <- debarcoderui_get_exprs()
    
    if(!is.null(bc.key) && !is.null(m)) {
      res <- premessa:::debarcode_data(m, bc.key)
      
      all.labels <- unique(res$labels)
      all.labels <- sort(all.labels[!is.na(all.labels)])
      updateSelectizeInput(session, "debarcoderui_selected_sample", choices = all.labels)
      
      updateSelectizeInput(session, "debarcoderui_xaxis", choices = res$bc.channels)
      updateSelectizeInput(session, "debarcoderui_yaxis", choices = res$bc.channels)
      
      return(res)
    }
  })
  
  output$debarcoderui_dialog_selected_key <- renderText({
    debarcoderui.reactive.values$bc.key.fname
  })
  
  output$debarcoderui_dialog_selected_fcs <- renderText({
    debarcoderui.reactive.values$fcs.fname
  })
  
  output$debarcoderui_plot1 <- renderPlot({
    if(!is.null(input$debarcoderui_plot_data) && input$debarcoderui_plot_data) {
      withProgress({
        session$sendCustomMessage(type = "plot_loading", "none")
        ret <- NULL
        isolate({
          bc.res <- debarcoderui_get_bc_results()
          if(!is.null(bc.res)) {
            if(input$debarcoderui_plot_type == "Separation")
              ret <- premessa:::plot_separation_histogram(bc.res)
            else if(input$debarcoderui_plot_type == "Event" || input$debarcoderui_plot_type == "Single biaxial"
                    || input$debarcoderui_plot_type == "All barcode biaxials") {
              mahal.dist <- debarcoderui_get_mahalanobis_distance()
              ret <- premessa:::plot_barcode_yields(bc.res, input$debarcoderui_separation_threshold,
                                                    input$debarcoderui_mahal_dist_threshold, mahal.dist)
            }
          }
        })
      }, message = "Plotting figure 1", value = 0.5)
      
      return(ret)
    }
  })
  
  output$debarcoderui_plot2 <- renderPlot({
    if(!is.null(input$debarcoderui_plot_data) && input$debarcoderui_plot_data) {
      withProgress({
        session$sendCustomMessage(type = "plot_loading", "none")
        ret <- NULL
        isolate({
          bc.res <- debarcoderui_get_bc_results()
          if(!is.null(bc.res)) {
            if(input$debarcoderui_plot_type == "Separation")
              ret <- premessa:::plot_barcode_separation(bc.res, input$debarcoderui_separation_threshold)
            else {
              m <- debarcoderui_get_exprs()
              mahal.dist <- debarcoderui_get_mahalanobis_distance()
              m <- cbind(m, mahal.dist = mahal.dist)
              sel.rows <- premessa:::get_sample_idx(input$debarcoderui_selected_sample,
                                                    bc.res, input$debarcoderui_separation_threshold, input$debarcoderui_mahal_dist_threshold, mahal.dist)
              m <- m[sel.rows, ]
              
              if(input$debarcoderui_plot_type == "Event")
                ret <- premessa:::plot_barcode_channels_intensities(m, bc.res$bc.channels, bc.res$m.normed[sel.rows,])
              else {
                if(input$debarcoderui_plot_type == "Single biaxial")
                  ret <- premessa:::plot_color_coded_biaxial(m, input$debarcoderui_xaxis,
                                                             input$debarcoderui_yaxis, "mahal.dist")
                else if(input$debarcoderui_plot_type == "All barcode biaxials")
                  ret <- premessa:::plot_all_barcode_biaxials(m, bc.res$bc.channels)
              }
            }
          }
        })
      }, message = "Plotting figure 2", value = 0.5)
      
      return(ret)
    }
  })
  
  output$debarcoderui_save_files = downloadHandler(
    filename = function() {
      paste0(input$project_name, '_debarcoded.zip')
    },
    content = function(file) {
      browser()
      withProgress({
        fcs <- debarcoderui_get_fcs()
        bc.key <- debarcoderui_get_bc_key()
        if(!is.null(fcs) && !is.null(bc.key)) {
          fcs.fname <- debarcoderui.reactive.values$fcs.fname
          out.dir <- file.path(dirname(fcs.fname), "debarcoded")
          dir.create(out.dir, recursive = T)
          premessa:::debarcode_fcs(fcs, bc.key, out.dir,
                                   tools::file_path_sans_ext(basename(fcs.fname))
                                   , input$debarcoderui_separation_threshold
                                   , input$debarcoderui_mahal_dist_threshold)
          
          files = dir(path = out.dir, pattern = ".*fcs$", full.names = T)
          zip(file, files)
        }
      }, message = "processing data")

    }
  )
  
  # observeEvent(input$debarcoderui_save_files, {
  #   browser()
  #   out.dir = "./"
  #   isolate({
  #     
  #   })
  # })
  
  observeEvent(input$debarcoderui_select_key, {
    isolate({
      debarcoderui.reactive.values$bc.key.fname <- input$debarcoderui_select_key$datapath
    })
  })
  
  observeEvent(input$debarcoderui_select_fcs, {
    isolate({
      debarcoderui.reactive.values$fcs.fname <- input$debarcoderui_select_fcs$datapath
  
      fcs = read.FCS(input$debarcoderui_select_fcs$datapath)
      pd <- fcs@parameters@data
      markers <- paste(pd$name, "<", pd$desc, ">", sep = "")
      updateSelectInput(session, "pd_used", choices = markers)
    })
  })
  #######################################
  
  #### normalizeation raw page ####
  #options(warn = -1)
  working.directory <- "./"
  # working.directory <- dirname(file.choose())
  normed.dir <- file.path(working.directory, "normed")
  beads.removed.dir <- file.path(normed.dir, "beads_removed")
  beadremovalui.plots.number <- 3
  
  output$normalizerUI <- render_normalizer_ui(working.directory, input, output, session)
  output$normalizerUI_plot_outputs <- generate_normalizerui_plot_outputs(5)
  output$beadremovalUI <- render_beadremoval_ui(working.directory, input, output, session)
  output$beadremovalUI_plot_outputs <- generate_beadremovalui_plot_outputs(beadremovalui.plots.number)
  
  #beadremovalUI functions
  
  get_beadremovalui_fcs <- reactive({
    ret <- NULL
    
    if(!is.null(input$beadremovalui_selected_fcs) && input$beadremovalui_selected_fcs != "")
      ret <- flowCore::read.FCS(file.path(normed.dir, input$beadremovalui_selected_fcs))
    
    return(ret)
  })
  
  observeEvent(input$beadremovalui_remove_beads, {
    isolate({
      dir.create(beads.removed.dir, recursive = T)
      showModal(modalDialog(
        title = "Normalizer report",
        "Bead removal started, please wait..."
      ))
      premessa::remove_beads_from_file(file.path(normed.dir, input$beadremovalui_selected_fcs),
                                       input$beadremovalui_cutoff, beads.removed.dir)
      showModal(modalDialog(
        title = "Normalizer report",
        sprintf("Beads removed from file: %s", input$beadremovalui_selected_fcs)
      ))
    })
  }
  )
  
  observeEvent(input$beadremovalui_remove_beads_all_files, {
    isolate({
      dir.create(beads.removed.dir, recursive = T)
      files.list <- list.files(normed.dir, pattern = "*.fcs$", ignore.case = T)
      showModal(modalDialog(
        title = "Normalizer report",
        "Bead removal started, please wait..."
      ))
      files.list <- lapply(files.list, function(f.name) {
        fcs <- flowCore::read.FCS(file.path(normed.dir, f.name))
        premessa::remove_beads_from_file(file.path(normed.dir, f.name), input$beadremovalui_cutoff, beads.removed.dir)
        return(f.name)
      })
      showModal(modalDialog(
        title = "Normalizer report",
        p("Beads removed from files:", br(),
          lapply(files.list, function(x) list(x, br()))
        )
      ))
    })
  }
  )
  
  observe({
    if(!is.null(input$beadremovalui_selected_fcs) && input$beadremovalui_selected_fcs != "") {
      fcs <- flowCore::read.FCS(file.path(normed.dir, input$beadremovalui_selected_fcs))
      beads.type <-  premessa:::get_beads_type_from_description(input$beadremovalui_beads_type)
      
      beads.cols.names <- premessa:::find_beads_channels_names(fcs, beads.type)
      combs <- rep(beads.cols.names, length.out = beadremovalui.plots.number * 2)
      
      m <- flowCore::exprs(fcs)
      
      lapply(seq(1, length(combs), 2), function(i) {
        plot.idx <- ceiling(i / 2)
        plot.output <- premessa:::plot_distance_from_beads(m, combs[i], combs[i + 1])
        output[[paste("beadremovalui_plot", plot.idx, sep ="")]] <- renderPlot(plot.output)
      })
    }
  })
  
  
  
  #normalizerUI functions
  
  beads.gates <- reactiveValues()
  normalizerui.baseline.dir <- NULL
  
  output$normalizerui_dialog1 <- renderText({
    if(!is.null(input$normalizerui_baseline) && length(grep("^Existing", input$normalizerui_baseline)) > 0) {
      normalizerui.baseline.dir <<- dirname(file.choose())
      return(normalizerui.baseline.dir)
    }
  })
  
  output$normalizerui_dialog2 <- renderText({
    paste(names(beads.gates), collapse = ", ")
  })
  
  
  observeEvent(input$normalizerui_apply_gates_all_files, {
    isolate({
      cur.gates <- get_beads_gates_for_current_file()
      files.list <- list.files(working.directory, pattern = "*.fcs$", ignore.case = T)
      lapply(files.list, function(f.name) {
        beads.gates[[f.name]] <- cur.gates
        NULL
      })
    })
  })
  
  get_beads_gates_for_current_file <- reactive({
    if(!input$normalizerui_selected_fcs %in% names(beads.gates)) {
      fcs <- get_fcs()
      beads.gates[[input$normalizerui_selected_fcs]] <<- premessa:::get_initial_beads_gates(fcs)
    }
    beads.gates[[input$normalizerui_selected_fcs]]
  })
  
  
  get_fcs <- reactive({
    ret <- NULL
    
    if(!is.null(input$normalizerui_selected_fcs) && input$normalizerui_selected_fcs != "")
      ret <- flowCore::read.FCS(file.path(working.directory, input$normalizerui_selected_fcs))
    
    return(ret)
  })
  
  get_exprs <- reactive({
    fcs <- get_fcs()
    m <- flowCore::exprs(fcs)
    m <- asinh(m / 5)
    if(nrow(m) > 50000) {
      m <- m[sample(1:nrow(m), 50000),]
    }
    return(m)
    
  })
  
  get_beads_type <- reactive({
    premessa:::get_beads_type_from_description(input$normalizerui_beads_type)
  })
  
  
  
  do_plot_outputs <- function(sel.beads = NULL) {
    fcs <- get_fcs()
    beads.type <- get_beads_type()
    beads.cols <- premessa:::find_bead_channels(fcs, beads.type)
    dna.col <- premessa:::find_dna_channel(fcs)
    
    gates <- isolate({get_beads_gates_for_current_file()})
    
    m <- get_exprs()
    colors <- rep("black", nrow(m))
    if(!is.null(sel.beads))
      colors[sel.beads] <- "red"
    
    #Needs to be in lapply to work
    #see https://github.com/rstudio/shiny/issues/532
    lapply(1:length(beads.cols), function(i) {
      xAxisName <- premessa:::get_parameter_name(fcs, beads.cols[i])
      yAxisName <- premessa:::get_parameter_name(fcs, dna.col)
      
      output[[paste("normalizerui_gateplot", i, sep ="")]] <- reactive({
        list(
          x = m[, beads.cols[i]],
          y = m[, dna.col],
          color = colors,
          xAxisName = xAxisName,
          yAxisName = yAxisName,
          file = input$normalizerui_selected_fcs,
          channelGates = gates[[xAxisName]]
        )
      })
    })
  }
  
  observe({
    fcs <- get_fcs()
    if(!is.null(fcs)) {
      do_plot_outputs()
    }
    
  })
  
  observeEvent(input$normalizerui_normalize_files, {
    isolate({
      beads.type <- get_beads_type()
      baseline <- NULL
      if(length(grep("^Existing", input$normalizerui_baseline)) > 0)
        baseline <- normalizerui.baseline.dir
      
      showModal(modalDialog(
        title = "Normalizer report",
        "Normalization started, please wait..."
      ))
      premessa::normalize_folder(working.directory, "normed",
                                 reactiveValuesToList(beads.gates), beads.type, baseline = baseline)
      updateSelectizeInput(session, input$beadremovalui_selected_fcs,
                           choices = c("", list.files(normed.dir, pattern = "*normalized.fcs$", ignore.case = T)))
      showModal(modalDialog(
        title = "Normalizer report",
        sprintf("Normalization complete! The output files are located in %s",
                normed.dir)
      ))
    })
    
  })
  
  observeEvent(input$normalizerui_visualize_beads, {
    withProgress({
      isolate({
        fcs <- get_fcs()
        m <- get_exprs()
        dna.col <- premessa:::find_dna_channel(fcs)
        
        gates <- get_beads_gates_for_current_file()
        beads.type <- get_beads_type()
        beads.cols <- premessa:::find_bead_channels(fcs, beads.type)
        beads.cols.names <- premessa:::get_parameter_name(fcs, beads.cols)
        
        sel <- premessa:::identify_beads(m, gates, beads.cols.names, dna.col)
        do_plot_outputs(sel)
      })
    }, message = "Visuallizing beads", value = 0.5)
    
  })
  
  observe({
    if(!is.null(input$normalizerui_gate_selected)) {
      isolate({
        #Change this to use get_beads_gates_for_current_file ??
        gate.data <- input$normalizerui_gate_selected
        temp <- beads.gates[[input$normalizerui_selected_fcs]][[gate.data$xAxisName]]
        
        temp$x <- unlist(gate.data$xLim)
        temp$y <- unlist(gate.data$yLim)
        
        beads.gates[[input$normalizerui_selected_fcs]][[gate.data$xAxisName]] <<- temp
        
      })
    }
    
  })
  #################################
  v = reactiveValues(
    # fcsFile = NULL
    markers = NULL
    , selected_markers = NULL
    , export_figure_selected = list()
    , data = NULL
    , sampleInfo = NULL
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
    old_name = sort(dir(path = cur_dir
                        , pattern = "*.fcs", full.names = T))
    file.rename(old_name, paste0(cur_dir, input$rawfcs$name))
    inputs$fcsFiles = paste0(cur_dir, input$rawfcs$name)
    # input$rawfcs$datapath = paste0(cur_dir, input$rawfcs$name)
    
    fcs <- suppressWarnings(read.FCS(inputs$fcsFiles[1]))
    pd <- fcs@parameters@data
    markers <- paste(pd$name, "<", pd$desc, ">", sep = "")
    inputs$markers = markers
    v$markers = markers
    updateSelectInput(session, "markers", choices = markers)
    updatePickerInput(session, 'pickermarker', choices = markers)
    if (length(markers) == 0) {
      shinyalert(title = "Error", text = "No markers found in the FCS file!", type = "error")
      # stop()
    }
  })
  
  observeEvent(input$select_all_markers, {
    updateSelectInput(session, "markers", selected = v$markers)
  })
  observeEvent(input$deselect_all_markers, {
    updateSelectInput(session, "markers", selected = "")
  })
  
  observeEvent(input$submit, {
    # browser()
    ### need to check data
    withProgress(
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
               resultDir = "./",
               saveResults = TRUE,
               saveObject = TRUE,
               l_w = as.numeric(inputs[["l_w"]]), 
               l_t = as.numeric(inputs[["l_t"]]), 
               l_m = as.numeric(inputs[["l_m"]]), 
               l_a = as.numeric(inputs[["l_a"]]))
      , value = 0.5, message = "Analysing")
  })
  
  output$download_analysis_res = downloadHandler(
    filename = function() {
      paste0(input$project_name, '.zip')
    },
    content = function(file) {
      # browser()
      cur_dir = './'
      files = dir(path = cur_dir, pattern = input$project_name)
      zip(file, files)
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
      cat(cytofkitObj$datapath)
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
  
  output$sampleSelect <- renderUI({
    if(is.null(v$data) || is.null(v$sampleInfo)){
      return(NULL)
    }else{
      sampleNames <- unique(as.character(v$sampleInfo$cellSample))
      checkboxGroupInput('samples', NULL, 
                         sampleNames, selected = sampleNames)
    }   
  })
  
  observeEvent(input$selectDeselectAll, {
    allSamp <- input$selectDeselectAll
    sampleNames <- unique(as.character(v$sampleInfo$cellSample))
    if(allSamp == TRUE){
      updateCheckboxGroupInput(session, 'samples', selected = sampleNames)
    }else{
      updateCheckboxGroupInput(session, 'samples', selected = character(0))
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
        ## check results saving path
        if(is.null(v$data$resultDir) || !dir.exists(v$data$resultDir)){
          v$data$resultDir <- path.expand("~")  ## default save to home if not specified
        }
        saveToFCS <- input$saveFCS
        if(is.null(v$data$rawFCSdir)){
          saveToFCS <- FALSE
          warning("Path for original FCS files is not provided, 
                  data cannnot be saved to new copies of FCS files.")
        }else if(!dir.exists(v$data$rawFCSdir)){
          saveToFCS <- FALSE
          warning(paste0("Path for original FCS files doesn't exist, 
                         data cannnot be saved to new copies of FCS files.", 
                         "Please check path: ", v$data$rawFCSdir))
        }
        
        ## NOTE: if samples are regrouped, then new FCS file cannot be saved
        incProgress(1/2, message = paste0("To ", v$data$resultDir))
        v$data$sampleInfo <- v$sampleInfo
        analysis_results <<- v$data
        cytof_writeResults(analysis_results,
                           saveToRData = input$saveRData,
                           saveToFCS = saveToFCS,
                           saveToFiles = input$saveCsv)
        incProgress(1/2)
        ## open the results directory
        opendir(v$data$resultDir)
      })
      
      cur_dir = './'
      files = dir(path = cur_dir, pattern = input$project_name)
      zip(file, files)
    }
  )
  
  observeEvent(input$OpenDir, {
    pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
    if(dir.exists(pdfDir)){
      opendir(pdfDir)
    }else{
      stop("PDF not created yet!")
    }
  })
  
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
  
  observeEvent(input$PDFClusterPlot, {
    if(!is.null(v$data)){
      withProgress(message="Downloading Clusterplot PDF files...", value=0, {
        print(getwd())
        dir.create(paste0("cytofkit_PDF_Plots_", Sys.Date()))
        pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
        filename1 <- paste0(pdfDir, .Platform$file.sep, "cytofkit_shinyAPP_Clusterplot_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename1)){
          filename1 <- paste0(pdfDir, .Platform$file.sep,
                              "cytofkit_shinyAPP_Clusterplot_",
                              Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        pdf(filename1, 
            width=as.integer(input$tab_w), 
            height=as.integer(input$tab_h))
        C_ScatterPlotInput()
        dev.off()
      })
    }
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
        if(ilabel == ""){
          clusterLabels[clusterLabels==clusteri] <- "Unknown"
        }else{
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
    # dev.copy2pdf(file = "cytofkit_shinyAPP_marker_heatmap.pdf",
    #              width=as.integer(input$tab_w), 
    #              height=as.integer(input$tab_h))
  })
  
  output$M_heatmapPlot <- renderPlot({
    M_heatmapPlotInput()
  }, height = 900, width = 950)
  
  observeEvent(input$PDFHeatmap, {
    if(!is.null(v$data)){
      withProgress(message="Downloading Marker Heatmap PDF files...", value=0, {
        print(getwd())
        dir.create(paste0("cytofkit_PDF_Plots_", Sys.Date()))
        pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
        filename1 <- paste0(pdfDir, .Platform$file.sep, "cytofkit_shinyAPP_Marker_Heatmap_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename1)){
          filename1 <- paste0(pdfDir, .Platform$file.sep,
                              "cytofkit_shinyAPP_Marker_Heatmap_",
                              Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        file.copy("cytofkit_shinyAPP_marker_heatmap.pdf", filename1)
      })
    }
  })
  
  session$onSessionEnded(function(){
    file.remove("cytofkit_shinyAPP_marker_heatmap.pdf")
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
  
  observeEvent(input$PDFExpPlot, {
    if(!is.null(v$data)){
      withProgress(message="Downloading Marker Expression Plot PDF files...", value=0, {
        print(getwd())
        dir.create(paste0("cytofkit_PDF_Plots_", Sys.Date()))
        pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
        filename1 <- paste0(pdfDir, .Platform$file.sep, "cytofkit_shinyAPP_Marker_Expression_Plot_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename1)){
          filename1 <- paste0(pdfDir, .Platform$file.sep,
                              "cytofkit_shinyAPP_Marker_Expression_Plot_",
                              Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        pdf(filename1, 
            width=as.integer(input$tab_w), 
            height=as.integer(input$tab_h))
        M_markerExpressionPlotInput()
        dev.off()
      })
    }
  })
  
  
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
  
  observeEvent(input$PDFHistogram, {
    if(!is.null(v$data)){
      withProgress(message="Downloading Stack Density Plot PDF files...", value=0, {
        print(getwd())
        dir.create(paste0("cytofkit_PDF_Plots_", Sys.Date()))
        pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
        filename1 <- paste0(pdfDir, .Platform$file.sep, "cytofkit_shinyAPP_Stack_Density_Plot_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename1)){
          filename1 <- paste0(pdfDir, .Platform$file.sep,
                              "cytofkit_shinyAPP_Stack_Density_Plot_",
                              Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        pdf(filename1, 
            width=as.integer(input$tab_w), 
            height=as.integer(input$tab_h))
        M_stackDensityPlotInput()
        dev.off()
      })
    }
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
      markerNames <- colnames(v$data$expressionData)
      newMarkerNames <- NULL
      for (i in 1:length(markerNames)){
        iName <- input[[paste0('marker_', i, '_name')]]
        newMarkerNames <- c(newMarkerNames, iName)
      }
      ## update new cluster colours
      colnames(v$data$expressionData) <- newMarkerNames
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
    
    dev.copy2pdf(file = "cytofkit_shinyAPP_cells_heatmap_plot_plot.pdf",
                 width=as.integer(input$tab_w), 
                 height=as.integer(input$tab_h))
  })
  
  output$S_heatmapPlot <- renderPlot({
    S_heatmapPlotInput()
  }, height = 900, width = 950)
  
  observeEvent(input$PDFSamHeat, {
    if(!is.null(v$data)){
      withProgress(message="Downloading Sample Heatmap PDF files...", value=0, {
        print(getwd())
        dir.create(paste0("cytofkit_PDF_Plots_", Sys.Date()))
        pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
        filename1 <- paste0(pdfDir, .Platform$file.sep, "cytofkit_shinyAPP_Sample_Heatmap_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename1)){
          filename1 <- paste0(pdfDir, .Platform$file.sep,
                              "cytofkit_shinyAPP_Sample_Heatmap_",
                              Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        file.copy("cytofkit_shinyAPP_cells_heatmap_plot_plot.pdf", filename1)
      })
    }
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
  
  observeEvent(input$PDFrateChange, {
    if(!is.null(v$data)){
      withProgress(message="Downloading Rate Change Plot PDF files...", value=0, {
        print(getwd())
        dir.create(paste0("cytofkit_PDF_Plots_", Sys.Date()))
        pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
        filename1 <- paste0(pdfDir, .Platform$file.sep, "cytofkit_shinyAPP_Rate_Change_Plot_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename1)){
          filename1 <- paste0(pdfDir, .Platform$file.sep,
                              "cytofkit_shinyAPP_Rate_Change_Plot_",
                              Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        pdf(filename1, 
            width=as.integer(input$tab_w), 
            height=as.integer(input$tab_h))
        S_rateChangePlotInput()
        dev.off()
      })
    }
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
  
  ## update sample groups
  observeEvent(input$updateSampleGroups, {
    if(!is.null(v$data) && !is.null(v$sampleInfo)){
      v$sampleInfo$originalCellSample <- v$sampleInfo$cellSample
      uniqueSampleNames <- sort(unique(v$sampleInfo$originalCellSample))
      
      sampleGroupNames <- NULL
      for(i in 1:length(uniqueSampleNames)){
        sampleGroupNames <- c(sampleGroupNames, input[[paste0("Sample", i)]])
        v$data$sampleNames[[i]] <- c(v$data$sampleNames[[i]], input[[paste0("Sample", i)]])
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
  
  observeEvent(input$PDFScatter, {
    if(!is.null(v$data)){
      withProgress(message="Downloading Progression Scatterplot PDF files...", value=0, {
        print(getwd())
        dir.create(paste0("cytofkit_PDF_Plots_", Sys.Date()))
        pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
        filename1 <- paste0(pdfDir, .Platform$file.sep, "cytofkit_shinyAPP_Scatterplot_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename1)){
          filename1 <- paste0(pdfDir, .Platform$file.sep,
                              "cytofkit_shinyAPP_Scatterplot_",
                              Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        pdf(filename1, 
            width=as.integer(input$tab_w), 
            height=as.integer(input$tab_h))
        P_scatterPlotInput()
        dev.off()
      })
    }
  })
  
  
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
  
  observeEvent(input$PDFmarkerPlot, {
    if(!is.null(v$data)){
      withProgress(message="Downloading Marker Plot PDF files...", value=0, {
        print(getwd())
        dir.create(paste0("cytofkit_PDF_Plots_", Sys.Date()))
        pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
        filename1 <- paste0(pdfDir, .Platform$file.sep, "cytofkit_shinyAPP_Marker_Plot_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename1)){
          filename1 <- paste0(pdfDir, .Platform$file.sep,
                              "cytofkit_shinyAPP_Marker_Plot_",
                              Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        pdf(filename1, 
            width=as.integer(input$tab_w), 
            height=as.integer(input$tab_h))
        P_markerPlotInput()
        dev.off()
      })
    }
  })
  
  
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
  
  
  ###### export figure page #############
  # browser()
  observeEvent(input$sbm, {
    # browser()
    if(input$sbm == "report_gen"){
      
      ### Cluster plot parameters
      updateSelectInput(session, inputId = "export_cluster_plot_vis_method"
                        , choices = visualizationMethods())
      updateSelectInput(session, inputId = "export_cluster_plot_cluster_by"
                        , choices = visualizationFunctions())
      clusterMethod <- input$c_PlotFunction
      clusterIDs <- sort(unique(v$data$clusterRes[[clusterMethod]]))
      updateSelectInput(session, inputId = "export_cluster_plot_cluster_filter"
                        , choices = clusterIDs, selected = clusterIDs)
      ### Expression heatmap parameters
      updateSelectInput(session, inputId = "export_expression_heatmap_cluster_method"
                        , choices = clusterMethods())
      sorted_markerNames <- colnames(v$data$expressionData)
      markerNames <- sorted_markerNames[order(sorted_markerNames)]
      initNum <- ifelse(length(markerNames) >=4, 4, 1)
      updateSelectInput(session, inputId = "export_expression_heatmap_select_markers"
                        , choices = markerNames, selected = markerNames[1:initNum])
      
      ### Expression level plot
      updateSelectInput(session, inputId = "export_expression_level_vis_method"
                        , choices = visualizationMethods())
      
      sorted_markers <- colnames(v$data$expressionData)
      sorted_markers <- sorted_markers[order(sorted_markers)]
      #markers <- c(sorted_markers, "All Markers", "All Markers(scaled)")
      updateSelectInput(session, inputId = "export_expression_level_markers"
                        , choices = sorted_markers, selected = sorted_markers[1])
      ### Expression histogram
      stackFactorChoice <- c(names(v$data$clusterRes), "sample") 
      updateSelectInput(session, inputId = "export_expression_hist_stack_factor"
                        , choices = stackFactorChoice, selected = stackFactorChoice[1])
      
      sorted_markerNames <- colnames(v$data$expressionData)
      markerNames <- sorted_markerNames[order(sorted_markerNames)]
      initNum <- ifelse(length(markerNames) >=4, 4, 1)
      updateSelectInput(session, inputId = "export_expression_hist_select_markers"
                        , choices = markerNames, selected = markerNames[1:initNum])
      
    }
    
  })
  
  ## reponse to add figure
  observeEvent(input$add_report_figure, {
    # browser()
    num = length(v$export_figure_selected) + 1
    v$export_figure_selected[[num]] = list()
    v$export_figure_selected[[num]][["figure"]] = input$export_fig_selection
  })
  
  observeEvent(v$export_figure_selected, {
    figures = c()
    for (i in 1:length(v$export_figure_selected)) {
      figures[i] = v$export_figure_selected[[i]]$figure
    }
    figures = paste(figures, collapse = "\n")
    updateTextAreaInput(session, 'figure_selected', value = figures)
    max_i <<- length(v$export_figure_selected)
  })
  
  output$report_preview <- renderUI({
    plot_output_list <- lapply(1:length(v$export_figure_selected), function(i) {
      plotname <- paste("plot", i, sep="")
      plotOutput(plotname, height = 280, width = 250)
    })
    
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
  })
  
  # Call renderPlot for each one. Plots are only actually generated when they
  # are visible on the web page.
  for (i in 1:max_i) {
    # Need local so that each item gets its own number. Without it, the value
    # of i in the renderPlot() will be the same across all instances, because
    # of when the expression is evaluated.
    local({
      my_i <- i
      plotname <- paste("plot", my_i, sep="")
      
      output[[plotname]] <- renderPlot({
        plot(1:my_i, 1:my_i,
             xlim = c(1, max_i),
             ylim = c(1, max_i),
             main = paste("1:", my_i, ".  n is ", input$n, sep = "")
        )
      })
    })
  }
  # output$figure_selected <- renderText({
  #   # browser()
  # 
  # })
}

# app = shinyApp(ui = shiny_UI, server = shinyServer)
# runApp(app, host = '0.0.0.0',
#        port = 5678)


































