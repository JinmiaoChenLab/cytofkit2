#' Function for launching the user friendly GUI client for \code{cytofkit-package}
#' 
#' This GUI provides an easy way to apply \code{cytofkit} package. 
#' Main parameters for running 'cytofkit' main function were integrated in this GUI, 
#' and each parameter has a help button to show the instructions. 
#' The \code{cytofkit} analysis will be automatically started after submitting.
#' 
#' @author Hao Chen
#' @return the GUI for \code{cytofkit-package}
#' @import tcltk
#' @export
#' @seealso \code{\link{cytofkit-package}}, \code{\link{cytofkit}}
#' @references \url{http://signbioinfo.github.io/cytofkit/}
#' @examples
#' #cytofkit_GUI()  # remove the hash symbol to run
cytofkit_GUI <- function() {
    
    ##--------------------------##
    ## parameter initialization ##
    ##--------------------------##
    
    fcsFiles <- ""
    cur_dir <- getwd()
    mergeMethods <- c("all", "min", "ceil", "fixed")
    transformMethods <- c("autoLgcl", "cytofAsinh", "logicle", "none")
    vizMethods <- c("pca", "isomap", "tsne", "NULL")
    clusterMethods <- c("Rphenograph", "ClusterX", "DensVM", "FlowSOM", "NULL")
    progressionMethods <- c("diffusionmap", "isomap", "NULL")
    fixedLgclParas = c(l_w = 0.5, l_t = 500000, l_m = 4.5, l_a = 0)
    
    rawFCSdir <- tclVar(cur_dir)
    fcsFile <- tclVar("")
    resDir <- tclVar(cur_dir)
    projectName <- tclVar("cytofkit")
    mergeMethod <- tclVar("ceil")
    fixedNum <- tclVar("5000")
    markers <- tclVar("")
    transformMethod <- tclVar("autoLgcl")
    progressionMethod <- tclVar("NULL")
    Rphenograph_k <- tclVar("30")
    tsne_perp <- tclVar("30")
    tsne_maxIter <- tclVar("1000")
    FlowSOM_k <- tclVar("40")
    seed <- tclVar("42")
    
    # logicle parameters
    l_w <- tclVar(fixedLgclParas[1])
    l_t <- tclVar(fixedLgclParas[2])
    l_m <- tclVar(fixedLgclParas[3])
    l_a <- tclVar(fixedLgclParas[4])
    
    clusterSelect <- c()
    i <- 1
    while (i <= length(clusterMethods)) {
        aux <- paste("clusterSelect", i, sep = "")
        clusterSelect <- c(clusterSelect, as.character(aux))
        tclvalue(clusterSelect[i]) <- "0"
        i <- i + 1
    }
    tclvalue(clusterSelect[1]) <- "1"  # default Rphenograph
    
    vizSelect <- c()
    i <- 1
    while (i <= length(vizMethods)) {
        aux <- paste("vizSelect", i, sep = "")
        vizSelect <- c(vizSelect, as.character(aux))
        tclvalue(vizSelect[i]) <- "0"
        i <- i + 1
    }
    tclvalue(vizSelect[3]) <- "1"  # default tsne 
    
    ret_var <- tclVar("")
    
    ##-------------------##
    ##  button functions ##
    ##-------------------##
    
    highCell_warning <- function() {
        #potentially_slow_methods <- c("all", "DensVM", "isomap")
        #if(any(potentially_slow_methods) %in% para){}
        tkmessageBox(title = "Warning",
                     message = "Please note that using the options DensVM or isomap with more than 10,000 cells will be very slow!")
    }
    
    reset_rawFCS_dir <- function() {
        rawFCS_dir <- ""
        rawFCS_dir <- tclvalue(tkchooseDirectory(title = "Choose your rawFCS direcetory ..."))
        if (rawFCS_dir != "") {
            tclvalue(rawFCSdir) <- rawFCS_dir
            tclvalue(resDir) <- rawFCS_dir
        }
    }
    
    reset_res_dir <- function() {
        res_dir <- ""
        res_dir <- tclvalue(tkchooseDirectory(title = "Choose your result directory ..."))
        if (res_dir != "") {
            tclvalue(resDir) <- res_dir
        }
    }
    
    reset_fcs_data <- function() {
        fnames <- ""
        fnames <- tk_choose.files(default = paste(tclvalue(rawFCSdir), 
            "fcs", sep = .Platform$file.sep), caption = "Select FCS files", 
            multi = TRUE, filters = matrix(c("{fcs files}", "{.fcs}"), 
                1, 2), index = 1)
        if (length(fnames) >= 1) {
            fnames <- fnames[!(grepl(paste0(.Platform$file.sep, 
                "fcs$"), fnames))]  # remove empty .fcs files
            tclvalue(fcsFile) <- paste(fnames, collapse = "}{")
        }
    }
    
    reset_para_data <- function() {
        
        if (tclvalue(fcsFile) == "" && tclvalue(rawFCSdir) == "") {
            tkmessageBox(title = "cytofkit: an Integrated Mass Cytometry Data Analysis Pipeline", 
                message = "Please input your \"rawFCSdirectory\" or \"fcsFile\".", 
                icon = "info", type = "ok")
        }
        if (tclvalue(fcsFile) != "") {
            fcsFiles <- strsplit(tclvalue(fcsFile), "}{", fixed = TRUE)[[1]]
        }
        fcsDir <- tclvalue(rawFCSdir)
        selectMarkers <- getParameters_GUI(fcsFiles, fcsDir)
        
        if (length(selectMarkers) > 0) {
            tclvalue(markers) <- paste(selectMarkers, collapse = "}{")
        }
    }
    
    reset_num2null <- function() {
        tclvalue(fixedNum) <- "NULL"
    }
    
    method_all_warning <- function(){
        reset_num2null()
        highCell_warning()
    }
    
    reset_num2any <- function() {
        tclvalue(fixedNum) <- "5000"
    }
    
    rawFCSdir_help <- function() {
        tkmessageBox(title = "rawFCSdir", message = "The directory that contains fcs files.", 
            icon = "info", type = "ok")
    }
    
    fcsFile_help <- function() {
        tkmessageBox(title = "fcsFiles", message = "The fcs files to be analyzed. One or multiple fcs files are allowed. When multiple fcs files are selected, cells from each fcs file are combined for analysis.", 
            icon = "info", type = "ok")
    }
    
    para_help <- function() {
        tkmessageBox(title = "markers", message = "Select the list of makers to be used for analysis.", 
            icon = "info", type = "ok")
    }
    
    projectName_help <- function() {
        tkmessageBox(title = "projectName", message = "A prefix that will be added to the names of result files.", 
            icon = "info", type = "ok")
    }
    
    mergeMethod_help <- function() {
        tkmessageBox(title = "mergeMethod", message = "When multiple fcs files are selected, cell expression data can be merged using one of the four different methods including \"ceil\",\"all\", \"min\",\"fixed\". \n\n\"ceil\" (the default option): up to a fixed number (specified by fixedNum) of cells are sampled without replacement from each fcs file and combined for analysis. \n\n\"all\": all cells from each fcs file are combined for analysis. \n\n\"min\": The minimum number of cells among all the selected fcs files are sampled from each fcs file and combined for analysis. \n\n\"fixed\": a fixed num (specified by fixedNum) of cells are sampled (with replacement when the total number of cell is less than fixedNum) from each fcs file and combined for analysis.", 
            icon = "info", type = "ok")
    }
    
    fixedNum_help <- function() {
        tkmessageBox(title = "fixedNum", message = "Up to fixedNum of cells from each fcs file are used for analysis.", 
            icon = "info", type = "ok")
    }
    
    rPk_help <- function() {
        tkmessageBox(title = "Rphenograph K", message = "Number of nearest neighbours to pass to Rphenograph.", 
                     icon = "info", type = "ok")
    }
    
    fSk_help <- function() {
        tkmessageBox(title = "FlowSOM K", message = "Number of clusters for meta clustering in FlowSOM.", 
                     icon = "info", type = "ok")
    }
    
    transformMethod_help <- function() {
        tkmessageBox(title = "transformationMethod", message = "Data Transformation method, including \"cytofAsinh\"(Customized Asinh transformation for CyTOF data), \"autoLgcl\"(automatic logicle transformation for CyTOF data), \"logicle\"(customize your own parameters for logicle transformation) and \"none\"(if your data is already transformed).", 
            icon = "info", type = "ok")
    }
    
    cluster_help <- function() {
        tkmessageBox(title = "clusterMethods", message = "The method(s) for clustering, including \"DensVM\", \"ClusterX\", \"Rphenograph\", and \"FlowSOM\". \n\nIf \"NULL\" was selected, no clustering will be performed.", 
            icon = "info", type = "ok")
    }
    
    visualizationMethods_help <- function() {
        tkmessageBox(title = "visualizationMethods", message = "The method(s) used for visualizing the clustering results, multiple selections are allowed. Including \"pca\", \"isomap\", \"tsne\". \n\nWARNING: \"tsne\" is the default selection, \"isomap\" may take long time.", 
            icon = "info", type = "ok")
    }
    
    progressionMethod_help <- function() {
        tkmessageBox(title = "progressionMethod", message = "The method used for cellular progression analysis including \"diffusion map\" and \"isomap\"\n\nIf \"NULL\" was selected, no progression estimation will be performed.", 
                     icon = "info", type = "ok")
    }
    
    resDir_help <- function() {
        tkmessageBox(title = "resDir", message = "The directory where result files will be generated", 
            icon = "info", type = "ok")
    }
    
    reset <- function() {
        tclvalue(rawFCSdir) = cur_dir
        tclvalue(fcsFile) = ""
        tclvalue(resDir) = cur_dir
        tclvalue(projectName) = "cytofkit"
        tclvalue(mergeMethod) = mergeMethods[3]
        tclvalue(fixedNum) = "5000"
        tclvalue(markers) = ""
        tclvalue(transformMethod) = "autoLgcl"
        tclvalue(clusterSelect[1]) = "0"
        tclvalue(clusterSelect[2]) = "1"
        tclvalue(clusterSelect[3]) = "0"
        tclvalue(clusterSelect[4]) = "0"
        tclvalue(vizSelect[1]) <- "0"
        tclvalue(vizSelect[2]) <- "0"
        tclvalue(vizSelect[3]) <- "1"
        tclvalue(progressionMethod) <- "NULL"
        tclvalue(Rphenograph_k) <- "30"
        tclvalue(FlowSOM_k) <- "40"
        tclvalue(tsne_perp) <- "30"
        tclvalue(tsne_maxIter) <- "1000"
        tclvalue(seed) <- "42"
    }
    
    submit <- function() {
        has_error = FALSE
        if (tclvalue(markers) == "") {
            tkmessageBox(title = "cytofkit: an Integrated Analysis Pipeline for Mass Cytometry Data", 
                message = "Please select the markers for your analysis.", 
                icon = "info", type = "ok")
            has_error = TRUE
        }
        
        if (has_error == FALSE) {
            tclvalue(ret_var) <- "OK"
            tkdestroy(tt)
        }
    }
    
    quit <- function() {
        tkdestroy(tt)
    }
    
    
    ##----------------##
    ##  build the GUI ##
    ##--------------- ##
    
    ## head line
    tt <- tktoplevel(borderwidth = 20)
    tkwm.title(tt, "cytofkit: an Integrated Analysis Pipeline for Mass Cytometry Data")
    
    if(.Platform$OS.type == "windows"){
        box_length <- 63
    }else{
        box_length <- 55 
    }
    cell_width <- 3
    bt_width <- 8
    #hb_width <- 8
    
    imgfile <- system.file("extdata", "help.png", package = "cytofkit")
    image1 <- tclVar()
    tkimage.create("photo", image1, file = imgfile)
    image2 <- tclVar()
    tkimage.create("photo", image2)
    tcl(image2, "copy", image1, subsample = 6)
    
    ## rawFCSdir
    rawFCSdir_label <- tklabel(tt, text = "Raw FCS Directory :")
    rawFCSdir_entry <- tkentry(tt, textvariable = rawFCSdir, width = box_length)
    rawFCSdir_button <- tkbutton(tt, text = " Choose... ", width = bt_width, command = reset_rawFCS_dir)
    rawFCSdir_hBut <- tkbutton(tt, image = image2, command = rawFCSdir_help)
    
    ## resDir
    resDir_label <- tklabel(tt, text = "Result Directory :")
    resDir_entry <- tkentry(tt, textvariable = resDir, width = box_length)
    resDir_button <- tkbutton(tt, text = " Choose... ", width = bt_width, 
        command = reset_res_dir)
    resDir_hBut <- tkbutton(tt, image = image2, command = resDir_help)
    
    ## fcsFiles
    fcsFile_label <- tklabel(tt, text = "FCS File(s) :")
    fcsFile_entry <- tkentry(tt, textvariable = fcsFile, width = box_length)
    fcsFile_button <- tkbutton(tt, text = " Select... ", width = bt_width, 
        command = reset_fcs_data)
    fcsFile_hBut <- tkbutton(tt, image = image2, command = fcsFile_help)
    
    ## markers
    markers_label <- tklabel(tt, text = "Markers :")
    markers_entry <- tkentry(tt, textvariable = markers, width = box_length)
    markers_button <- tkbutton(tt, text = " Select... ", width = bt_width, 
        command = reset_para_data)
    markers_hBut <- tkbutton(tt, image = image2, command = para_help)
    
    ## projectName
    projectName_label <- tklabel(tt, text = "Project Name :")
    projectName_entry <- tkentry(tt, textvariable = projectName, width = box_length)
    projectName_hBut <- tkbutton(tt, image = image2, command = projectName_help)
    
    ## mergeMethod && fixedNum
    mergeMethod_label <- tklabel(tt, text = "Merge Method :")
    mergeMethod_hBut <- tkbutton(tt, image = image2, command = mergeMethod_help)
    merge_method_rbuts <- tkframe(tt)
    tkpack(tklabel(merge_method_rbuts, text = ""), side = "left")
    tkpack(tkradiobutton(merge_method_rbuts, text = mergeMethods[1], 
        variable = mergeMethod, value = mergeMethods[1], command = method_all_warning), 
        side = "left")
    tkpack(tkradiobutton(merge_method_rbuts, text = mergeMethods[2], 
        variable = mergeMethod, value = mergeMethods[2], command = reset_num2null), 
        side = "left")
    tkpack(tkradiobutton(merge_method_rbuts, text = mergeMethods[3], 
        variable = mergeMethod, value = mergeMethods[3], command = reset_num2any), 
        side = "left")
    tkpack(tkradiobutton(merge_method_rbuts, text = mergeMethods[4], 
        variable = mergeMethod, value = mergeMethods[4], command = reset_num2any), 
        side = "left")
    tkpack(tkentry(merge_method_rbuts, textvariable = fixedNum, 
                   width = 9), side = "right")
    tkpack(tklabel(merge_method_rbuts, text = "Fixed Number :"), 
        side = "right")
    tkpack(tklabel(merge_method_rbuts, text = "                 "), 
        side = "left")
    fixedNum_hBut <- tkbutton(tt, image = image2, command = fixedNum_help)
    
    ## transformMethod
    transformMethod_label <- tklabel(tt, text = "Transformation Method :")
    transformMethod_hBut <- tkbutton(tt, image = image2,
        command = transformMethod_help)
    transformMethod_rbuts <- tkframe(tt)
    tkpack(tklabel(transformMethod_rbuts, text = ""), side = "left")
    tkpack(tkradiobutton(transformMethod_rbuts, text = transformMethods[1], 
        variable = transformMethod, value = transformMethods[1]), side = "left")
    tkpack(tkradiobutton(transformMethod_rbuts, text = transformMethods[2],
        variable = transformMethod, value = transformMethods[2]), side = "left")
    tkpack(tkradiobutton(transformMethod_rbuts, text = "Fixedlogicle",
        command = function(){ 
            fixedLgclParas <- fixedLogicleParameters_GUI(fixedLgclParas) 
            tclvalue(l_w) <- fixedLgclParas[1] 
            tclvalue(l_t) <- fixedLgclParas[2] 
            tclvalue(l_m) <- fixedLgclParas[3] 
            tclvalue(l_a) <- fixedLgclParas[4] 
            },
        variable = transformMethod, value = transformMethods[3]), side = "left")
    tkpack(tkradiobutton(transformMethod_rbuts, text = transformMethods[4],
        variable = transformMethod, value = transformMethods[4]), side = "left")
    
    ## cluster method
    cluster_label <- tklabel(tt, text = "Cluster Method(s) :")
    cluster_hBut <- tkbutton(tt, image = image2, command = cluster_help)
    
    clusterMethods_cbuts <- tkframe(tt)
    tkpack(tklabel(clusterMethods_cbuts, text = ""), side = "left")
    tkpack(tkcheckbutton(clusterMethods_cbuts, text = clusterMethods[1], 
                         variable = eval(clusterSelect[1])), side = "left")
    tkpack(tkcheckbutton(clusterMethods_cbuts, text = clusterMethods[2], 
                         variable = eval(clusterSelect[2])), side = "left")
    tkpack(tkcheckbutton(clusterMethods_cbuts, text = clusterMethods[3], 
                         variable = eval(clusterSelect[3]), command = highCell_warning), side = "left")
    tkpack(tkcheckbutton(clusterMethods_cbuts, text = clusterMethods[4], 
                         variable = eval(clusterSelect[4])), side = "left")
    tkpack(tkcheckbutton(clusterMethods_cbuts, text = clusterMethods[5], 
                         variable = eval(clusterSelect[5])), side = "left")
    
    ## cluster param (Rphenograph_k and FlowSOM_k)
    rphenoK_label <- tklabel(tt, text = "Rphenograph_k:")
    rphenoK_hBut <- tkbutton(tt, image = image2, command = rPk_help)
    cluster_Param <- tkframe(tt)
    tkpack(tklabel(cluster_Param, text = " "), side = "left")
    tkpack(tkentry(cluster_Param, textvariable = Rphenograph_k, width = 4), side = "left")
    tkpack(tklabel(cluster_Param, text = "                 "), side = "left")
    tkpack(tklabel(cluster_Param, text = "tsne Perplexity"), side = "left")
    tkpack(tkentry(cluster_Param, textvariable = tsne_perp, width = 4), side = "left")
    tkpack(tklabel(cluster_Param, text = "tsne Max Iterations"), side = "left")
    tkpack(tkentry(cluster_Param, textvariable = tsne_maxIter, width = 4), side = "left")
    tkpack(tkbutton(cluster_Param, image = image2, command = fSk_help), side = "right")
    tkpack(tkentry(cluster_Param, textvariable = FlowSOM_k, width = 4), side = "right")
    tkpack(tklabel(cluster_Param, text = "FlowSOM_k:"), side = "right")
    
    ## visualizationMethods
    visualizationMethods_label <- tklabel(tt, text = "Visualization Method(s) :")
    visualizationMethods_hBut <- tkbutton(tt, image = image2,
        command = visualizationMethods_help)
    visualizationMethods_cbuts <- tkframe(tt)
    tkpack(tklabel(visualizationMethods_cbuts, text = ""), side = "left")
    tkpack(tkcheckbutton(visualizationMethods_cbuts, text = vizMethods[1],
                         variable = eval(vizSelect[1])), side = "left")
    tkpack(tkcheckbutton(visualizationMethods_cbuts, text = vizMethods[2],
                         variable = eval(vizSelect[2]), command = highCell_warning), side = "left")
    tkpack(tkcheckbutton(visualizationMethods_cbuts, text = vizMethods[3],
                         variable = eval(vizSelect[3])), side = "left")
    tkpack(tkcheckbutton(visualizationMethods_cbuts, text = vizMethods[4],
                         variable = eval(vizSelect[4])), side = "left")
    tkpack(tklabel(visualizationMethods_cbuts, text = "                 "), side = "left")
    tkpack(tklabel(visualizationMethods_cbuts, text = "Seed"), side = "left")
    tkpack(tkentry(visualizationMethods_cbuts, textvariable = seed, width = 4), side = "left")
    
    ## progressionMethod
    progressionMethod_label <- tklabel(tt, text = "Cellular Progression :")
    progressionMethod_hBut <- tkbutton(tt, image = image2, 
                                     command = progressionMethod_help)
    progressionMethod_rbuts <- tkframe(tt)
    tkpack(tklabel(progressionMethod_rbuts, text = ""), side = "left")
    tkpack(tkradiobutton(progressionMethod_rbuts, text = progressionMethods[1], 
                         variable = progressionMethod, value = progressionMethods[1]), side = "left")
    tkpack(tkradiobutton(progressionMethod_rbuts, text = progressionMethods[2], 
                         variable = progressionMethod, value = progressionMethods[2], command = highCell_warning),
           side = "left")
    tkpack(tkradiobutton(progressionMethod_rbuts, text = progressionMethods[3], 
                         variable = progressionMethod, value = progressionMethods[3]), side = "left")
    
    ## submit / reset / quit
    submit_button <- tkbutton(tt, text = "Submit", command = submit)
    reset_button <- tkbutton(tt, text = "Reset", command = reset)
    quit_button <- tkbutton(tt, text = "Quit", command = quit)
    
    ## display GUI
    tkgrid(rawFCSdir_label, rawFCSdir_hBut, rawFCSdir_entry, rawFCSdir_button, 
        padx = cell_width)
    tkgrid.configure(rawFCSdir_label, rawFCSdir_entry, rawFCSdir_button, 
        sticky = "e")
    tkgrid.configure(rawFCSdir_hBut, sticky = "e")
    
    tkgrid(fcsFile_label, fcsFile_hBut, fcsFile_entry, fcsFile_button, 
        padx = cell_width)
    tkgrid.configure(fcsFile_label, fcsFile_entry, fcsFile_button, 
        sticky = "e")
    tkgrid.configure(fcsFile_hBut, sticky = "e")
    
    tkgrid(markers_label, markers_hBut, markers_entry, markers_button, 
        padx = cell_width)
    tkgrid.configure(markers_label, markers_entry, markers_button, 
        sticky = "e")
    tkgrid.configure(markers_hBut, sticky = "e")
    
    tkgrid(resDir_label, resDir_hBut, resDir_entry, resDir_button, 
        padx = cell_width)
    tkgrid.configure(resDir_label, resDir_entry, resDir_button, 
        sticky = "e")
    tkgrid.configure(resDir_hBut, sticky = "e")
    
    tkgrid(projectName_label, projectName_hBut, projectName_entry, padx = cell_width)
    tkgrid.configure(projectName_label, projectName_entry, sticky = "e")
    tkgrid.configure(projectName_hBut, sticky = "e")
    
    tkgrid(mergeMethod_label, mergeMethod_hBut, merge_method_rbuts, 
        fixedNum_hBut, padx = cell_width)
    tkgrid.configure(mergeMethod_label, sticky = "e")
    tkgrid.configure(mergeMethod_hBut, sticky = "e")
    tkgrid.configure(merge_method_rbuts, sticky = "w")
    tkgrid.configure(fixedNum_hBut, sticky = "w")
    
    tkgrid(transformMethod_label, transformMethod_hBut, transformMethod_rbuts,
        padx = cell_width)
    tkgrid.configure(transformMethod_label, sticky = "e")
    tkgrid.configure(transformMethod_rbuts, sticky = "w")
    tkgrid.configure(transformMethod_hBut, sticky = "e")
    
    tkgrid(cluster_label, cluster_hBut, clusterMethods_cbuts, padx = cell_width)
    tkgrid.configure(cluster_label, sticky = "e")
    tkgrid.configure(clusterMethods_cbuts, sticky = "w")
    tkgrid.configure(cluster_hBut, sticky = "e")
    
    tkgrid(rphenoK_label, rphenoK_hBut, cluster_Param, padx = cell_width)
    tkgrid.configure(rphenoK_label, rphenoK_hBut, sticky = "e")
    tkgrid.configure(cluster_Param, sticky = "w")
    
    tkgrid(visualizationMethods_label, visualizationMethods_hBut, 
        visualizationMethods_cbuts, padx = cell_width)
    tkgrid.configure(visualizationMethods_label, sticky = "e")
    tkgrid.configure(visualizationMethods_cbuts, sticky = "w")
    tkgrid.configure(visualizationMethods_hBut, sticky = "e")
    
    tkgrid(progressionMethod_label, progressionMethod_hBut, progressionMethod_rbuts, 
           padx = cell_width)
    tkgrid.configure(progressionMethod_label, sticky = "e")
    tkgrid.configure(progressionMethod_rbuts, sticky = "w")
    tkgrid.configure(progressionMethod_hBut, sticky = "e")
    
    tkgrid(tklabel(tt, text = "\n"), padx = cell_width)  # leave blank line
    
    tkgrid(reset_button, tklabel(tt, text = ""), submit_button, 
        quit_button, padx = cell_width)
    tkgrid.configure(reset_button, sticky = "e")
    tkgrid.configure(quit_button, sticky = "w")
    
    tkwait.window(tt)
    
    
    ##-------------------##
    ## Return parameters ##
    ##-------------------##
    
    if (tclvalue(ret_var) != "OK") {
        okMessage <- "Analysis is cancelled."
    }else{
        fcsFiles <- strsplit(tclvalue(fcsFile), "}{", fixed = TRUE)[[1]]
        parameters <- strsplit(tclvalue(markers), "}{", fixed = TRUE)[[1]]
        
        clusterCheck <- c()
        i <- 1
        while (i <= length(clusterMethods)) {
            v <- as.numeric(tclvalue(clusterSelect[i])) > 0
            clusterCheck <- c(clusterCheck, v)
            i <- i + 1
        }
        
        vizCheck <- c()
        i <- 1
        while (i <= length(vizMethods)) {
            v <- as.numeric(tclvalue(vizSelect[i])) > 0
            vizCheck <- c(vizCheck, v)
            i <- i + 1
        }
        
        inputs <- list()
        inputs[["fcsFiles"]] <- fcsFiles
        inputs[["markers"]] <- parameters
        inputs[["mergeMethod"]] <- tclvalue(mergeMethod)
        inputs[["fixedNum"]] <- suppressWarnings(as.numeric(tclvalue(fixedNum)))
        inputs[["transformMethod"]] <- tclvalue(transformMethod)
        inputs[["dimReductionMethod"]] <- "tsne"
        inputs[["clusterMethods"]] <- clusterMethods[clusterCheck]
        inputs[["visualizationMethods"]] <- vizMethods[vizCheck]
        inputs[["progressionMethod"]] <- tclvalue(progressionMethod)
        inputs[["Rphenograph_k"]] <- tclvalue(Rphenograph_k)
        inputs[["tsne_perp"]] <- tclvalue(tsne_perp)
        inputs[["tsne_maxIter"]] <- tclvalue(tsne_maxIter)
        inputs[["FlowSOM_k"]] <- tclvalue(FlowSOM_k)
        inputs[["seed"]] <- tclvalue(seed)
        inputs[["projectName"]] <- tclvalue(projectName)
        inputs[["resultDir"]] <- tclvalue(resDir)
        
        inputs[["l_w"]] <- tclvalue(l_w)
        inputs[["l_t"]] <- tclvalue(l_t)
        inputs[["l_m"]] <- tclvalue(l_m)
        inputs[["l_a"]] <- tclvalue(l_a)
        
        
        # pass the parameters and run the cytofkit function
        cytofkit(fcsFiles = inputs[["fcsFiles"]],
                 markers = inputs[["markers"]],
                 projectName = inputs[["projectName"]],
                 mergeMethod = inputs[["mergeMethod"]],
                 fixedNum = inputs[["fixedNum"]],
                 transformMethod = inputs[["transformMethod"]],
                 dimReductionMethod = inputs[["dimReductionMethod"]],
                 clusterMethods = inputs[["clusterMethods"]],
                 visualizationMethods = inputs[["visualizationMethods"]],
                 progressionMethod = inputs[["progressionMethod"]],
                 Rphenograph_k = as.numeric(inputs[["Rphenograph_k"]]),
                 FlowSOM_k = as.numeric(inputs[["FlowSOM_k"]]),
                 seed = as.numeric(inputs[["seed"]]),
                 clusterSampleSize = 500,
                 resultDir = inputs[["resultDir"]],
                 saveResults = TRUE,
                 saveObject = TRUE,
                 l_w = as.numeric(inputs[["l_w"]]), 
                 l_t = as.numeric(inputs[["l_t"]]), 
                 l_m = as.numeric(inputs[["l_m"]]), 
                 l_a = as.numeric(inputs[["l_a"]]))
        
        okMessage <- paste0("Analysis done, results are saved under ",
                            inputs[["resultDir"]])
        RData_path <- paste0(inputs[["resultDir"]], .Platform$file.sep, inputs[["projectName"]], ".RData")
    }
    
    launchShinyAPP_GUI(message = okMessage, dir = inputs[["resultDir"]], obj = RData_path)
}


#' GUI for launching shiny APP
#' 
#' A shiny APP for interactive exploration of analysis results
#' 
#' @param message Message when asking if user wants to open the shiny APP
#' @param dir Result directory.
#' @param obj The RData piped from cytofkit function to the Shiny App
#' 
#' @return Window asking user if they wish to open shinyApp directly
#' @export
#' @examples
#' # launchShinyAPP_GUI()
launchShinyAPP_GUI <- function(message="cytofkit", dir = getwd(), obj = NULL){
    
    if(message == "Analysis is cancelled."){
        message("Analysis is cancelled!")
    }else{
        ifAPP <- tclVar("n")
        ss <- tktoplevel(borderwidth = 10)
        tkwm.title(ss, "cytofkit: Analysis Done")
        
        onYes <- function() {
            tclvalue(ifAPP) <- "y"
            tkdestroy(ss)
        }
        
        onNo <- function() {
            tclvalue(ifAPP) <- "n"
            tkdestroy(ss)
        }
        yesBut <- tkbutton(ss, text = " Yes ", command = onYes)
        noBut <- tkbutton(ss, text = " No ", command = onNo)
        openDirBut <- tkbutton(ss, text = "Open", command = function(){opendir(dir)})
        okBut <- tkbutton(ss, text = "OK", command = function(){tkdestroy(ss)})
        tkgrid(tklabel(ss, text = message))
        
        tkgrid(openDirBut)
        tkgrid(tklabel(ss, text = "\n"))
        tkgrid(tklabel(ss, text = "Launch Shiny APP to check your results:"))
        tkgrid(noBut, tklabel(ss, text = "    "), yesBut)
        tkgrid.configure(noBut, sticky = "e")
        tkgrid.configure(yesBut, sticky = "e")
        tkwait.window(ss)
        
        if(tclvalue(ifAPP) == "y"){
            cytofkitShinyAPP(obj)
        }
    }
}


## function for opening the results directory
opendir <- function(dir = getwd()){
    if (.Platform['OS.type'] == "windows"){
        shell.exec(dir)
    } else {
        system(paste(Sys.getenv("R_BROWSER"), dir))
    }
}


#' GUI for marker selection 
#' 
#' Extract the markers from the fcsfiles
#' 
#' @param fcsFile The name of the FCS file
#' @param rawFCSdir The path of the FCS file
#' @return List of markers for ddimension reduction and clustering
#' @examples 
#' #getParameters_GUI()
getParameters_GUI <- function(fcsFile, rawFCSdir) {
    
    if (missing(fcsFile)) {
        fcsFile <- list.files(path = rawFCSdir, pattern = ".fcs$", full.names = TRUE)
    }
    
    fcs <- suppressWarnings(read.FCS(fcsFile[1]))
    pd <- fcs@parameters@data
    markers <- paste(pd$name, "<", pd$desc, ">", sep = "")
    channels <- paste(pd$name, "<", pd$desc, ">", sep = "")
    
    if (length(markers) == 0) {
        stop("No markers found in the FCS file!")
    }
    
    # GUI
    markerChoice <- tclVar("")
    mm <- tktoplevel()
    tkwm.title(mm, "cytofkit: Marker Selection")
    scr <- tkscrollbar(mm, repeatinterval = 5, command = function(...) tkyview(tl, 
        ...))
    tl <- tklistbox(mm, height = 30, width = 40, selectmode = "multiple", yscrollcommand = function(...) tkset(scr, 
        ...), background = "white")
    OnOK <- function() {
        tclvalue(markerChoice) <- paste(markers[as.numeric(tkcurselection(tl)) + 
            1], collapse = "}{")
        tkdestroy(mm)
    }
    OK.but <- tkbutton(mm, text = " OK ", command = OnOK)
    tkgrid(tklabel(mm, text = "Please select your markers:"))
    tkgrid(tl, scr)
    tkgrid.configure(scr, rowspan = 4, sticky = "nsw")
    for (i in (1:length(markers))) {
        tkinsert(tl, "end", markers[i])
    }
    tkgrid(OK.but)
    tkwait.window(mm)
    
    # return parameters
    paras <- strsplit(tclvalue(markerChoice), "}{", fixed = TRUE)[[1]]
    paras <- channels[match(paras, markers)]
    return(paras)
} 


#' GUI for gettting parameter for logicle transformation
#' 
#' Extract the parameter for fixed logicle transformation
#' 
#' @param fixedLgclParas parameters vector containing w, t, m, a
#' @return Parameters for fixed logicle transformation
#' @examples 
#' #fixedLogicleParameters_GUI
fixedLogicleParameters_GUI <- function(fixedLgclParas=c(0.5, 500000, 4.5, 0)) {
    
    # logicle parameters
    l_w <- tclVar(fixedLgclParas[1])
    l_t <- tclVar(fixedLgclParas[2])
    l_m <- tclVar(fixedLgclParas[3])
    l_a <- tclVar(fixedLgclParas[4])
    
    cell_width <- 3
    box_length <- 10
    
    # GUI
    mm <- tktoplevel(borderwidth = 20)
    tkwm.title(mm, "Parameters for fixed logicle transformation")
    
    w_label <- tklabel(mm, text = "w :")
    w_entry <- tkentry(mm, textvariable = l_w, width = box_length)
    
    t_label <- tklabel(mm, text = "t :")
    t_entry <- tkentry(mm, textvariable = l_t, width = box_length)
    
    m_label <- tklabel(mm, text = "m :")
    m_entry <- tkentry(mm, textvariable = l_m, width = box_length)
    
    a_label <- tklabel(mm, text = "a :")
    a_entry <- tkentry(mm, textvariable = l_a, width = box_length)
    
    OnOK <- function(){
        tkdestroy(mm)
    }
    
    OK.but <- tkbutton(mm, text = " OK ", command = OnOK)
    
    tkgrid(tklabel(mm, text = "Specify your parameters\nfor logical transformation:"), columnspan=2)
    
    tkgrid(tklabel(mm, text = "\n"), padx = cell_width)  # leave blank line
    
    tkgrid(w_label, w_entry, padx = cell_width)
    tkgrid.configure(w_label, sticky = "e")
    tkgrid.configure(w_entry, sticky = "w")
    
    tkgrid(t_label, t_entry, padx = cell_width)
    tkgrid.configure(t_label, sticky = "e")
    tkgrid.configure(t_entry, sticky = "w")
    
    tkgrid(m_label, m_entry, padx = cell_width)
    tkgrid.configure(m_label, sticky = "e")
    tkgrid.configure(m_entry, sticky = "w")
    
    tkgrid(a_label, a_entry, padx = cell_width)
    tkgrid.configure(a_label, sticky = "e")
    tkgrid.configure(a_entry, sticky = "w")

    tkgrid(tklabel(mm, text = "\n"), padx = cell_width)  # leave blank line
    
    tkgrid(tklabel(mm, text = ""), OK.but, padx = cell_width)
    tkgrid.configure(OK.but, sticky = "w")
    tkwait.window(mm)
    
    # return parameters
    paras <- c(tclvalue(l_w), tclvalue(l_t), tclvalue(l_m), tclvalue(l_a))
    paras <- as.numeric(paras)
    return(paras)
} 

#' A Shiny APP to interactively visualize the analysis results 
#' 
#' Take the the RData object file saved by cytofkit as input, automatically load the data and allow exploration of the analysis results with interactive control
#'
#'
#' @param RData Either the RData object file or data object, if missing, RData file need to be loaded on the ShinyAPP
#' @param onServer Logical value, if \verb{TRUE}, sets shinyApp host to 0.0.0.0 for other clients to access, otherwise defaults to 127.0.0.1 (local host)
#' 
#' @return Opens shinyApp session for data visualisation
#' @import shiny
#' @import shinyFiles
#' @importFrom grDevices dev.copy2pdf
#' @importFrom graphics plot
#' @author Hao Chen
#' @export
#' @examples 
#' d <- system.file('extdata', package = 'cytofkit2')
#' Rdata <- list.files(d, pattern = '.RData$', full.names = TRUE)
#' #only for interactive sessions, remove hash to run
#' #cytofkitShinyAPP(Rdata)
cytofkitShinyAPP <- function(RData = NULL, onServer = FALSE) {
    
    source(system.file('shiny', "global.R", package = 'cytofkit2'))
  
    analysis_results <- NULL
    sampleInformation <- NULL
    progCluster <- NULL
    serverObj <- NULL
    roots <- c(wd=getwd())
    
    if(!missing(RData)){
        if(class(RData) == "character"){
            if(file.exists(RData)){
              if(tools::file_ext(RData) == "RData"){
                load(RData)
                direct_analysis_results <- analysis_results
                message(".RData loaded!")
              }else{
                stop("Argument is not .RData file!")
              }
            }else{
                stop("RData file doesn't exist! Please check your obj file")
            }
        }else{
            analysis_results <- RData
        }
        
        if(is.null(analysis_results$projectName)){
            analysis_results$projectName <- "cytofkit_shinyAPP_output"
        }
        
        if(!is.null(analysis_results$progressionRes)){
            ## default the first cluster results are used for progression analysis
            progCluster <- names(analysis_results$clusterRes)[1]
        }
        
        sampleInformation <- data.frame(cellID = row.names(analysis_results$expressionData),
                                        cellSample = factor(sub("_[0-9]*$", "", row.names(analysis_results$expressionData))),
                                        stringsAsFactors = FALSE)
        analysis_results$sampleInfo <- sampleInformation
    }
    
    if(isTRUE(onServer)){
      host <- "0.0.0.0"
    }else{
      host <- "127.0.0.1"
    }
    
    #shiny::runApp(system.file('shiny', package = 'cytofkit2'))
    options(shiny.launch.browser = TRUE, shiny.port = 4455, shiny.host = host, shiny.maxRequestSize=1024^10)
    shinyApp(
        ui = fluidPage(
          titlePanel("Interactive Exploration of cytofkit Analysis Results"),
          hr(),
          fluidRow(
            ## side panel--take 1/4 space
            column(3,
                   h4('Load cytofkit RData:'),
                   wellPanel(
                     fileInput(inputId = 'cytofkitObj',
                               label = NULL,
                               multiple = FALSE,
                               accept = c('text/RData', '.RData')),
                     shinyFilesButton('serverObj', label = "Server File Select", title = "Please select your RData", multiple = FALSE),
                     textOutput("rdata_desc"),
                     fluidRow(
                       column(6,
                              actionButton("goButton", "Submit", icon = icon("hand-o-right"))
                       ),
                       column(6,
                              actionButton("reset", "Reset Data", icon = icon("repeat"))
                       )
                     )
                   ),
                   
                   hr(),
                   
                   conditionalPanel(" input.main_panel == 'C_panel' && input.C_clusterTabs == 'C_tab1' ",
                                    h4("Plot Control:"),
                                    wellPanel(
                                      checkboxInput(inputId = "C_addLabel", label = "Add Cluster Labels", value = TRUE),
                                      checkboxInput(inputId = "C_labelRepel", label = "Repel Cluster Labels", value = FALSE),
                                      checkboxInput(inputId = "C_facetPlot", label = "Separate Plot by Samples", value = FALSE)
                                    ),
                                    actionButton("PDFClusterPlot", "Download Cluster Plot in PDF", icon = icon("download"))
                                    ),
                   conditionalPanel(" input.main_panel == 'M_panel' && input.M_markerTabs == 'M_tab1' ",
                                    h4("Plot Control:"),
                                    wellPanel(
                                      selectInput('M_heatmap_dendrogram', strong('Heatmap Dendrogram:'), 
                                                  choices = c("both","row","column","none"), 
                                                  selected = "both", width = "100%"),  
                                      selectInput('M_heatmap_colorPalette', strong('Color Palette:'), 
                                                  choices = c("bluered", "greenred", "spectral1", "spectral2"), 
                                                  selected = "bluered", width = "100%")
                                    ),
                                    actionButton("PDFHeatmap", "Download Marker Heatmap in PDF", icon = icon("download"))
                                    ),
                   conditionalPanel(" input.main_panel == 'M_panel' && input.M_markerTabs == 'M_tab2' ",
                                    h4("Plot Control:"),
                                    wellPanel(
                                      actionButton("PDFExpPlot", "Download Exp Plot in PDF", icon = icon("download"))
                                    )),
                   conditionalPanel(" input.main_panel == 'M_panel' && input.M_markerTabs == 'M_tab3' ",
                                    h4("Plot Control:"),
                                    wellPanel(
                                      actionButton("PDFHistogram", "Download Histogram in PDF", icon = icon("download"))
                                    )),
                   conditionalPanel(" input.main_panel == 'S_panel' && input.S_sampleTabs == 'S_tab1' ",
                                    h4("Plot Control:"),
                                    wellPanel(
                                      selectInput('S_heatmap_dendrogram', strong('Heatmap Dendrogram:'), 
                                                  choices = c("both","row","column","none"), 
                                                  selected = "both", width = "100%"),  
                                      selectInput('S_heatmap_colorPalette', strong('Color Palette:'), 
                                                  choices = c("bluered", "greenred", "spectral1", "spectral2"), 
                                                  selected = "bluered", width = "100%")
                                    ),
                                    actionButton("PDFSamHeat", "Download Sample Heatmap in PDF", icon = icon("download"))
                                    ),
                   conditionalPanel(" input.main_panel == 'S_panel' && input.S_sampleTabs == 'S_tab2' ",
                                    h4("Plot Control:"),
                                    actionButton("PDFrateChange", "Download Rate Change Plot in PDF", icon = icon("download"))
                                    ),
                   conditionalPanel(" input.main_panel == 'P_panel' && input.P_progressionTabs == 'P_tab1' ",
                                    h4("Plot Control:"),
                                    wellPanel(
                                      checkboxInput(inputId = "P_addLabel", label = "Add Cluster Labels", value = TRUE),
                                      checkboxInput(inputId = "P_labelRepel", label = "Repel Cluster Labels", value = FALSE),
                                      checkboxInput(inputId = "P_facetPlot", label = "Separate Plot by Samples", value = FALSE)
                                    ),
                                    actionButton("PDFScatter", "Download Scatterplot in PDF", icon = icon("download"))
                                    ),
                   conditionalPanel(" input.main_panel == 'P_panel' && input.P_progressionTabs == 'P_tab2' ",
                                    h4("Plot Control:"),
                                    wellPanel(
                                      checkboxInput(inputId = "P_addLabel2", label = "Add Cluster Labels", value = TRUE)
                                    ),
                                    actionButton("PDFmarkerPlot", "Download Marker Plot in PDF", icon = icon("download"))
                                    ),
                   br(),
                   fluidRow(
                     column(6,
                            sliderInput(inputId="tab_w", label = "PDF width(in):", 
                                        min=3, max=20, value=8, width=100, ticks=FALSE)
                     ),
                     column(6, 
                            sliderInput(inputId="tab_h", label = "PDF height(in):", 
                                        min=3, max=20, value=8, width=100, ticks=FALSE)
                     )),
                   
                   actionButton("OpenDir", "Open download folder", icon = icon("folder")),
                   
                   hr(),
                   h4("Sample Filter:"),
                   wellPanel(uiOutput("selectAll"),
                             uiOutput("sampleSelect")),
                   
                   hr(),
                   h4("Data Summary:"),
                   wellPanel(
                     h5("Expression Data:"),
                     textOutput("summaryText1"),
                     h5("Markers used for dimension reduction and clustering:"),
                     textOutput("summaryText5"),
                     h5("Cluster Method(s):"),
                     textOutput("summaryText2"),
                     h5("Visualization Method(s):"),
                     textOutput("summaryText3"),
                     h5("Progression Method(s):"),
                     textOutput("summaryText4")
                   ),
                   
                   hr(),
                   h4("Save results:"),
                   h5("Outputs to save"),
                   fluidRow(
                     column(4,
                            checkboxInput(inputId = "saveFCS", label = "FCS", value = TRUE)
                     ),
                     column(4,
                            checkboxInput(inputId = "saveRData", label = "RData", value = TRUE)
                     ),
                     column(4,
                            checkboxInput(inputId = "saveCsv", label = "csv", value = FALSE)
                     )
                   ),
                   actionButton("saveButton", "Save Data", icon = icon("download")),
                   
                   hr(),
                   h4(tags$a(href="mailto:jinmiao@gmail.com,a0124008@u.nus.edu?subject=[cytofkit-question]", 
                             "Contact Us")),
                   imageOutput("logo", height = "60px")
            ),
            ## main panel--take 3/4 space
            column(9,
                   tabsetPanel(id="main_panel", type = "pills",
                               tabPanel(title="Cluster Panel", value="C_panel", fluidPage(
                                 hr(),
                                 tabsetPanel(id="C_clusterTabs", type = "tabs",
                                             tabPanel(title="Cluster Plot", value="C_tab1", 
                                                      br(),
                                                      fluidRow(
                                                        column(3,
                                                               uiOutput("C_PlotMethod")
                                                        ),
                                                        column(3, 
                                                               uiOutput("C_PlotFunction")
                                                        ),
                                                        column(3,
                                                               numericInput("C_PointSize", "Point Size:", value = 1)
                                                        ),
                                                        column(3, 
                                                               numericInput("C_LabelSize", "Label Size:", value = 12)
                                                        )
                                                      ),
                                                      uiOutput("C_clusterSelect"),
                                                      hr(),
                                                      plotOutput("C_ScatterPlot", width = "100%")
                                             ),
                                             tabPanel(title="Change Cluster Color", value="C_tab2",
                                                      br(),
                                                      wellPanel(
                                                        uiOutput("C_colourCluster")
                                                      ),
                                                      hr(),
                                                      lapply(1:100, function(i) {
                                                        uiOutput(paste0('Cluster_', i, '_col'))
                                                      }),
                                                      hr(),
                                                      fluidRow(
                                                        column(3,
                                                               actionButton("C_updateClusterColor", "Update Cluster Color", 
                                                                            icon = icon("hand-o-right"), width = "100%")
                                                        ),
                                                        column(3, 
                                                               actionButton("C_revertClusterColor", "Revert to default", 
                                                                            icon = icon("hand-o-right"), width = "100%")
                                                        ),
                                                        column(6)
                                                      ),
                                                      hr()),
                                             tabPanel(title="Annotate Clusters", value="C_tab3",
                                                      br(),
                                                      wellPanel(
                                                        uiOutput("C_labelCluster"),
                                                        uiOutput("C_labelCluster_name")
                                                      ),
                                                      hr(),
                                                      lapply(1:100, function(i) {
                                                        uiOutput(paste0('Cluster', i))
                                                      }),
                                                      hr(),
                                                      actionButton("updatelabel", "Submit Cluster Label", icon = icon("hand-o-right")),
                                                      hr()),
                                             tabPanel(title="Run FlowSOM", value="C_tab4",
                                                      br(),
                                                      h4("FlowSOM Clustering Setup:"),
                                                      hr(),
                                                      wellPanel(
                                                        numericInput("C_FlowSOM_k", "Cluster k", value = 10, width = "30%"),
                                                        uiOutput("C_markerSelect")
                                                      ),
                                                      hr(),
                                                      actionButton("C_runFlowSOM", "Run FlowSOM", icon = icon("hand-pointer-o")))
                                 )
                               )),
                               
                               tabPanel(title = "Marker Panel", value = "M_panel", fluidPage(
                                 hr(),
                                 
                                 tabsetPanel(id="M_markerTabs", type = "tabs",
                                             tabPanel(title="Expression Heat Map", value="M_tab1", 
                                                      br(),
                                                      fluidRow(
                                                        column(4, 
                                                               uiOutput("M_plotCluster")
                                                        ),
                                                        column(2,
                                                               selectInput('M_plotMethod', strong('Heatmap Type:'), 
                                                                           choices = c("mean", "median"), 
                                                                           selected = "mean", width = "100%")
                                                        ),
                                                        column(2,
                                                               selectInput('M_scaleMethod', strong('Scale Data:'), 
                                                                           choices = c("none", "row", "column"), 
                                                                           selected = "none", width = "100%")
                                                        ),
                                                        column(2,
                                                               numericInput("M_rowLabelSize", "Row Label Size:", value = 1, step = 0.5)
                                                        ),
                                                        column(2, 
                                                               numericInput("M_colLabelSize", "Col Label Size:", value = 1, step = 0.5)
                                                        )
                                                      ),
                                                      fluidRow(
                                                        column(10,
                                                               uiOutput("M_heatmapmarkerSelect")
                                                        ),
                                                        column(2,
                                                               actionButton("M_heatmapSelectAll", "All Markers"),
                                                               actionButton("M_updateHeatmap", "Update Plot")
                                                        )
                                                      ),
                                                      hr(),
                                                      plotOutput("M_heatmapPlot", width = "100%")),
                                             tabPanel(title="Expression Level Plot", value="M_tab2",
                                                      br(),
                                                      fluidRow(
                                                        column(3,
                                                               uiOutput("M_PlotMethod")
                                                        ),
                                                        column(3,
                                                               numericInput("M_PointSize", "Point Size:", value = 1),
                                                               sliderInput("M_Alpha", "Transparency:", value = 1, min = 0, max = 1, step = 0.1)
                                                        ),
                                                        column(3,
                                                               selectInput('M_colorPalette', label = "Color Palette:", 
                                                                           choices = c("bluered", "spectral1", "spectral2", "heat"), 
                                                                           selected = "bluered", width = "100%")
                                                        ),
                                                        column(3,
                                                               selectInput('M_ScaleOptions', label = "Scaling Range:", 
                                                                           choices = c("Local", "Global"), 
                                                                           selected = "Local", width = "100%"),
                                                               selectInput('M_scaledData', label = "Centering:", 
                                                                           choices = c("Un-centered", "Centered"), 
                                                                           selected = "Un-centered", width = "100%")
                                                        )
                                                      ),
                                                      fluidRow(
                                                        column(10,
                                                               uiOutput("M_PlotMarker")
                                                        ),
                                                        column(2,
                                                               actionButton("M_chooseAllMarker", "All Markers"),
                                                               actionButton("M_updateExPlot", "Update Plot")
                                                        )
                                                      ),
                                                      hr(),
                                                      plotOutput("M_markerExpressionPlot", width = "100%")), 
                                             tabPanel(title="Expression Histogram", value="M_tab3", 
                                                      br(),
                                                      fluidRow(
                                                        column(4,
                                                               uiOutput("M_stackFactor")
                                                        ),
                                                        column(2,
                                                               numericInput("M_markerTextSize", "Marker Text Size:", 
                                                                            value = 12, step = 1, min=1, max=15)
                                                        ),
                                                        column(2,
                                                               numericInput("M_xlab_size", "x Label Size:", 
                                                                            value = 2, step = 1, min=1, max=10)
                                                        ),
                                                        column(2,
                                                               numericInput("M_legendTextSize", "Legend Size:", 
                                                                            value = 1, step = 0.5, min=1, max=10)
                                                        ),
                                                        column(2,
                                                               numericInput("M_legendRow", "Legend Row:", 
                                                                            value = 2, step = 1, min=1, max=10)
                                                        )
                                                      ),
                                                      fluidRow(
                                                        column(10,
                                                               uiOutput("M_markerSelect")
                                                        ),
                                                        column(2,
                                                               actionButton("M_histSelectAll", "All Markers")
                                                        )
                                                      ),
                                                      hr(),
                                                      actionButton("M_updateDensityPlot", "Update Plot", icon = icon("hand-pointer-o")),
                                                      plotOutput("M_stackDensityPlot", width = "100%")),
                                             tabPanel(title="Update Marker Names", value="M_tab4", 
                                                      h5('Type in Your New Name for Each Marker:'),
                                                      hr(),
                                                      lapply(1:100, function(i) {
                                                        uiOutput(paste0('Marker_', i, "_name"))
                                                      }),
                                                      hr(),
                                                      actionButton("C_updateMarkerNames", "Update Marker Name", icon = icon("hand-pointer-o")))
                                 )
                               )),
                               
                               tabPanel(title = "Sample Panel", value = "S_panel", fluidPage(
                                 hr(),
                                 tabsetPanel(id="S_sampleTabs", type = "tabs",
                                             tabPanel(title="Cell Percentage Heatmap", value="S_tab1", 
                                                      br(),
                                                      fluidRow(
                                                        column(4, 
                                                               uiOutput("S_plotCluster")
                                                        ),
                                                        column(2,
                                                               selectInput('S_plotMethod', strong('Heatmap Type:'), 
                                                                           choices = c("percentage"), 
                                                                           selected = "percentage", width = "100%")
                                                        ),
                                                        column(2,
                                                               selectInput('S_scaleMethod', strong('Scale Data:'), 
                                                                           choices = c("none", "row", "column"), 
                                                                           selected = "none", width = "100%")
                                                        ),
                                                        column(2,
                                                               numericInput("S_rowLabelSize", "Row Label Size:", value = 1, step = 0.5)
                                                        ),
                                                        column(2, 
                                                               numericInput("S_colLabelSize", "Col Label Size:", value = 1, step = 0.5)
                                                        )
                                                      ),
                                                      hr(),
                                                      plotOutput("S_heatmapPlot", width = "100%")
                                             ),
                                             tabPanel(title="Cell Percentage Line Chart", value="S_tab2", 
                                                      br(),
                                                      uiOutput("S_clusterMethod2"),
                                                      uiOutput("S_clusterFilter"),
                                                      hr(),
                                                      plotOutput("S_rateChangePlot", width = "100%")
                                             ),
                                             tabPanel(title="Regroup Samples", value="S_tab3",
                                                      br(),
                                                      h4("Type in the Group Name for Each Sample:"),
                                                      lapply(1:100, function(i) {
                                                        uiOutput(paste0('S_sample', i))
                                                      }),
                                                      hr(),
                                                      textInput("sampleGroupLevels", "Group Name Levels: (to order the group names)", 
                                                                value = "", width = "100%",
                                                                placeholder = "Type in group names in order, seperated by semicolon(;)"),
                                                      hr(),
                                                      fluidRow(
                                                        column(3,
                                                               actionButton("updateSampleGroups", "Submit New Sample Groups", icon = icon("hand-o-right"))
                                                        ),
                                                        column(3, 
                                                               actionButton("revertSampleNames", "Revert to Old Sample Names", icon = icon("hand-o-right"))
                                                        ),
                                                        column(6)
                                                      ),
                                                      hr())
                                 )
                               )),
                               
                               tabPanel(title="Progression Panel", value = "P_panel", fluidPage(
                                 hr(),
                                 tabsetPanel(id="P_progressionTabs", type = "tabs",
                                             tabPanel(title="Subset Relationship Plot", value="P_tab1", 
                                                      br(),
                                                      fluidRow(
                                                        column(3,
                                                               uiOutput("P_xlab")
                                                        ),
                                                        column(3,
                                                               uiOutput("P_ylab")
                                                        ),
                                                        column(3,
                                                               numericInput("P_PointSize", "Point Size:", value = 3)
                                                        ),
                                                        column(3, 
                                                               numericInput("P_LabelSize", "Label Size:", value = 8)
                                                        )
                                                      ),
                                                      plotOutput("P_scatterPlot", width = "80%")), 
                                             tabPanel(title="Marker Expression Profile", value="P_tab2", 
                                                      br(),
                                                      fluidRow(
                                                        column(3,
                                                               uiOutput("P_orderBy")
                                                        ),
                                                        column(2,
                                                               numericInput("P_LabelSize2", "Label Size:", value = 5)
                                                        ),
                                                        column(7,
                                                               uiOutput("P_clusterSelect")
                                                        )
                                                      ),
                                                      hr(),
                                                      uiOutput("P_markerSelect"),
                                                      hr(),
                                                      fluidRow(
                                                        column(3,
                                                               actionButton("P_updateRegressionPlot", "Update Plot", icon = icon("hand-pointer-o"))
                                                        ),
                                                        column(2,
                                                               checkboxInput("P_reverseOrder", label = "Reverse Order", value = FALSE)
                                                        ),
                                                        column(3,
                                                               checkboxInput("P_combineTrends", label = "Combine Trend Lines", value = FALSE)
                                                        ),
                                                        column(4)
                                                        
                                                      ),
                                                      plotOutput("P_markerPlot", width = "100%")), 
                                             tabPanel(title="Run Diffusion Map", value="P_tab3",
                                                      br(),
                                                      h4("Diffusionmap Setup:"),
                                                      
                                                      wellPanel(
                                                        h5("Cluster-based down-sampling to remove subset aboundance heterogeneity"),
                                                        
                                                        fluidRow(
                                                          column(4,
                                                                 uiOutput("P_clusterMethod")
                                                          ),
                                                          column(4,
                                                                 numericInput("P_clusterSampleSize", "Cluster Sample Size", value = 500, 
                                                                              min = 10, max = 1000, step = 5, width = "100%")
                                                          ),
                                                          column(4,
                                                                 selectInput('P_sampleMethod', 'Downsample Method:', choices = c("ceil", "all", "fixed", "min"), 
                                                                             selected = "ceil", width = "100%")
                                                          )
                                                        ),
                                                        
                                                        tableOutput('P_clusterTable'),
                                                        
                                                        uiOutput("P_clusterFilter"),
                                                        hr(),
                                                        
                                                        h5("Diffusionmap Parameters"),
                                                        fluidRow(
                                                          column(6,
                                                                 selectInput('P_distMethod', 'Distance calculation Method:', choices = c("euclidean"), 
                                                                             selected = "euclidean", width = "100%")
                                                          ),
                                                          column(6,
                                                                 numericInput("P_outDim", "Output Dimensionality:", value = 4, 
                                                                              min = 1, max = 6, step = 1, width = "100%")
                                                          )
                                                        )
                                                      ),
                                                      hr(),
                                                      actionButton("P_runDiffusionmap", "Run Diffusionmap", icon = icon("hand-pointer-o")))
                                 )
                               )) 
                   )
            )
          )
        ),
        
        server = function(input, output, session) {
          
          ##------------------Reactive Values and Reactive Objects-------------------
          
          #if?
          v <- reactiveValues(data = NULL, sampleInfo = NULL)
          c <- reactiveValues(clusterCol = list())
          p <- reactiveValues(progressionCluster = NULL)
          
          if(!is.null(analysis_results)) {
            v$data <- analysis_results
            v$sampleInfo <- data.frame(cellID = row.names(analysis_results$expressionData),
                                       cellSample = factor(sub("_[0-9]*$", "", row.names(analysis_results$expressionData))),
                                       stringsAsFactors = FALSE)
            p$progressionCluster <- names(analysis_results$clusterRes)[1]
          }
          
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
          
          ## Save and parse cytofkit RData object
          observeEvent(input$saveButton, {
            if (!is.null(v$data)){
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
              }
          })
          
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
                colourInput(inputId=paste0('cluster_', i, '_col'), 
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
            dev.copy2pdf(file = "cytofkit_shinyAPP_marker_heatmap.pdf",
                         width=as.integer(input$tab_w), 
                         height=as.integer(input$tab_h))
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
        }
    )
    
    
    
}

