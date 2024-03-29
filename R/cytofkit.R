#' cytofkit2: an integrated mass cytometry data analysis pipeline
#' 
#' This package is designed to facilitate the analysis workflow of mass cytometry data with automatic subset identification and mapping of cellular progression. Both command line and a GUI client are provided for executing the workflow easily.
#' 
#' This package integrates merging methods of multiple FCS files, dimension reduction methods (PCA, t-SNE and ISOMAP) and clustering methods (DensVM, ClusterX, and Rphenograph) for rapid subset detection. Analysis results can be visualized and explored interactively using a specially designed shiny web APP, see \code{\link{cytofkitShinyAPP}}. Moreover, the method isomap is provided to map the cellular progression. 
#' This workflow can be easily executed with the main function.
#' 
#' Pre-processing
#' 
#' Using function \code{\link{cytof_exprsMerge}}, one or multiple FCS files will be loaded via the *read.FCS* function in the *flowCore* package. Then transformation will be applied to the expression value of each FCS file. Transformation methods include \code{autoLgcl}, \code{cytofAsinh}, \code{logicle} and \code{arcsinh}, where \code{cytofAsinh} is the default.Then mutilple FCS files are merged using one of the merging methods \code{all}, \code{min}, \code{fixed} or \code{ceil}.
#' 
#' Dimensionality reduction
#' 
#' Using function \code{\link{cytof_dimReduction}}, t-Distributed Stochastic Neighbor Embedding (\code{tsne}) is suggested for dimensionality reduction with selected markers, although we also provide methods like \code{isomap} and \code{pca}.
#' 
#' Cluster 
#' 
#' Using function \code{\link{cytof_cluster}}, three cluster method are provided, \code{DensVM}, \code{ClusterX}, \code{Rphenograph} and \code{FlowSOM}. \code{DensVM}, \code{densityClustX} are performend on the dimension reduced data, while \code{Rphenograph} works directed on the high dimensional expression data. Method \code{FlowSOM} is integrated from FlowSOM package (https://bioconductor.org/packages/release/bioc/html/FlowSOM.html). 
#' 
#' Post-processing
#' 
#'  - Using function \code{\link{cytof_clusterPlot}} to visualize the cluster results in a scatter plot, in which dots represent cells, colours indicate their assigned clusters and point shapes represent their belonging samples.
#' 
#'  - Using function \code{\link{cytof_heatmap}} to generate heat map to visualize the mean expression of every marker in every cluster. This heat maps is useful to interrogate marker expression to identify each cluster's defining markers. 
#' 
#'  - Using function \code{\link{cytof_progressionPlot}} to visualize the expression pattern of selected markers against the estimated cellular progression order.
#'  
#'  - Using function \code{\link{cytof_addToFCS}} to add any dimension reduced data, cluster results, progression data into the original FCS files, new FCS files will be saved for easy checking with other softwares like FlowJo. 
#' 
#' All the above post processing can be automatically implemented and saved using one function \code{\link{cytof_writeResults}}.
#' 
#' @author Hao Chen, Jinmiao Chen
#' 
#' 
#' ## Run on shiny_dashboard
#' #cytofkit_shiny_dashboard()  # remove the hash symbol to launch the GUI
#' 
#' ## Run on command
#' dir <- system.file('extdata',package='cytofkit2')
#' file <- list.files(dir, pattern='.fcs$', full=TRUE)
#' parameters <- list.files(dir, pattern='.txt$', full=TRUE)
#' ## remove the hash symbol to run the following command
#' #results <- cytofkit(fcsFile = file, markers = parameters)   
#' 
#' ## Checking the vignettes for more details 
#' if(interactive()) browseVignettes(package = 'cytofkit2')
#' 
#' @references \url{https://github.com/JinmiaoChenLab/cytofkit2}
#' @docType package
#' @name cytofkit2-package
#' 
NULL



#' cytofkit2: an integrated mass cytometry data analysis pipeline
#' 
#' The main function to drive the cytofkit2 workflow.
#' 
#' cytofkit2 works as the main funciton to perform the analysis of one or multiple FCS files. 
#' The workflow contains data merging from multiple FCS file, expression data transformation, 
#' dimensionality reduction with \code{PCA}, \code{isomap} or \code{tsne} (default), clustering 
#' analysis with methods includes \code{DensVM}, \code{ClusterX}, \code{Rphenograph)} and \code{FlowSOM} for 
#' subpopulation detection, and estimation of cellular progression using \code{isomap}. The analysis 
#' results can be visualized using scatter plot, heatmap plot or progression plot. Dimension reduced 
#' data and cluster labels will be saved back to new copies of FCS files. By default the analysis 
#' results will be automatically saved under \code{resultDir} for further annotation. Moreover An 
#' interactive web application is provided for interactive exploration of the analysis results, 
#' see \code{cytofkitShinyAPP}.
#' 
#' 
#' @param fcsFiles It can be either the path where your FCS files are stored or a vector of FCS file names. 
#' @param markers It can be either a text file that containing markers to be used for analysis or a vector of the marker names.
#' @param projectName A prefix that will be added to the names of all result files.
#' @param ifCompensation Boolean value, to apply compensation contained in FCS, or a compensation matrix.
#' @param transformMethod Data Transformation method, including \code{autoLgcl}, \code{cytofAsinh}, \code{logicle} and \code{arcsinh}, or \code{none} to avoid transformation.
#' @param mergeMethod When multiple fcs files are selected, cells can be combined using 
#' one of the four different methods including \code{ceil}, \code{all}, \code{min}, \code{fixed}. 
#' The default option is \code{ceil}, up to a fixed number (specified by \code{fixedNum}) of cells are sampled 
#' without replacement from each fcs file and combined for analysis.
#' \code{all}: all cells from each fcs file are combined for analysis. 
#' \code{min}: The minimum number of cells among all the selected fcs files are sampled from each fcs file and combined for analysis. 
#' \code{fixed}: a fixed num (specified by fixedNum) of cells are sampled (with replacement when the total number of cell is less than 
#' fixedNum) from each fcs file and combined for analysis.
#' @param fixedNum The fixed number of cells to be extracted from each FCS file.
#' @param dimReductionMethod The method used for dimensionality reduction, including \code{tsne}, \code{pca} and \code{isomap}.
#' @param umap_neighbor This parameter controls how UMAP balances local versus global structure in the data.
#' @param umap_min_dist Controls how tightly UMAP is allowed to pack points together.
#' @param clusterMethods The clustering method(s) used for subpopulation detection, including \code{DensVM}, \code{ClusterX}, \code{Rphenograph} and \code{FlowSOM}. Multiple selections are accepted.
#' @param visualizationMethods The method(s) used for visualize the cluster data, including \code{tsne}, \code{pca} and \code{isomap}. Multiple selections are accepted.
#' @param progressionMethod Use the first ordination score of \code{isomap} to estimated the progression order of cells, choose \code{NULL} to ignore.
#' @param Rphenograph_k Integer number of nearest neighbours to pass to Rphenograph.
#' @param FlowSOM_k Number of clusters for meta clustering in FlowSOM.
#' @param seed Integer to set a seed for reproducible results.
#' @param clusterSampleSize The uniform size of each cluster.
#' @param resultDir The directory where result files will be generated.
#' @param saveResults Save the results, and the post-processing results including scatter plot, heatmap, and statistical results.
#' @param saveObject Save the results into RData objects for loading back to R for further analysis
#' @param openShinyAPP Opens the shinyAPP automatically when the analysis was done, default FALSE.
#' @param ... Other arguments passed to \code{cytof_exprsExtract}
#' 
#' @return a list containing \code{expressionData}, \code{dimReductionMethod}, \code{visualizationMethods}, \code{dimReducedRes}, \code{clusterRes}, \code{progressionRes}, \code{projectName}, \code{rawFCSdir} and \code{resultDir}. If choose 'saveResults = TRUE', results will be saved into files under \code{resultDir}.
#' @author Hao Chen, Jinmiao Chen
#' @references \url{https://github.com/JinmiaoChenLab/cytofkit2}
#' @useDynLib cytofkit2
#' @export
#' @examples
#' dir <- system.file('extdata',package='cytofkit2')
#' file <- list.files(dir, pattern='.fcs$', full=TRUE)
#' parameters <- list.files(dir, pattern='.txt$', full=TRUE)
#' ## remove the hash symbol to run the following command
#' #cytofkit(fcsFile = file, markers = parameters) 
cytofkit <- function(fcsFiles = getwd(), 
                     markers = "parameter.txt", 
                     projectName = "cytofkit", 
                     ifCompensation = FALSE, 
                     transformMethod = c("autoLgcl", "cytofAsinh", "logicle", "arcsinh", "none"), 
                     mergeMethod = c("ceil", "all", "min", "fixed"), 
                     fixedNum = 10000, 
                     dimReductionMethod = c("umap", "tsne", "pca", "isomap"), 
                     umap_neighbor = 30,
                     umap_min_dist = 0.3,
                     clusterMethods = c("Rphenograph", "ClusterX", "DensVM", "FlowSOM", "NULL"), 
                     visualizationMethods = c("umap", "tsne", "pca", "isomap", "NULL"), 
                     progressionMethod = c("NULL", "diffusionmap", "isomap"),
                     Rphenograph_k = 30,
                     FlowSOM_k = 40,
                     seed = NULL,
                     clusterSampleSize = 500,
                     resultDir = getwd(), 
                     saveResults = TRUE, 
                     saveObject = TRUE, 
                     openShinyAPP = FALSE, ...) {
    
    ## arguments checking
    if (is.null(fcsFiles) || is.na(fcsFiles) || is.nan(fcsFiles)){
        stop("Wrong input fcsFiles!")
    }else if (length(fcsFiles) == 1 && file.info(fcsFiles)$isdir) {
        fcsFiles <- list.files(path = fcsFiles, pattern = ".fcs$", 
            full.names = TRUE)
        rawFCSdir <- dirname(fcsFiles)
    }else{
        if(dirname(fcsFiles[1]) == "."){
            rawFCSdir <- getwd()
        }else{
            rawFCSdir <- dirname(fcsFiles[1])  
        }
    }
    setwd(rawFCSdir)
    if(length(fcsFiles) < 1)
        stop("No FCS file found, please select your fcsFiles!")
    if(!all(file.exists(fcsFiles)))
        stop("Can not find file(s):", fcsFiles[which(!file.exists(fcsFiles))])
    
    if (length(markers) == 1 && file.exists(markers)) {
        markers <- as.character(read.table(markers, sep = "\t", 
                                           header = TRUE)[, 1])
    }
    if (is.null(markers) || length(markers) < 1) 
        stop("no marker selected!")
  
    mergeMethod <- match.arg(mergeMethod)
    
    if (!is.null(fixedNum) && !(is.numeric(fixedNum))) 
        stop("clusterSampleSize must be a numeric number!")
    
    transformMethod <- match.arg(transformMethod)
    dimReductionMethod <- match.arg(dimReductionMethod) 
    
    if(missing(clusterMethods)){
        clusterMethods <- "Rphenograph"
    }else{
        clusterMethods <- match.arg(clusterMethods, several.ok = TRUE)
    }
    
    if(missing(visualizationMethods)){
        visualizationMethods <- "tsne"
    }else{
        visualizationMethods <- match.arg(visualizationMethods, several.ok = TRUE)
    }
    
    progressionMethod <- match.arg(progressionMethod)
    
    if (!(is.numeric(clusterSampleSize))) 
        stop("clusterSampleSize must be a numeric number!")
    
    
    ## print arguments for user info
    message("Input arguments:")
    message("* Project Name: ")
    message(projectName, "\n")
    message("* Input FCS files for analysis:\n ")
    message(paste0("  -", basename(fcsFiles), "\n"))
    message("* Markers:\n ")
    message(paste0("  -", markers, "\n"))
    message("* Data merging method: ")
    message(mergeMethod, "\n")
    message("* Data transformation method: ")
    message(transformMethod, "\n")
    message("* Dimensionality reduction method: ")
    message(dimReductionMethod, "\n")
    message("* Data clustering method(s): ")
    message(clusterMethods, "\n")
    message("* Data visualization method(s): ")
    message(visualizationMethods, "\n")
    message("* Subset progression analysis method: ")
    message(progressionMethod, "\n\n")
    
    #set.seed(seed)
    ## get transformed, combined exprs data
    message("Extract expression data...")
    exprs_data <- cytof_exprsMerge(fcsFiles, comp = ifCompensation, verbose = FALSE, 
                                   transformMethod = transformMethod, 
                                   markers = markers, 
                                   mergeMethod = mergeMethod, fixedNum = fixedNum, ...)
    message("  ", nrow(exprs_data), " x ", ncol(exprs_data), " data was extracted!\n")
    
    
    ## dimension reduced data, a list
    message("Dimension reduction...")
    alldimReductionMethods <- unique(c(visualizationMethods, dimReductionMethod))
    allDimReducedList <- lapply(alldimReductionMethods, 
                                cytof_dimReduction
                                , data = exprs_data
                                , markers = markers
                                , umap_neighbor = umap_neighbor
                                , umap_min_dist = umap_min_dist
                                , tsneSeed = seed)
    names(allDimReducedList) <- alldimReductionMethods
    
    
    ## cluster results, a list
    message("Run clustering...")
    #set.seed(seed)
    cluster_res <- lapply(clusterMethods, cytof_cluster, 
                          ydata = allDimReducedList[[dimReductionMethod]], 
                          xdata = exprs_data[, markers],
                          Rphenograph_k = Rphenograph_k,
                          FlowSOM_k = FlowSOM_k,
                          flowSeed = seed)
    names(cluster_res) <- clusterMethods
    
    
    ## progression analysis results, a list  
    ## NOTE, currently only the first cluster method results 
    ## are used for progression visualization(by default: cluster_res[[1]])
    message("Progression analysis...")   
    progression_res <- cytof_progression(data = exprs_data, 
                                         cluster = cluster_res[[1]], 
                                         method = progressionMethod,
                                         out_dim = 4,
                                         clusterSampleSize = clusterSampleSize)
    
    ## overall list of markers
    ## using cytof_exprsMerge to generate list instead of
    ## read.FCS, for consistent marker naming
    message("Listing markers used for dimension reduction...")
    markerlist <- markers
    
    ## original fcs sample names
    message("Stashing sample names...")
    names <- sub("^.+/", "", unique(sub(".fcs$", "", fcsFiles)))
    samples <- as.list(NULL)
    for(i in seq_along(names)){
      samples[[i]] <- names[i]
    }
    
    ## wrap the results
    message("Wrapping results...")
    analysis_results <- list(expressionData = exprs_data,
                             dimReductionMethod = dimReductionMethod,
                             visualizationMethods = alldimReductionMethods, #visualizationMethods,
                             dimReducedRes = allDimReducedList,
                             clusterRes = cluster_res, 
                             progressionRes = progression_res,
                             projectName = projectName,
                             rawFCSdir = rawFCSdir,
                             resultDir = resultDir,
                             dimRedMarkers = markers,
                             sampleNames = samples)
     
    
    ## save the results
    message("Analysis DONE, saving the results...") 
    cytof_writeResults(analysis_results = analysis_results,
                       saveToRData = saveObject,
                       saveToFCS = saveResults,
                       saveToFiles = saveResults)
    
    if(openShinyAPP){
        cytofkitShinyAPP(RData = analysis_results)
    }
    
    invisible(analysis_results)
}


#' check the package update news
#' 
#' @return Opens .Rd file of package update news
#' 
#' @export
cytofkitNews <- function() 
{
    newsfile <- file.path(system.file(package = "cytofkit2"),
                          "NEWS.Rd")
    file.show(newsfile)
}


