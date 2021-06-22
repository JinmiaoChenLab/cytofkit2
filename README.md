cytofkit2: an integrated mass cytometry data analysis pipeline
============

**NOTE**: <u>This is the development version of cytofkit2 package</u>

### cytofkit2

This package is designed to facilitate the analysis workflow of mass cytometry data with automatic subset identification and mapping of cellular progression. Both command line and a GUI client are provided for executing the workflow easily.

### Installation

#### 1. Install R and Rstudio

If you have never used R, please install R and Rstudio following the steps below:

- Download the proper R version for your operation system from [R download page](http://cran.stat.nus.edu.sg).

- Double-click the downloaded R installation file, install with all the defaults.

- Download the proper Rstudio version for your operation system from [here](https://www.rstudio.com/products/rstudio/download/).

- Double-click the downloaded Rstudio installation file, install with all the defaults.


<u>**Special Notes for Mac Users**</u>

For Mac OS X 10.8 or later, you need to install XQuartz to support the GUI:

* Download the disk image (dmg) file for [XQuartz](http://xquartz.macosforge.org).

* Open the dmg file by double-clicking on it, you'll find XQuartz.pkg, double-click on it to run the installer, clicking through all the defaults.

* After the installer runs, you'll have to **restart your mac computer**.

#### 2. Install python

#### 3. Install cytofkit2 package

The offical and stable version, please refer to 

- [Bioconductor](https://www.bioconductor.org/packages/cytofkit2/)
- [GitHub](https://github.com/JinmiaoChenLab/cytofkit2)

Install the stable version from Bioconductor, use:

``` r
## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("cytofkit2")
```

Install this development version, use:

``` r
install.packages("reticulate")

#### check if python is installed
vv <- system("pyv=\"$(python -V)\" | echo $pyv| grep \"Python\"")
if(vv){
 print("Python is installed")
}

#### install python package umap-learn
library(reticulate)
py_install("umap-learn")

#### install cytofkit2 from github
if(!require(devtools)){
  install.packages("devtools") # If not already installed
}
devtools::install_github("JinmiaoChenLab/cytofkit2")
```

### Usage

After successfully installing the cytofkit2 package, run the following codes to open the cytofkit GUI:

``` r
library("cytofkit2")
cytofkit_shiny_dashboard()
```

<u>Check the following vignettes for more details:</u>

- [cytofkit: Analysis Pipeline](https://www.bioconductor.org/packages/release/bioc/vignettes/cytofkit/inst/doc/cytofkit_workflow.html)    
- [cytofkit: Quick Start](https://www.bioconductor.org/packages/release/bioc/vignettes/cytofkit/inst/doc/cytofkit_example.html)   
- [cytofkit: ShinyAPP Tutorial](https://www.bioconductor.org/packages/release/bioc/vignettes/cytofkit/inst/doc/cytofkit_shinyAPP.html)    




