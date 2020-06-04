# FLOWMAPR

This repository houses the code for the FLOW-MAP algorithm, which was developed in R and originally published in Zunder et al., "A Continuous Molecular Roadmap to iPSC Reprogramming Through Progression Analysis of Single Cell Mass Cytometry." Cell Stem Cell 2015.

This code has been reformatted and compiled into a single R package known as **FLOWMAPR**, for extensible analysis of single-cell datasets including FCS files generated from mass/flow cytometry or data.frames of other single-cell data processed in R. This README provides full instructions for how to install, set-up, and use FLOWMAPR. Other details, including a comparison with other visualization tools for single-cell data, will be published imminently and linked when available.

Version 2 released June 2020: 
1. Better scalability to large cell numbers (when using UMAP layout)
2. Removes bugs specified in issues from v1
3. Removes need for dependencies that were no longer maintained & causing install issues

# Navigation
<!-- [Code Status]() -->
[Getting Started: Installing FLOWMAPR](#install)  
[Updating FLOWMAPR](#update)  
[Running FLOWMAPR: Starting from FCS Files](#FLOWMAPR-FCS)  
[Running FLOWMAPR: Starting from a Dataframe in R](#FLOWMAPR-DF)  
[Example Code for FLOWMAP()](#example-code)  
[Example Data](#example-data)  
[Example Output](#example-output)  
[Example Code for FLOWMAPfromDF()](#example-code2)  
[Practical Guidelines for Running FLOWMAPR](#FLOWMAPR-guide)  
[Practical Guidelines for Post-Processing in Gephi](#gephi-guide)  
[Using the GUI](#gui-guide)  
[Authors and License](#info)  


## Getting Started
FLOWMAPR visualizations can be run on any 32- and 64-bit computer with at least 2.2-GHz processor running Windows or Mac OS X with ≥4 GB of RAM (16 GB preferred). FLOWMAPR was developed and tested on a 64-bit computer with a 2.7-GHz processor running Mac OS X (Version 10.11.6) with 16 GB of RAM.

The instructions below demonstrate how to install this package directly from Github to get the latest release.

### Software and Package Prerequisites:
Install version 3.3.0 or later of R. Users can install R by downloading the appropriate R-x.y.z.tar.gz  file from http://www.r-project.org and following the system-specific instructions. FLOWMAPR was developed and tested on version 3.3.0 of R. As of this release, we recommend using version 3.3.0.

FLOWMAPR depends on the following R libraries: igraph (version 1.2.1), Rclusterpp (version 0.2.3), SDMTools (version 1.1.221), robustbase (version 0.92.8), shiny (version 1.0.5), tcltk (version 3.4.3), rhandsontable (version 0.3.6), spade (version 1.10.4) and flowCore (version 1.44.2) from Bioconductor, and scaffold (version 0.1) published on the Nolan Lab GitHub. The versions provided are the R package versions for which this FLOWMAPR code has been tested.

In order to install a package from github, you will need the devtools package. You can install this package with the following commands:

```
install.packages("devtools")
library(devtools)
```

FLOWMAPR package depends on several packages, which can be installed using the below commands:

```
install.packages("igraph")
install.packages("robustbase")
install_github("cran/SDMTools") 
install_github("ParkerICI/vite")
source("http://bioconductor.org/biocLite.R")
biocLite("flowCore")
```

The GUI has package dependencies for multiple Shiny-related packages and Rhandsontable. Install these packages with:

```
install.packages("shiny")
install.packages("shinythemes")
install.packages("shinyFiles")
install.packages("shinydashboard")
install.packages("shinyalert")

library(devtools)
install_version("rhandsontable", version = "0.3.4", repos = "http://cran.us.r-project.org")
```

Lastly, FLOWMAPR utilizes the R/C++ implementation of ForceAtlas2 as made available in the scaffold package from the Nolan Lab. Instructions to install that package are available here: https://github.com/nolanlab/scaffold. After installing all dependencies, you can install the package with the following commands:

```
install_github("nolanlab/scaffold")

```

<a name="install"></a>
### Installing FLOWMAPR:

To currently get the FLOWMAPR R package up and working on your computer, once you have installed all package dependencies (see above):

1. Make a github account if you haven't already.
2. Get access to repo zunderlab/FLOWMAP for your github account.
3. Create a token by going to this site when logged into your github account: https://github.com/settings/tokens/new
4. Check off "repo" in the settings for your token.
5. Click generate token and copy/save the provided code (your PAT) somewhere.
6. Open R studio and load devtools using `library(devtools)`. If you don't have devtools you may have to install it with `install.packages("devtools")` and then use `library(devtools)`.
7. Type the following into R studio: `install_github(repo = "zunderlab/FLOWMAP", auth_token = "PAT")` but replace PAT in quotations with your code in quotations. This should start installing all library dependencies so it may take a bit to finish. Check that it finishes without ERROR messages, though it may print WARNINGS.

Typical installation time should take no more than 5 minutes for the most up-to-date FLOWMAPR package. However, total installation time will vary depending on the installation time of other required packages and the speed of your internet connectoin.

<a name="update"></a>
### Updating FLOWMAPR:

To quickly update your FLOWMAPR R package and get the latest version from GitHub:

1. Open R studio and load devtools using `library(devtools)`.
2. Type the following into R studio: `install_github(repo = "zunderlab/FLOWMAP", auth_token = "PAT")` but replace PAT in quotations with your code in quotations.
3. Load FLOWMAPR using `library(FLOWMAPR)`.

If the above commands run without error, you should have the latest version of FLOWMAPR.

## Running FLOWMAPR
<a name="FLOWMAPR-FCS"></a>

### Starting from FCS Files:

To run a FLOW-MAP analysis on your data set if you are using FCS files or an example data set:

0. **Make your data available and parseable by FLOWMAPR.**
* For MultiFLOW-MAP, you must specify the "files" variable as a directory wherein each subfolder represented samples at the same time. If FCS files in the same subdirectory that come from different time points (e.g. "ConditionA-d01.fcs" with "ConditionB-d02.fcs"), FLOWMAPR will pick one time label arbitrarily.
* Please make sure the time labels can be parsed and sorted with proper labels (e.g. "01", "02", "04", "06", "10" vs. "1", "2", "4", "6", "10" where it would sort as "1" and "10" first instead of "10" last).
* To properly label each condition within the timepoint, please put the Condition as the first part of the file name separated by "-" or "." characters (e.g. "ConditionA-d01.fcs" where "ConditionA" will be the condition label).
* Do not use any digits (i.e. 0-9) in the name of the FCS file unless they specify time. Change any labels for the conditions in the FCS file name to be alphabetical characters. Ex: Condition1-t24.fcs should be renamed to ConditionOne-t24.fcs or else the time label will be parsed as "124" instead of "24" for this file.
* Please note that when your FCS files are loaded into FLOWMAPR, any "Time" variables already in the data will be removed and overwritten with the "time" of each FCS file.
1. Once you have successfully installed and loaded FLOWMAPR using `library(FLOWMAPR)`, if you are working in R Studio, you should see `FLOWMAPR::FLOWMAP()` autocomplete if you type it into the command line.
2. Establish variable names (you can copy the way they are assigned from the FLOWMAP_run.R file to declare each variable).  Some variables you have to assign are:
  * `mode` - what type of FLOW-MAP you want to run, this can be "single" - one condition, multiple timepoints, "multi" - multiple conditions, multiple timepoints or "one" - one condition, one timepoint
  * `files` - the directory where you can find the FCS files to be used
  * `name.sort` - sort FCS files according to name in alphabetical/numerical order, default is set to TRUE for sorting
  * `seed.X` - an integer that sets the seed, can be re-used to reproduce results, default is set to 1
  * `var.remove` - any channels you want completely excluded from analysis, you can auto-generate a suggested vector of variables for removal using the `FLOWMAPR::SuggestVarRemove()` function and supplying `var.annotate` as well as an optional `var.to.remove` vector that describes a character string for channels that should be removed (for example, blank channels in FCS files for mass cytometry may retain the "Di" or "Dd" substring in the desc, like in "Ba138Di"), default is to remove these channels with "Di" or "Dd" in the desc
  * `var.annotate` - rename channels as you see fit, the names you provide will the ones used to print out the PDFs, you can auto-generate a suggested `var.annotate` based on the `desc` attribute in your FCS files by using the `FLOWMAPR::ConstructVarAnnotate()` function and supplying a single FCS file name (full path)
  * `clustering.var` - which channels to use to influence the graph shape, you can generate a set of suggested `clustering.var` using the `FLOWMAPR::SuggestClusteringVar()` function, supplying the variables `fcs.file.names` (complete set of FCS files to be used in analysis), `mode` (as in FLOWMAPR mode), `var.annotate`, `var.remove`, `top.num` (which specifies how many clustering variables you want to use, less than the total number of variables in the dataset)
  * `cluster.numbers` - how many clusters to generate from each subsampled file, recommended ratio 1:2 from subsample (if subsample = 1000, recommended cluster.numbers = 500), default is set to 100
  * `distance.metric` - choose "manhattan" or "euclidean" for most cases, default is set to "manhattan"
  * `subsamples` - how many cells to randomly subsample from each FCS file, default is set to 200
  * `downsample` - use density-dependent downsampling, in which case you may want to specify and pass optional variables `exclude.pctile`, `target.pctile`, `target.number`, `target.percent`, default is set to FALSE for downsampling
  * `k` - number of nearest neighbors to use for density calculation to determine number of edges per vertex; lower values provide higher resolution and enable better separation of rare populations, but are more sensitive to noise whereas higher values are less sensitive to noise but will only capture the more robust populations
  * `minimum` - minimum number of edges allotted based on density, affects connectivity, recommended default is 2
  * `maximum` - maximum number of edges allotted based on density, affects connectivity, recommended default is 5
  * `save.folder` - where you want the output files to be saved to, default is set to current directory or getwd() result
  * `savePDFs` - produce PDF files or only produce graphml files, default is set to TRUE for printing all results
  * `graph.out` - "ForceDirected" or "UMAP" layout for vizualization of FLOWMAP graph
  * `umap.n.neighbors` - number of neighbors to use for UMAP layout, must be <= `k`
  * `which.palette` - optional argument for savePDFs functionality, can be “jet” (rainbow) or “bluered” or “CB” (the colorblind-friendly option), default is set to "bluered"


3. Run `FLOWMAPR::FLOWMAP()` as a command in R Studio, but pass the variables that you assigned into FLOWMAP() function. A full example is provided below.
4. Check that it saves an output folder with reasonable looking graph and layout files.

<a name="FLOWMAPR-DF"></a>
### Starting from a Dataframe in R:

To run a FLOW-MAP analysis and generate FLOW-MAP graphs from data that you need to load/preprocess in R (differently than how FCS files are handled): 

0. **Complete all preprocessing steps on your data and make sure it is parseable by FLOWMAPR.** Here are the expected formats of your data for each FLOW-MAP run mode.
* **mode "one"** = one data.frame object;
* **mode "single"** = a list of data.frame objects where each element contains cells from a different timepoint;
* **mode "multi"** = a list of lists of data.frame objects where the first level of each list corresponds to different timepoints and sublists correspond to different conditions within that timepoint
* If your data starts as one single large data.frame, you can use the `FLOWMAPR::RestructureDF()` function to reformat your data. Pass `time.col.label` (default is "Time") and `condition.col.label` (default is NULL, only required for MultiFLOW-MAP), which it will use to separate the data into a list or list of lists of data.frames basaed on time and/or condition.
* **You must preprocess your data, including all subsampling, renaming or removing of channels, and transformation necessary for your data type.** `FLOWMAPfromDF()` has the option to cluster your data, but it does not have options to subsample/downsample, change channels, or transform data like `FLOWMAP()` does.
* Please note that you can use data that has already been clustered in this FLOW-MAP analysis (each row is a cluster with the median expression values). However, Counts/percent.total will not be accurately generated from the data. You can include a channel that reflects the size of each cluster (e.g. "cluster.counts"), but it will not be used to adjust node size in the graph in the autogenerated PDFs.

1. Once you have successfully installed and loaded FLOWMAPR using `library(FLOWMAPR)`, if you are working in R Studio, you should see `FLOWMAPR::FLOWMAPfromDF()` autocomplete if you type it into the command line.
2. Establish variable names (you can copy the way they are assigned from the FLOWMAP_run.R file to declare each variable).  Some variables you have to assign are:
  
  * `mode` - what type of FLOW-MAP you want to run, this can be "single" - one condition, multiple timepoints, "multi" - multiple conditions, multiple timepoints or "one" - one condition, one timepoint
  * `project.name` - a text label that will be appended to some of the files generated as results from the FLOW-MAP run
  * `df` - your data as a data.frame format object, a list of data.frame objects, or a list of lists of data.frame objects in R, it is expected that first level of each list corresponds to different timepoints and sublists correspond to different conditions (if applicable)
  * `time.col.label` - required variable, function will use the column with this label as the time label for each cell, default is set to "Time"
  * `condition.col.label` - variable that is only required for MultiFLOW-MAP runs to distinguish data from different conditions/treatments/timecourses, function will use the column with this label as the condition label for each cell, default is set to NULL
  * `clustering.var` - which channels to use to influence the graph shape
  * `distance.metric` - choose "manhattan" or "euclidean" for most cases, default is set to "manhattan"
  * `k` - number of nearest neighbors to use for density calculation to determine number of edges per vertex; lower values provide higher resolution and enable better separation of rare populations, but are more sensitive to noise whereas higher values are less sensitive to noise but will only capture the more robust populations
  * `minimum` - minimum number of edges allotted based on density, affects connectivity, recommended default is 2
  * `maximum` - maximum number of edges allotted based on density, affects connectivity, recommended default is 5
  * `save.folder` - where you want the output files to be saved to, default is set to current directory or getwd() result
  * `name.sort` - sort timepoints according to time label in alphabetical/numerical order, default is set to TRUE for sorting
  * `clustering` - cluster within each timepoint, in which case you will want to specify optional variable `cluster.numbers`, default is set to FALSE for clustering
  * `seed.X` - an integer that sets the seed, can be re-used to reproduce results, default is set to 1
  * `savePDFs` - produce PDF files or only produce graphml files, default is set to TRUE for printing all results
  * `graph.out` - "ForceDirected" or "UMAP" layout for vizualization of FLOWMAP graph
  * `which.palette` - optional argument for savePDFs functionality, can be “jet” (rainbow) or “bluered” or “CB” (the colorblind-friendly option), default is set to "bluered"

3. Run `FLOWMAPR::FLOWMAPfromDF()` as a command in R Studio, but pass the variables that you assigned into FLOWMAPfromDF() function. A full example is provided below.
4. Check that it saves an output folder with reasonable-looking graph and layout files. **From previous update: reproducible run.R files are NOT generated for FLOW-MAPs from matrix/dataframe.**

<a name="example-code"></a>
### Example Code for FLOWMAP():

```
library(FLOWMAPR)
seed.X <- 1
files <- "/Users/mesako/Desktop/SingleFLOWMAP"
name.sort <- FALSE
mode <- "single"
save.folder <- "/Users/mesako/Desktop"
var.annotate <- list("marker1" = "marker1", "marker2" = "marker2")
var.remove <- c()
clustering.var <- c("marker1", "marker2")
downsample <- FALSE
subsamples <- 200
cluster.numbers <- 100
distance.metric <- "manhattan"
k <- 5
minimum <- 2
maximum <- 5
graph.out <- 'UMAP'
savePDFs <- TRUE
which.palette <- "bluered"

FLOWMAPR::FLOWMAP(mode = mode, files = files, var.remove = var.remove, var.annotate = var.annotate,
                  clustering.var = clustering.var, cluster.numbers = cluster.numbers,
                  distance.metric = distance.metric, minimum = minimum, maximum = maximum,
                  save.folder = save.folder, subsamples = subsamples, k = k, graph.out = graph.out,
                  name.sort = name.sort, downsample = downsample, seed.X = seed.X,
                  savePDFs = savePDFs, which.palette = which.palette)
```

### Example Code for FLOWMAP() with SPADE downsampling:

```
library(FLOWMAPR)
files <- "/Users/mesako/Desktop/SingleFLOWMAP"
mode <- "single"
save.folder <- "/Users/mesako/Desktop"
var.annotate <- list("marker1" = "marker1", "marker2" = "marker2")
var.remove <- c()
clustering.var <- c("marker1", "marker2")
cluster.numbers <- 100
distance.metric <- "manhattan"
k<- 5
minimum <- 2
maximum <- 5
seed.X <- 1
subsamples <- FALSE
exclude.pctile <- 0.01
target.pctile <- 0.99
target.number <- NULL
target.percent <- NULL
name.sort <- FALSE
downsample <- TRUE
graph.out <- 'UMAP'
savePDFs <- TRUE
which.palette <- "bluered"

FLOWMAPR::FLOWMAP(mode = mode, files = files, var.remove = var.remove, var.annotate = var.annotate,
                  clustering.var = clustering.var, cluster.numbers = cluster.numbers,
                  distance.metric = distance.metric, minimum = minimum, maximum = maximum,
                  save.folder = save.folder, subsamples = subsamples,
                  name.sort = name.sort, downsample = downsample, seed.X = seed.X,
                  savePDFs = savePDFs, which.palette = which.palette,
                  exclude.pctile = exclude.pctile, target.pctile = target.pctile,
                  target.number = target.number, target.percent = target.percent)
```

#### Template for a FLOWMAPR Run using FLOWMAP() function in R

You can also access example FLOWMAPR runs as .R files, which you can then modify to work with your own data.

To access a SingleFLOW-MAP run without downsampling (random subsampling):
```
file.edit(system.file("tools/run_FLOWMAPR.R", package = "FLOWMAPR"))
```

To access a SingleFLOW-MAP run with SPADE downsampling:
```
file.edit(system.file("tools/run_downsample_FLOWMAPR.R", package = "FLOWMAPR"))
```
You will need to change the `files` and the `save.folder` to point to folders on your own computer. Furthermore, you will need to update the `var.annotate`, `var.remove`, and `clustering.var` parameters to reflect channels in your own data. You can change all other parameters based on your analysis needs.

<a name="example-data"></a>
### Example Data:

An example (synthetic) data set is available as raw FCS files with the FLOWMAPR package for testing purposes. You can access these data sets by finding their directory on your computer using the following commands after you have installed and loaded FLOWMAPR.

To access the SingleFLOW-MAP data (one condition, one time course data in a single folder of FCS files):
```
files <- system.file("extdata/SingleFLOWMAP", package = "FLOWMAPR")
```

To access the MultiFLOW-MAP data (two conditions, one time course data in a folder that contains subfolders of FCS files, wherein each subfolder is numbered according to sequential timepoints):
```
files <- system.file("extdata/MultiFLOWMAP", package = "FLOWMAPR")
```

Supply this `files` variable as the `files` parameter in a `FLOWMAPR::FLOWMAP()` command.

<a name="example-output"></a>
### Example Output:
Running the provided example files for SingleFLOW-MAP and MultiFLOW-MAP with the provided code should generate a folder of results in the `save.folder` specified by the user. The results folder will have the naming convention "YYYY-MM-DD_HH.MM.SS_ChosenFLOW-MAPMode_run" with the values filled in accordingly.

Upon successful completion of one of the main `FLOWMAP` functions, FLOWMAPR will generate outputs specified by `graph.out` parameter:
- "ForceDirected" selection will produce three graphml files: one is the initial graph produced by the algorithm without force-directed layout applied, one produced after the iterative force-directed layout with x-y coordinates for each node, and a final graphml with x-y coordinates and a relabeling of timepoints/conditions to the names scraped from the FCS files (if relevant). If the user chooses to save PDF outputs from FLOWMAPR, a folder with the same name as the final graphml generated with the added suffix "\_pdf" will be created that contains PDFs of the final graph colored by each parameter.
- "UMAP" selection will generate one graphml file produced by the algorithm without force-directed layout applied, one PDF containing the layout colored by all parameters in a grid, and one file containing the layout colored by time.
- In both cases, FLOWMAPR creates a text file named "YYYY-MM-DD_HH.MM.SS_FLOW-MAPR_run_settings_summary.txt" that will summarize all provided settings as well as a file named "run_FLOWMAPR.R" that can be used to reproduce the same FLOWMAPR results when opened and executed in R.

<a name="example-code2"></a>
### Example Code for FLOWMAPfromDF():
```
library(FLOWMAPR)
library(readxl)
mode <- "single"
save.folder <- "/Users/mesako/Desktop"
project.name <- "Example_FLOWMAP_Run"
clustering.var <- c("marker1", "marker2")
distance.metric <- "manhattan"
k <- 5
minimum <- 2
maximum <- 5
seed.X <- 1
name.sort <- FALSE
clustering <- FALSE
savePDFs <- TRUE
which.palette <- "bluered"

time.col.label <- "Time"
condition.col.label <- NULL

file <- "/Users/mesako/Downloads/Example-dataset.xlsx"
df <- read_excel(file)
df <- as.data.frame(df)
df.keep <- subset.data.frame(df, select = c("Time"))
df.transform <- subset.data.frame(df, select = setdiff(colnames(df), c("Time")))
df.transform <- apply(df.transform, 2, log) 
final.df <- cbind(df.keep, df.transform)
df <- FLOWMAPR::RestructureDF(final.df, time.col.label = time.col.label, 
                              condition.col.label = condition.col.label)$new.df

FLOWMAPR::FLOWMAPfromDF(mode = mode, df = df, project.name = project.name,
                        time.col.label = time.col.label, condition.col.label = condition.col.label,
                        clustering.var = clustering.var, distance.metric = distance.metric,
                        minimum = minimum, maximum = maximum, save.folder = save.folder,
                        name.sort = name.sort, clustering = clustering,
                        seed.X = seed.X, savePDFs = savePDFs, which.palette = which.palette)

```

<a name="FLOWMAPR-guide"></a>
### Practical Guidelines for Running FLOWMAPR:

FLOWMAPR is an R package for visualization of high-dimensional data, though ultimately producing visualizations is a subjective process. If you are not sure where to start or what settings to use, here are some basic guidelines for your analysis:

0. **Consider what your question is for your dataset.** FLOWMAPR can just be used to explore a given timecourse dataset, but generally your choice of settings will be simplified if you have a directed question. For example, you may be interested in how a biological process changes with different treatment conditions, or how different subpopulations change in a given set of markers during a process. 
1. Analyze your timecourse data by some conventional means (heatmaps, histogram, dotplots, contour plots) to get an intuition for the data and to be able to answer, at least partially, the following questions:
+ **What are some of the different subpopulations in my data, especially those of interest to my question?** What fraction of the total number of cells do each of these subpopulations represent? This info can inform whether you choose to downsample (`downsample <- TRUE`) or randomly subsample, or your clustering ratio (`subsample` : `cluster.numbers`).
+ **What markers do or do not change across the timecourse?** Generally, you will not want to include completely uninformative markers as clustering variables (`clustering.var`). This choice can be guided by the data (what you observe changing) and biological expert knowledge (what you know should change in the process you profiled). We have implemented a function `FLOWMAPR::SuggestClusteringVar()` that can generate a set of clustering variables based on variation within and between timepoints. We also suggest using PCA/tSNE and other visualizations to find variables that either contribute the most to each principal component or that clearly drive clusters/spread in the data.
+ **What are some expected trajectories?** That is, what are some cell subpopulations you expect to observe, and what changes should you observe in these cells over time. Any expert knowledge can help you know what to look for (confirmation in your results) after a FLOWMAPR run before you investigate novel findings.
+ **Are there any parameters that are useless or unrelated to your biological question?** You can try not removing any markers, but we generally recommend you include markers with no relevance (e.g. DNA, Event_length, Eubeads) in `var.remove` so that they do not carry through the analysis and no PDFs are generated for them. This step will save small amounts of time and make the graph less inconvenient to navigate in Gephi.
2. **Start off with a small number of clusters and generally keep the clustering ratio larger (`subsample` close to, either equal to or slightly less than, `cluster.numbers`).** Though the resulting figures may not be fully representative of the variation in the data, these settings will allow you to quickly iterate through different configurations of edge settings and different choices of clustering variables. Given that a graph with about 1200 total nodes or clusters takes 2 minutes to run, try starting with cluster numbers set to approx. 1200 / # of files. For example, if you have a single timecourse with 5 FCS files, try setting `cluster.numbers <- 250` and `subsamples <- 500`.
3. **You will need to do several iterations of FLOWMAPR to try to arrive at appropriate settings for `minimum`, `maximum`, and `clustering.var`.** The order in which you proceed depends on the results you see, so you may need to reverse the order of the steps (4-5) below. For example, the default edge settings may be wrong for any `clustering.var` you try, in which case you should tweak edge settings and then compare different `clustering.var`.
4. **Try using the default edge settings for `minimum` and `maximum` with different options for `clustering.var`.** These results will show you how informative different sets of markers are. You should try to narrow down to a particular marker set that you can use to refine the edge settings.
5. **Once you select `clustering.var`, you can change edge settings `minimum` and `maximum` to try to arrive at the maximal separation within your data.** Generally you want to arrive at a graph that best resolves difference and allows for spread of different trajectories in the data, so that it captures as much info from the high-dimensional shape of the data as possible. Here are some guides for how to tweak these edge settings:
+ If the graph is **too interconnected (hairball-like)**, try reducing `maximum`. You can try reducing `minimum` to 1, but generally we recommend that `minimum` is at least 2. Try moving `maximum` to being at most `minimum` + 1.
+ If the graph is **not interconnected enough (spiky, single nodes radiating out)**, try increasing the `minimum` and/or `maximum`.
+ Be careful as it's possible for graphs to essentially become tangled as they are processed with a force-directed layout. If results do not look useful, check for these tangles that can be resolved in Gephi. Additionally, the force-directed layout step is a computationally intensive and time-consuming step, so it is possible within the R package that the process does not complete. These graphs can be resolved to a stable shape in Gephi.
6. **Once you arrive at a FLOW-MAP with the "best" settings, repeat the analysis with multiple settings of `seed.X` to produce "technical replicates" of your analysis.** For steps that involve randomness, such as random subsampling of your FCS files, you will want to try to reproduce your FLOW-MAP figures with different samplings.

From anecdotal evidence, most datasets work well with setting `minimum` to 2 and `maximum` to anything from 5 to 20. Some datasets will show a "saturation point" where more edges allotted (a higher `maximum`) does not significantly change the graph shape so you can start with `maximum` set to 20. If you need more cohesiveness, increase `minimum`. If you need less cohesiveness, reduce `maximum`.

#### Timing for FLOWMAPR Runs

In a SingleFLOW-MAP with no downsampling (uses random subsampling), 1200 total nodes takes about 2 min to produce results (including PDFs). In comparison, 3000 total nodes takes about 6 min to produce all results, 6000 total nodes takes about 21 min, and 12000 total nodes takes about 59 min. These all ran with a `subsample` : `cluster.numbers` ratio of 2:1.

For users familiar with R who wish to only iterate through edge settings for the FLOW-MAP algorithm with the same starting dataset, we recommend a shortcut instead of starting with the `FLOWMAPR::FLOWMAP()`. Preprocessing your FCS files in R and saving the intermediate (as an .rds or .rdata object) to load into FLOWMAPR (instead using the `FLOWMAPR::FLOWMAPfromDF()` function) can save valuable time.

<a name="gephi-guide"></a>
### Practical Guidelines for Post-Processing in Gephi:

Producing aesthetically pleasing graphs is easier in Gephi. FLOWMAPR autogenerates PDF results so that the user can quickly scan through the resulting graphs and iterate through different settings. However, Gephi allows for greater customization of visual settings.

0. **Download and install Gephi.** Gephi is a free, open-source program available online at [http://www.gephi.org](http://www.gephi.org). We recommend using Gephi version "0.9.2-SNAPSHOT" as Gephi 0.9.1 has a bug that prevents users from changing node size of FLOW-MAP graphs using the "percent.total" parameter.

1. **Open your FLOW-MAP graphml file in Gephi.** We recommend you use the resulting graphml file that contains the substring "xy_orig_time" in the file name. A window will pop up called the "Import report," just hit the OK button.

2. **Set the node size in your graph.** The nodes will need to be the final size you intend for the graph before you run the force-directed layout.
* Go to the appearance panel usually in the top left. Make sure you have "Nodes" selected and click the button that shows three concentric circles inside one another (hovering should show the label "Size").
* Click on "Ranking" and use "percent.total" if you want the nodes to reflect the relative sizes of the clusters. If you did not cluster, you can pick one single size using "Unique" instead of "Ranking."
* Pick any size or window of sizes you prefer, as this may depend on how many nodes you have in your graph.

3. **Resolve your FLOW-MAP graph further using Gephi's ForceAtlas2 algorithm.** To do this, go to the Layout panel usually in the bottom left after opening your graphml file. Click the "Choose a Layout" field and select "ForceAtlas 2." We generally recommend not changing any of the default setting. Press Run to start resolving the graph in a force-directed layout.

* To speed up this process (and make it less computationally intensive), you can make edges invisible. Do so by pressing the button below the graph viewing window that looks like a line segment. If you hover over the button, it is called "Show Edges." Click to toggle it on and off.

* You may see regions of the graph that appear tangled. As long as ForceAtlas 2 is actively running (the Run button changes into a Stop button), you can click on nodes in the graph and while holding down the mouse button, move and manipulate them. The graph should respond as you move nodes.

* ForceAtlas 2 does not have a stopping time. You can continue to resolve the graph as long as you like, though it is recommended that you do so until there are no tangles and the graph stops changing.

* **Once you are happy with your graph layout**, we recommend that you toggle on the "Dissuade Hubs" option under Layout. This should be pressed while ForceAtlas 2 is still running. Wait until the graph finishes moving and then toggle on the "Prevent Overlap" option. This will spread the nodes so you can clearly see them as opposed to having them stack on top of each other. Once nodes approximately stop moving (they will continue to jiggle a little), hit the Stop button.

* **We recommend you save and export this new graph with the finalized layout** so you can consistently generate figures from this one graph. Gephi has its own file format (.gephi) and can also export .graphml files. We recommend saving both. If you export a graphml file, be sure to check the Options button and make sure that "Position (x,y)" is checked off.

4. **Use consistent visual settings to color the graph by marker expression.** You can set the color scheme by clicking on the paint palette in the Appearance pane. Make sure that "Nodes" are highlighted.

* "Partition" can be used to color by categorical variables (e.g. different timepoints or conditions) while "Ranking" can be used to color by marker expression or cluster size. For "Partition," you can set each separate value's color by clicking on each colored square and dragging it on the color wheel. In contrast, "Ranking" variables are colored on a gradient set by a few values. If you hover over the color bar, you can see arrows appear at the points where the color is defined. Doubleclicking on any arrow allows you to reassign the color using a color wheel, RGB values, or a Hex value. 

* We recommend using the bluered palette used in the FLOWMAPR package to color by marker expression. Cells with lowest expression levels will be blue, intermediate will be grey, and high will be red. You can set this up by clicking on the left, middle, and right arrows and assigning them to Hex values 1500FB, C3C3C7, D10100 respectively. **Press the Apply button to have the graph update with these new color settings.**

* However, to produce more accessible figures, you may want to consider using a **colorblind-friendly palette**. The colorblind option used in FLOWMAPR uses a blue-yellow scale. The hex color code for each arrow from left to right is 0072B2, C3C3C7, E69F00.

5. **Save PDFs or image files of your graph colored by each marker.** Once you have your desired visual settings, click on the "Preview" button on the top bar. **Click the Refresh button on the bottom left to make your graph appear.** This graph will essentially be a duplicate of what you produced in the "Overview" mode with a few exceptions.

* Preview mode shows your graph with edges, in particular curved edges. We recommend saving PDFs with edges invisible as the appearance of many grey edges will make your figure look messy. To turn off edges, click on "Show Edges" option under "Edges" on the left and click "Refresh" to update the graph. If you choose to display edges, you can adjust your edge settings (thickness, color, curve, opacity, etc.) in that section.

* We recommend that you also print your PDF with no borders around the nodes. We find that node borders generally makes the graph look messy and the colors harder to distinguish. To turn this off, click the value next to "Border Width" (generally defaults to 1.0) under the "Nodes" header. You can change this value to 0 and click "Refresh" to update the graph.

* Once you are happy with the graph appearance (having clicked "Refresh"), you can then click on the Export "SVG/PDF/PNG" to generate any of these file types.

<a name="gui-guide"></a>
## Using the GUI
0. Make sure all FCS files to be used are in one folder.
* If you plan to analyze files using FLOW-MAP mode "multi" (multiple conditions, multiple timepoints), you will need to set up a CSV file that outlines the FCS file paths (see section "CSV Format for Multiple Analysis") and select this folder instead of the FCS file folder.
* Regardless of the mode, make sure that your FCS files meet the naming conventions described at the beginning of the section ["Running FLOWMAPR: Starting from FCS Files"](#FLOWMAPR-FCS).
1. Run the command `FLOWMAPR::LaunchGUI()` and a dialogue box with the header "FLOWMAPR" should appear. 
2. Enter all of the relevant information which pertains to the type of experiment that is being analyzed, including if the data should be analyzed by FLOWMAPR mode "one", "single", or "multi".
3. When you press "Submit", a new window should appear that runs in the Shiny interface.

**The usage from this point on differs depending on what mode (e.g. "multi" or "single" or "one") is used.**

**For mode "one" (one condition, one timepoint):**
1. Select the FCS file to be analyzed in "Uploaded Order".
2. Press "Generate Parameters".
3. An interactive table will appear with all the parameters as well as options for selecting and deselecting them as clustering or removed variables. Removed variables will not be in the final generated graph, such as in the graphml file or the image PDFs.
4. You must check at least one or more of the parameters for clustering. These variables are used both for clustering (calculation of similarity) and for calculating edge distances during the graph building steps. If you want to rename a parameter, you can click on the name under the "annotate" column and type a new name. 
5. Press "Run FLOWMAPR" once the appropriate parameters have been checked and renamed to run the FLOW-MAP algorithm and generate all requested FLOWMAPR results (PDFs, graphml files, etc. in a new folder).

**For mode "single" (one condition, multiple timepoints):**
1. Enter in the order of the FCS files that you wish to use. Generally, files will be used in an alphanumerical order by time, but here you can specify the ordering if the naming system does not reflect the order you want. 
2. Press "Generate Parameters".
3. Two things will happen: an interactive table will appear with all the parameters and options for selecting how parameters should be used for analysis, and the menus for "Similar Fields" and "Different Fields" will autopopulate as an aid to help you process channels between the files.
4. If any channel needs to be merged, select the files from the "Different Fields" window, enter the new merged name in "Select New Merge Name", and press "Merge Selected Diff". This will automatically remove the channels from "Different Fields", add the merged name to "Similar Fields", and will update the table with new annotations.
5. The different parameters will by default be checked for removal. You must check at least one or more of the parameters for clustering. If you want to rename a parameter, click on the name under "annotate" and type a new name. 
6. Press "Run FLOWMAPR" once the appropriate parameters have been checked and renamed to run the FLOW-MAP algorithm and generate all requested FLOWMAPR results (PDFs, graphml files, etc. in a new folder).

**For mode "multi" (multiple conditions, multiple timepoints):**
1. Select the CSV file that has the corresponding FCS file paths. How the CSV file should be arranged (i.e. what information is put in the columns/rows) will be shown in the following section.
2. Press "Input CSV" once the CSV is selected in the box.
3. If any channel needs to be merged, select the files from the "Different Fields" window, enter the new merged name in "Select New Merge Name", and press "Merge Selected Diff". This will automatically remove the channels from "Different Fields", add the merged name to "Similar Fields", and will update the table with new annotations.
4. The different parameters will by default be checked for removal. You must check at least one or more of the parameters for clustering. If you want to rename a parameter, click on the name under "annotate" and type a new name. 
5. Press "Run FLOWMAPR" once the appropriate parameters have been checked and renamed to run the FLOW-MAP algorithm and generate all requested FLOWMAPR results (PDFs, graphml files, etc. in a new folder).

**Mode "one-special" (multiple conditions, one timepoint) is not currently implemented in the GUI. Please used the package in R.**

## CSV Format for Multiple Analysis
The CSV file should be formatted as shown below in cells. Assume the timepoints are to the left of the cells, and increase downwards starting with the first timepoint. The labels under time will be used to label the timepoints in the final FLOW-MAP graph. 

| Time |            ConditionA             |             ConditionB            |
|:----:|:---------------------------------:|:---------------------------------:|
|  01  |  ConditionA Time01 FCS File Path  |  ConditionB Time01 FCS File Path  |
|  02  |  ConditionA Time02 FCS File Path  |                                   |
|  03  |                                   |  ConditionB Time03 FCS File Path  |
|  04  |  ConditionA Time04 FCS File Path  |  ConditionB Time04 FCS File Path  |

## Troubleshooting
Here are some common issues and suggestions for how to fix them:

**1. The program crashes during ForceAtlas2.**

**Solution:** Usually these crashes originate during the ForceAtlas2 algorithm stage, which is programmed in C++ called from R. The source of this bug is still unclear so we recommend trying to circumvent the error by changing the seed of the FLOWMAPR analysis.

**2. In the `FLOWMAPfromDF()` function - the software does not recognize input.**

**Solution:** This error will appear if the provided input (the dataframe in R) does not match the mode specified by the user. We suggest that you doublecheck that the mode of analysis is what you intended and also check that the input is one of the accepted inputs for that mode. 

**3. In the `FLOWMAP()` function -- the software does not recognize input.**

**Solution:** This error will appear if the provided input (the full path of the folder or of the FCS files) does not match the mode specified by the user. We suggest that you doublecheck that the mode of analysis is what you intended and also check that the input is one of the accepted inputs for that mode. 

**4. The graph from the graphml file or the PDFs have unexpected labels (especially for time or condition).**

**Solution:** If you are performing a FLOWMAPR run using the `FLOWMAP()` function, check that your FCS files (and folders, if applicable) are named according to the acceptable naming convention. Condition and time names are scraped from these file paths. If you are performing a FLOWMAPR run using the `FLOWMAPfromDF()` function, check that you correctly specify the names of the Condition and Time columns in the dataframe, and that the labels contained in those columns are correct.

<a name="info"></a>
## Authors

* **Eli Zunder** - *Initial work*
* **Melissa Ko** - *R Package Development*
* **Corey Williams - R Package Development
* **Sarah Goggin - GUI, Version 2 
* **Rohit Rustagi - Initial GUI work

See also the list of [contributors](https://github.com/zunderlab/FLOWMAP/graphs/contributors) who participated in this project.

## License

This project is licensed under the GNU GPLv3 License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments

* Our package uses codes from the SPADE R package and scaffold R package. We thank authors Benedict Anchang and Federico Gherardini for their guidance in using their packages.
