# FLOWMAP

This repository houses the FLOWMAP algorithm code, which was developed in R and originally published in Zunder et al. A Continuous Molecular Roadmap to iPSC Reprogramming Through Progression Analysis of Single Cell Mass Cytometry. Cell Stem Cell. 2015.

## Code Status
### Known Issues (“Bugs”)
* fix percent/count after SPADE downsampling in final graph

### Planned Features
* improve “load” to allow loading cluster tables/matrix 
* ~~improve specifying mode “single” or “multi” type of FLOW-MAP~~
* ~~fix “detecting type” method to instead test inputs and throw error if not correct type for provided data~~
* ~~introduce SPADE downsampling variables to specify by user~~
* ~~update summary print method to include new variables~~

### Requested Features
* edges ranked within Gephi to remove/add to scale connectivity as we go
* add density rank of the vertices to nodes so that amount of downsampling can be changed, to remove/add nodes
* density rank can also be used to change amount of outlier removal
* adjust ratio of edges within same time point to edges between time points
* save template (panel, experimental set-up), at least save FLOWMAPR_run.R with set variables

## Getting Started

The instructions below demonstrate how to install this package directly from Github to get the latest release.

### Software and Package Prerequisites:
Install version 3.3.0 or later of R. As of this release, we recommend using version 3.3.0.

In order to install a package from github, you will need the devtools package. You can install this package with the following commands:

```
install.packages("devtools")
library(devtools)
```

FLOWMAPR package depends on several packages, which can be installed using the below commands:

```
install.packages("Rclusterpp") 
install.packages("SDMTools") 
install.packages("igraph")
install.packages("robustbase")
source("http://bioconductor.org/biocLite.R")
biocLite("flowCore")
```

Lastly, FLOWMAPR utilizes the R/C++ implementation of ForceAtlas2 as made available in the scaffold package from the NOlan Lab. Instructions to install that package are available here: https://github.com/nolanlab/scaffold. After installing all dependencies, you can install the package with the following commands:

```
library(devtools)
install_github("nolanlab/scaffold")
```

### Installing FLOWMAPR:

To currently get the FLOWMAPR R package up and working on your computer:

1. Make a github account if you haven't already.
2. Get access to repo zunderlab/FLOWMAP for your github account.
3. Create a token by going to this site when logged into your github account: https://github.com/settings/tokens/new
4. Check off "repo" in the settings for your token.
5. Click generate token and copy/save the provided code (your PAT) somewhere.
6. Open R studio and load devtools using `library(devtools)`. If you don't have devtools you may have to install it with `install.packages("devtools")` and then use `library(devtools)`.
7. Type the following into R studio: `install_github(repo = "zunderlab/FLOWMAP", auth_token = "PAT")` but replace PAT in quotations with your code in quotations. This should start installing all library dependencies so it may take a bit to finish. Check that it finishes without ERROR messages, though it may print WARNINGS.

### Installing the GUI

1. The GUI has a package dependancy for Shiny, TclTk, and Rhandsontable. Install these with:
```
install.packages("shiny")
install.packages("tcltk")
install.packages("rhandsontable")
```
2. The GUI is accessed by running:
```
FLOWMAPR::launch_GUI()
```

## Running FLOW-MAP

To run a FLOWMAP analysis on your data set or an example data set:

0. Make your data available and parseable by FLOWMAP. For MultiFLOWMAP, you must specify the "files" variable as a directory wherein each subfolder represented samples at the same time. Please label times sequentially from 1 ... n, even if that does not reflect the actual experimental timepoints. To properly label each condition within the timepoint, please put the Condition as the first part of the file name separated by "-" or "." characters.
1. Once you have successfully installed and loaded FLOWMAPR using `library(FLOWMAPR)`, if you are working in R Studio, you should see `FLOWMAPR::FLOWMAP()` autocomplete if you type it into the command line.
2. Establish variable names (you can copy the way they are assigned from the FLOWMAP_run.R file to declare each variable).  Some variables you have to assign are:
  * `files` - the directory where you can find the FCS files to be used
  * `save.folder` - where you want the output files to be saved to
  * `mode` - what type of FLOW-MAP you want to run, this can be "single" - one condition, multiple timepoints, "multi" - multiple conditions, multiple timepoints or "one" - one condition, one timepoint
  * `var.annotate` - rename channels as you see fit, the names you provide will the ones used to print out the PDFs
  * `var.remove` - any channels you want completely excluded from analysis
  * `per` - affects connectivity, recommended default is 1
  * `minimum` - minimum number of edges allotted based on density, affects connectivity, recommended default is 2
  * `maximum` - maximum number of edges allotted based on density, affects connectivity, recommended default is 3
  * `distance.metric` - choose manhattan or euclidean
  * `subsamples` - how many cells to randomly subsample from each FCS file, RECOMMENDATION PENDING
  * `cluster.numbers` - how many clusters to generate from each subsampled file, recommended ratio 1:2 from subsample (if subsample = 1000, recommended cluster.numbers = 500)
  * `seed.X` - set this for reproducibility
  * `clustering.var` - which channels to use to influence the graph shape
  * `shuffle` - shuffle cells used in each sample
  * `name.sort` - sort FCS files according to name in alphabetical/numerical order
  * `downsample` - use SPADE density-dependent downsampling, in which case you may want to specify and pass optional variables `exclude.pctile`, `target.pctile`, `target.number`, `target.percent`


3. Run `FLOWMAPR::FLOWMAP()` as a command in R Studio, but pass the variables that you assigned into FLOWMAP function. A full example is provided below.
4. Check that it saves an output folder with reasonable looking PDFs and graphml files.

### Example Code:

```
files <- "/Users/mesako/Desktop/FLOWMAP-Synthetic-Data-20161216/SingleFLOWMAP"
mode <- "single"
save.folder <- "/Users/mesako/Desktop"
var.annotate <- list("marker1" = "marker1", "marker2" = "marker2")
var.remove <- c()
per <- 1
minimum <- 2
maximum <- 3
distance.metric <- "manhattan"
subsamples <- 200
cluster.numbers <- 100
seed.X <- 1
clustering.var <- c("marker1", "marker2")
set.seed(seed.X)
FLOWMAPR::FLOWMAP(files = files, var.remove = var.remove, var.annotate = var.annotate,
                  clustering.var = clustering.var, cluster.numbers = cluster.numbers,
                  subsamples = subsamples, distance.metric = distance.metric,
                  minimum = minimum, maximum = maximum, per = per,
                  save.folder = save.folder, mode = mode,
                  shuffle = TRUE, name.sort = FALSE, downsample = FALSE)
```

### Example Data:

An example (synthetic) data set is available as raw FCS files with the FLOWMAPR package for testing purposes. You can access these data sets by finding their directory on your computer using the following commands after you have installed and loaded FLOWMAPR.

To access the SingleFLOWMAP data (one condition, one time course data in a single folder of FCS files):
```
files <- system.file("extdata/SingleFLOWMAP", package = "FLOWMAPR")
```

To access the MultiFLOWMAP data (two conditions, one time course data in a folder that contains subfolders of FCS files, wherein each subfolder is numbered according to sequential timepoints):
```
files <- system.file("extdata/MultiFLOWMAP", package = "FLOWMAPR")
```

Supply this `files` variable as the `files` parameter in a `FLOWMAPR::FLOWMAP()` command.

## Using the GUI
0. Make sure all FSC files that are to be tested are within one folder.
1. When `FLOWMAPR::launch_GUI()` is run, a dialogue box of the header "FLOWMAP" should appear. 
2. Enter in all of the relevant information which pertains to the type of experiment that is being analyzed.
3. When "Submit" is pressed, a new window should appear which runs with Shiny.
4. Press "Generate Order List", which will display all the FCS files in the folder.
5. Select the order of the FCS files that you wish to use. If there is only 1 FCS file, choose the number of the file that is shown.
6. Press "Generate Parameters."
7. Two things will now happen: an interactive table will appear with all the parameters and options for selecting and deselecting them for analysis, and the menus for "Similar Fields" and "Different Fields" will autopopulate as an aid to help you choose relevant channels.
8. The different parameters will by default be checked for removal, and the user must check the parameters for clustering. If the user wishes to rename a parameter, they can click on the name under "annotate" and type a new name. 
9. Press "Write File" once the appropriate parameters have been checked and renamed to run the FLOW-MAP algorithm and generate the FLOW-MAP results (PDFs, graphml files, etc. in a new folder).

## Contributing

????

## Versioning

???

## Authors

* **Eli Zunder** - *Initial work*
* **Melissa Ko** - *???*

See also the list of [contributors](https://github.com/zunderlab/FLOWMAP/graphs/contributors) who participated in this project.

## License

This project is licensed under the GNU GPLv3 License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments

* ???
