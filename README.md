# FLOWMAP

This repository is to clean up (and hopefully improve) the FLOWMAP algorithm code, which was developed in R and published in Zunder et al. A Continuous Molecular Roadmap to iPSC Reprogramming Through Progression Analysis of Single Cell Mass Cytometry. Cell Stem Cell. 2015.

## Getting Started

To get the FLOWMAPR R package up and working on your computer:

1. Make a github account if you haven't already.
2. Get access to repo zunderlab/FLOWMAP for your github account.
3. Create a token by going to this site when logged into your github account: https://github.com/settings/tokens/new
4. Check off "repo" in the settings for your token.
5. Click generate token and copy/save the provided code (your PAT) somewhere.
6. Open R studio and load devtools using `library(devtools)`. If you don't have devtools you may have to install it with `install.packages("devtools")` and then use `library(devtools)`.
7. Type the following into R studio: `install_github(repo = "zunderlab/FLOWMAP", auth_token = “PAT”)` but replace PAT in quotations with your code in quotations. This should start installing all library dependencies so it may take a bit to finish. Check that it finishes without ERROR messages, though it may print WARNINGS.

## Running FLOW-MAP

To run a FLOWMAP analysis on your data set or an example data set (INSTRUCTIONS INCOMING):

0. Make your data available and parseable by FLOWMAP. For MultiFLOWMAP, you must specify the "files" variable as a directory wherein each subfolder represented samples at the same time. Please label times sequentially from 1 ... n, even if that does not reflect the actual experimental timepoints. To properly label each condition within the timepoint, please put the Condition as the first part of the file name separated by "-" or "." characters.
1. Once you have successfully loaded FLOWMAPR, if you are working in R Studio, you should see `FLOWMAPR::FLOWMAP()` autocomplete if you type it into the command line.
2. Establish variable names (you can copy the way they are assigned from the FLOWMAP_run.R file to declare each variable).  Some variables you have to assign are:
  * files - the directory where you can find the FCS files to be used
  * save.folder - where you want the output files to be saved to
  * var.annotate - rename channels as you see fit, the names you provide will the ones used to print out the PDFs
  * var.remove - any channels you want completely excluded from analysis
  * per - affects connectivity, recommended default is 1
  * minimum - minimum number of edges allotted based on density, affects connectivity, recommended default is 2
  * maximum - maximum number of edges allotted based on density, affects connectivity, recommended default is 3
  * distance.metric - choose manhattan or euclidean
  * subsamples - how many cells to randomly subsample from each FCS file, RECOMMENDATION PENDING
  * cluster.numbers - how many clusters to generate from each subsampled file, recommended ratio 1:2 from subsample (if subsample = 1000, recommended cluster.numbers = 500)
  * seed.X - set this for reproducibility
  * clustering.var - which channels to use to influence the graph shape

3. Run `FLOWMAPR::FLOWMAP()` as a command in R Studio, but pass the variables that you assigned into FLOWMAP function. A full example is provided below.
4. Check that it saves an output folder with reasonable looking PDFs and graphml files.

### Example Code:

```
files <- "/Users/mesako/Desktop/FLOWMAP-Synthetic-Data-20161216/SingleFLOWMAP"
save.folder <- "/Users/mesako/Desktop"
file.format <- "*.fcs"
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
FLOWMAPR::FLOWMAP(files = files, file.format = file.format, var.remove = var.remove,
                  var.annotate = var.annotate, clustering.var = clustering.var,
                  cluster.numbers = cluster.numbers, subsamples = subsamples,
                  distance.metric = distance.metric, minimum = minimum, maximum = maximum,
                  per = per, save.folder = save.folder, shuffle = TRUE, name.sort = FALSE)
```

## Contributing

????

## Versioning

???

## Authors

* **Eli Zunder** - *Initial work*
* **Melissa Ko** - *???*

See also the list of [contributors](https://github.com/zunderlab/FLOWMAP/graphs/contributors) who participated in this project.

## License

This project is licensed under the ??? License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* ???
