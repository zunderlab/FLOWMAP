# FLOWMAP
====================

To execute the FLOW-MAP analysis, download a clone of the current commit.

Open the FLOWMAP_run.R in RStudio and tweak all variables as needed within that file including:
- prefolder
- files
- save.folder
- var.annotate
- var.remove
- per
- minimum
- maximum
- distance.metric
- subsamples
- cluster.numbers
- seed.X
- clustering.var

Then run the entire FLOWMAP_run.R file. You will not need to run any of the other R scripts in the directory as run.R should call them as needed.
====================


This repository is to clean up (and hopefully improve) the FLOWMAP algorithm code, which was developed in R and published in Zunder et al. A Continuous Molecular Roadmap to iPSC Reprogramming Through Progression Analysis of Single Cell Mass Cytometry. Cell Stem Cell. 2015.

Tasks:

1. update for compatibility with new R packages (igraph0 to igraph, etc.)
2. make code modular
3. make code standalone so it can run from start to finish in R
	a. remove dependence on SPADE
	b. remove dependence on GEPHI - or modify GEPHI to make it FLOWMAP specialized???
4. improve downsampling/clustering
5. simplify graph building??
6. build GUI with R/Shiny (and GEPHI/Java????)
