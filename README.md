# FLOWMAP
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
