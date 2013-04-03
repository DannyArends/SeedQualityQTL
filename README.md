## Seed Quality QTL mapping pipeline with MQM in R
This repository contains the code used to produce the results for the paper:

- Visualizing the genetic landscape of Arabidopsis seed performance: 
Ronny V. L. Joosen*, Danny Arends*, Leo Willems, Wilco Ligterink, Henk Hilhorst, Ritsert C. Jansen - Plant Physiology (2011)

This repository contains a single script.R which can be sourced into R 

    setwd("/path/to/script")
    source("script.R")

The analysis is then performed using a configuration file, and is started by the command:

    DoAnalysisOn("config.txt","/path/to/config")

This will run through several steps:

 - 1.0 Read the configFILE
 - 1.1 Check if the required parameters are supplied 
 - 2.0 Remove Outliers
 - 2.1 Creating a normalized cross object
 - 3.0 Main analysis loop
 - 3.1 PLot of basic statistics, Genetic Map, Distributions, Etc
 - 4.0 QTL by MQM
 - 4.1 QTL by scanone and scan.two
 - 4.2 QTL plots
 - 5.0 Based on MQM results -> Circleplots, Heatmaps, Effect plots, Grouping
 - 5.1 PCA analysis on the first 3 PCs
 - 5.2 Output of SIF network images

Additionally a log.txt will be created.

