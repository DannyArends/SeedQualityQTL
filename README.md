## Seed Quality QTL mapping pipeline with R/qtl & MQM

This repository contains the updated code used to produce the results for the paper:

- Visualizing the genetic landscape of Arabidopsis seed performance: 
Ronny V. L. Joosen*, Danny Arends*, Leo Willems, Wilco Ligterink, Henk Hilhorst, Ritsert C. Jansen - Plant Physiology (2011)

This repository contains a single script.R which can be sourced into R 

    setwd("/path/to/script")
    source("script.R")

The analysis is then performed using a configuration file, and is started by the command:

    results <- DoAnalysisOn("config.txt","/path/to/config")

This will run through several steps (depending on what is enabled in the config.txt):

 - 1.0 Read the configuration file
 - 1.1 Check if the required parameters are supplied 
 - 2.0 Remove outliers from data
 - 2.1 Normalization of phenotypes and creating a new cross object
 - 3.0 Main analysis loop
 - 3.1 Plots of basic statistics, genetic map, phenotype distributions, etc
 - 4.0 QTL mapping by MQM
 - 4.1 QTL mapping by scanone and scan.two
 - 4.2 QTL plots
 - 5.0 Based on MQM results -> Circleplots, Heatmaps, Effect plots, Grouping
 - 5.1 PCA analysis on the first 3 PCs
 - 5.2 Output of SIF network images

During the analysis a file log.txt will also be created.

### Want to contribute? Great!

1. Clone it.
2. Compile it.
3. Run it.
4. Modify some code. (Search -> 'TODO')
5. Go back to 2, or
6. Submit a patch

You can also just post comments on code / commits.
Danny Arends

### Disclaimer

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License,
version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License, version 3, for more details.

A copy of the GNU General Public License, version 3, is available
at [http://www.r-project.org/Licenses/GPL-3](http://www.r-project.org/Licenses/GPL-3 "GPL-3 Licence")

Copyright (c) 2009-2014 [Danny Arends](http://www.dannyarends.nl), with modificationd from: Ronny Joosen

