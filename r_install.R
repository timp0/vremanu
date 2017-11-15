##Install needed R pkgs


install.packages(c("tidyverse", "googlesheets", "ggjoy"), repo="http://cran.wustl.edu/")

source("https://bioconductor.org/biocLite.R") 

biocLite(c("ontoCAT", "GenomicAlignments"), ask=FALSE)

