
# 1. Unzip the ItemS1.zip file and set path2GMDF to be the path to the unzipped ItemS1.zip directory.

setwd(path2GMDF)

# 2. Make sure to install the dependencies
library(plyr)
library(rliger)
source("GMDF_wrapper.R")

# 3. To reproduce the CD8 T cell pan-cancer programs:
rslts.pancancer<-GMDF_combine_pancancer()

# 4. Run the toy example below

input<-readRDS("ToyExample.rds")
rslts1<-GMDF_wrapper(E = input$E, a = input$a, k = 5, k1 = 2,N1 = 1)
dir.create("GMDF.toy.example/")
rslts<-GMDF_wrapper(E = input$E, a = input$a, k = 5, k1 = 2,N1 = 5,outputdir = "GMDF.toy.example/")

# 5. To use GMDF for other datasets:
# rslts<-GMDF_wrapper(E, a, k, k1,N1,outputdir)

