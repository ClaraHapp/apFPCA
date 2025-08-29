### Application to Seismology Data (Northridge Earthquake Simulations) ###

### load required packages / utility functions
library(tidyfun) # for preprocessing
library(tidyverse) # for preprocessing
library(furrr) # for preprocessing
library(mgcv) # for preprocessing
library(funData) # for funData representation
library(MFPCA) # for MFPCA calculation
library(ggplot2) # for plotting
source("utils.R") # for warping utilities

### Preprocessing

if(file.exists("../data/seissol_tweedie_smooth.rds"))
{
  seissol <- readRDS("../data/seissol_tweedie_smooth.rds")
}  else {
  source("app_smooth.R")
} 

# select only those functions that are close to the epicenter (< 40 km)
ind <- which(seissol$hypo.dist < 40000)

# raw data
x <- seq(0, 30, by=.5)
raw <- funData(argvals = x, X = as.matrix(log1p(seissol$Bodenbewegung[ind, ])))


### Calculate SRVF warping & MFPCA

if(file.exists("../data/SRVF_seis_tw_smooth_rawPCA.Rdata"))
{
  load("../data/SRVF_seis_tw_smooth_rawPCA.Rdata")
} else {
  source("app_warp_MFPCA.R")
}
   
   
### Final analysis and plotting
   
source("app_plot.R")