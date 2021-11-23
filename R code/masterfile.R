# Masterfile for Aggregating Distributional Treatment Effects 

# Toggle for MCMC
rerun_MCMC <- FALSE # only turn this onto true if you have a very powerful computer or a lot of time at your disposal 
if(rerun_MCMC == TRUE){
  source()
  
  
}



# whether or not you run the MCMC yourself or use the supplied rstan files, you now need to make graphics and tables
# counterintuitively we must start with the graphics because I processed the model output into quantiles within the graphics files 
# if you want things to be easy, here is the order you must run them in 

source("R code/graphics and computation of quantiles from log normal models business variables.R")
source("R code/graphics for no pooling models and site specific output lognormal.R") 
source("R code/graphics and computation of quantiles from log normal consumption variables.R") # rm token is here. the above is fully checked. 