# Masterfile for Aggregating Distributional Treatment Effects 

# Toggle for MCMC
rerun_MCMC <- FALSE # only turn this onto true if you have a very powerful computer or a lot of time at your disposal 
if(rerun_MCMC == TRUE){
  source()
  
  
}



# whether or not you run the MCMC yourself or use the supplied rstan files, you now need to make graphics and tables
# counterintuitively we must start with the graphics because I processed the model output into quantiles within the graphics files 
# if you want things to be easy, here is the order you must run them in 

## FIGURES

# produces figures 1 and 2
source("R code/graphics and computation of quantiles from log normal models business variables.R", print.eval  = TRUE)
source("R code/graphics for no pooling models and site specific output lognormal.R", print.eval  = TRUE) 
source("R code/graphics and computation of quantiles from log normal consumption variables.R", print.eval  = TRUE)
# produces figure 3
source("R code/posterior-predictive-graphic-comparison-lognormal-pareto-full-sim.R", print.eval  = TRUE)  
# produces figure 4
source("R code/graphics and computation of quantiles from PLN model.R", print.eval  = TRUE)
source("R code/graphics and computation of quantiles composite tail model.R", print.eval  = TRUE) 
# produces figures 5 and 6 
source("R code/graphics and computation of quantiles from tailored hierarchical pdf lognormal models consumption types pb split.R", print.eval  = TRUE)
source("R code/graphics and computation of quantiles from tailored hierarchical pdf lognormal models pb split.R", print.eval  = TRUE)



## TABLES
# produces tables 1 and 2 
source("R code/tables-full-results-profit-consumption.R", print.eval = TRUE) # rm token is here. the above is fully checked. 
