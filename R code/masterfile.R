# Masterfile for Aggregating Distributional Treatment Effects 


##### MAIN PAPER #####

## MCMC Models 

# Toggle below for rerunning all the MCMC scripts for the main paper
# This is possible, but not necessary: I have provided all MCMC output you need to make the tables and figures in the repository 
# Only switch this to TRUE if you have Rstan working AND a very powerful computer OR a lot of time at your disposal. Caveat Emptor! 
rerun_MCMC <- FALSE
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
source("R code/tables-full-results-profit-consumption.R", print.eval = TRUE) 
# produces table 3
source("R code/bayesian tests of equality.R", print.eval = TRUE) 

# Rm token is here. it is time for an omnibus test of the above. 

##### APPENDICES #####