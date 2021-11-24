# posterior predictive comparison of lognormal and pareto models 
# Rachael Meager
# First Version: June 2018 

### Notes ###

# Analysis needs to be done before you can run this, obviously
# I sincerely apologise for the atrocious indenting here. I would like to especially apologise to Jonathan Huggins who taught me better.

### Preliminaries and Data Intake ###

# clear the workspace to avoid gremlins and past globals from past irresponsible scripts
# but we can't do this if the masterfile is being used to run the script, so we check that first:
if(exists("masterfile_run") == "FALSE"){
  rm(list = ls())
}


installation_needed  <- FALSE
loading_needed <- TRUE
package_list <- c('ggplot2', 'rstan','reshape','reshape2','coda','xtable', 'dplyr', 'Runuran', 'testthat',
                  "MCMCpack", "gtools", 'gPdtest', 'fBasics',"PtProcess", "VGAM", "stargazer")
if(installation_needed){for (p in package_list) {install.packages(p, repos='http://cran.us.r-project.org')}}
if(loading_needed){lapply(package_list, require, character.only = TRUE)}


fig_scale = 0.4

# Load data
load("output/microcredit_profit_lognormal_tailored_hierarchical_pdf_output_5000_iters.RData")
stan_fit_table_profit <- stan_fit_table
stan_fit_profit <- stan_fit
codafit_stan_draws_profit <- as.matrix(stan2coda(stan_fit_profit))
data_split_profit <- data_split
data_profit <- data
load("output/microcredit_expenditures_lognormal_tailored_hierarchical_pdf_output_5000_iters.RData")
stan_fit_table_expenditures <- stan_fit_table
stan_fit_expenditures <- stan_fit
codafit_stan_draws_expenditures <- as.matrix(stan2coda(stan_fit_expenditures))
data_split_expenditures <- data_split
data_expenditures <- data
load("output/microcredit_revenues_lognormal_tailored_hierarchical_pdf_output_5000_iters.RData")
stan_fit_table_revenues <- stan_fit_table
stan_fit_revenues <- stan_fit
codafit_stan_draws_revenues <- as.matrix(stan2coda(stan_fit_revenues))
data_split_revenues <- data_split
data_revenues <- data


### FUNCTIONS I WILL NEED


convert_logit_2_prob <- function(vector_betas){
  K <- length(vector_betas)
  vector_betas <- vector_betas - vector_betas[K]
  probs <- rep(NA,K)
  denominator_excluding_K <- (1+sum(exp((vector_betas[1:(K-1)]))))
  probs[K] <- 1/denominator_excluding_K
  for(k in 1:(K-1)){
    probs[k] <- exp(vector_betas[k])/denominator_excluding_K
  }
  
  return(probs)
}


get_posterior_intervals_function <- function(vector_draws){
  posterior_interval <- quantile(vector_draws, c(0.025, 0.25, 0.5, 0.75,0.975))
  return(posterior_interval)
  print(posterior_interval)
}


## DATA FOR THE 1 TAIL OUTCOMES


## expenditures

codafit_stan_draws <- codafit_stan_draws_expenditures
S <- dim(codafit_stan_draws)[1]
M <- 2 # mixture components
component_probs <- matrix(NA,nrow=S,ncol=M)
quantile_vec <- seq(0.05, 0.95,0.1)
N <- length(quantile_vec)
quantile_value_output_control <- matrix(NA, nrow=S, ncol=N)


for(s in 1:S){
  control_logmean <-rnorm(1,codafit_stan_draws[s,"mu[1]"],codafit_stan_draws[s,"sd_mu[1]"]) 
  control_logsd <- exp(rnorm(1, codafit_stan_draws[s,"sigma_control[1]"],codafit_stan_draws[s,"sd_sigma_control[1]"]) )
  control_logit_params <- c(rnorm(1,codafit_stan_draws[s,"beta_full[1,1]"], codafit_stan_draws[s, "sigma[1,1]"]), rnorm(1,codafit_stan_draws[s,"beta_full[2,1]"], codafit_stan_draws[s, "sigma[2,1]"]))
  control_probs <- convert_logit_2_prob(control_logit_params)
  
  # a question for me: does this structure REALLY propagate like the Gaussian structure does? 
  # One might ask "Should I have been doing everything on the log scale?" but this will not save you in a mixture model. 

  for(n in 1:N){

  u <- quantile_vec[n]
  if(u < control_probs[1]){ quantile_value_output_control[s,n] <-  0 }
  if(u > control_probs[1]){ quantile_value_output_control[s,n] <-  qlnorm(((u-control_probs[1])/control_probs[2]), meanlog = control_logmean, sdlog = control_logsd, lower.tail = TRUE, log.p = FALSE) }
  
  } # close forloop indexed by N 

}

expenditures_lognormal_quantile_value_output_control_median <- apply(quantile_value_output_control,2,median)
expenditures_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)


## revenues! 

codafit_stan_draws <- codafit_stan_draws_revenues
S <- dim(codafit_stan_draws)[1]
M <- 2 # mixture components
component_probs <- matrix(NA,nrow=S,ncol=M)
quantile_vec <- seq(0.05, 0.95,0.1)
N <- length(quantile_vec)
quantile_value_output_control <- matrix(NA, nrow=S, ncol=N)



for(s in 1:S){
  control_logmean <-rnorm(1,codafit_stan_draws[s,"mu[1]"],codafit_stan_draws[s,"sd_mu[1]"]) 
  control_logsd <- exp(rnorm(1, codafit_stan_draws[s,"sigma_control[1]"],codafit_stan_draws[s,"sd_sigma_control[1]"]) )
  control_logit_params <- c(rnorm(1,codafit_stan_draws[s,"beta_full[1,1]"], codafit_stan_draws[s, "sigma[1,1]"]), rnorm(1,codafit_stan_draws[s,"beta_full[2,1]"], codafit_stan_draws[s, "sigma[2,1]"]))
  control_probs <- convert_logit_2_prob(control_logit_params)
  
  # a question for me: does this structure REALLY propagate like the Gaussian structure does? 
  # One might ask "Should I have been doing everything on the log scale?" but this will not save you in a mixture model. 
  
  for(n in 1:N){
    
    u <- quantile_vec[n]
    if(u < control_probs[1]){ quantile_value_output_control[s,n] <-  0 }
    if(u > control_probs[1]){ quantile_value_output_control[s,n] <-  qlnorm(((u-control_probs[1])/control_probs[2]), meanlog = control_logmean, sdlog = control_logsd, lower.tail = TRUE, log.p = FALSE) }
    
  } # close forloop indexed by N 
}
revenues_lognormal_quantile_value_output_control_median <- apply(quantile_value_output_control,2,median)
revenues_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)



### FOR PROFIT


codafit_stan_draws <- codafit_stan_draws_profit
S <- dim(codafit_stan_draws)[1]
M <- 3 # mixture components
component_probs <- matrix(NA,nrow=S,ncol=M)
quantile_vec <- seq(0.05, 0.95,0.1)
N <- length(quantile_vec)
quantile_value_output_control <- matrix(NA, nrow=S, ncol=N)



for(s in 1:S){

  control_logmean_1 <- rnorm(1,codafit_stan_draws[s,"mu[1]"],codafit_stan_draws[s,"sd_mu[1]"])
  control_logsd_1 <- exp(rnorm(1,codafit_stan_draws[s,"sigma_control[1]"],codafit_stan_draws[s,"sd_sigma_control[1]"]))
  
  control_logmean_2 <- rnorm(1,codafit_stan_draws[s,"mu[2]"],codafit_stan_draws[s,"sd_mu[2]"])
  control_logsd_2 <- exp(rnorm(1,codafit_stan_draws[s,"sigma_control[2]"],codafit_stan_draws[s,"sd_sigma_control[2]"]))
  
  control_logit_params <- rnorm(3, codafit_stan_draws[s,c("beta_full[1,1]", "beta_full[2,1]", "beta_full[3,1]")], codafit_stan_draws[s,c("sigma[1,1]", "sigma[2,1]", "sigma[3,1]")]  )
  control_probs <- convert_logit_2_prob(control_logit_params)

  
  
  for(n in 1:N){

    u <- quantile_vec[n]
    if(u < control_probs[1]){ quantile_value_output_control[s,n]  <-  - qlnorm((u/control_probs[1]), meanlog = control_logmean_1, sdlog = control_logsd_1, lower.tail = TRUE, log.p = FALSE) }
    if(control_probs[1] < u & u < (control_probs[1]+control_probs[2])){ quantile_value_output_control[s,n]  <-  0 }
    if(u > (control_probs[1]+control_probs[2])){ quantile_value_output_control[s,n]  <-  qlnorm(((u-(control_probs[1]+control_probs[2]))/(1-control_probs[1]-control_probs[2])),  meanlog = control_logmean_2, sdlog = control_logsd_2, lower.tail = TRUE, log.p = FALSE)  }

 
  } # close forloop indexed by n 
  
} # closes forloop indexed by s 

profit_lognormal_quantile_value_output_control_median <- apply(quantile_value_output_control,2,median)
profit_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)




# now the paretos! 

# Load data
load("output/microcredit_profit_tailored_hierarchical_pdf_output.RData")
stan_fit_table_profit <- stan_fit_table
stan_fit_profit <- stan_fit
codafit_stan_draws_profit <- as.matrix(stan2coda(stan_fit_profit))
data_split_profit <- data_split
data_profit <- data
load("output/microcredit_expenditures_tailored_hierarchical_output.RData")
stan_fit_table_expenditures <- stan_fit_table
stan_fit_expenditures <- stan_fit
codafit_stan_draws_expenditures <- as.matrix(stan2coda(stan_fit_expenditures))
data_split_expenditures <- data_split
data_expenditures <- data
load("output/microcredit_revenues_tailored_hierarchical_output.RData")
stan_fit_table_revenues <- stan_fit_table
stan_fit_revenues <- stan_fit
codafit_stan_draws_revenues <- as.matrix(stan2coda(stan_fit_revenues))
data_split_revenues <- data_split
data_revenues <- data

### FUNCTIONS I WILL NEED


convert_logit_2_prob <- function(vector_betas){
  K <- length(vector_betas)
  vector_betas <- vector_betas - vector_betas[K]
  probs <- rep(NA,K)
  denominator_excluding_K <- (1+sum(exp((vector_betas[1:(K-1)]))))
  probs[K] <- 1/denominator_excluding_K
  for(k in 1:(K-1)){
    probs[k] <- exp(vector_betas[k])/denominator_excluding_K
  }
  
  return(probs)
}


pareto_quantile_function <- function( quantile, shape, min_location ){
  denominator <- ( 1-quantile)^(1/shape)
  q <- min_location/denominator
  return(q)
  print(q)
  
}

get_posterior_intervals_function <- function(vector_draws){
  posterior_interval <- quantile(vector_draws, c(0.025, 0.25, 0.5, 0.75,0.975))
  return(posterior_interval)
  print(posterior_interval)
}

## FOR THE 1 TAIL OUTCOMES

## revenues! 

codafit_stan_draws <- codafit_stan_draws_revenues
S <- dim(codafit_stan_draws)[1]
M <- 2 # mixture components
component_probs <- matrix(NA,nrow=S,ncol=M)
quantile_vec <- seq(0.05, 0.95,0.1)
N <- length(quantile_vec)
quantile_value_output_control <- matrix(NA, nrow=S, ncol=N)
data_location <- min(data_split_revenues$revenues)



for(s in 1:S){
  

  control_shape <- exp(rnorm(1, codafit_stan_draws[s,"control_shape"],codafit_stan_draws[s,"control_sigma"]) )
  control_logit_params <- c(rnorm(1,codafit_stan_draws[s,"beta_full[1,1]"], codafit_stan_draws[s, "sigma[1,1]"]), rnorm(1,codafit_stan_draws[s,"beta_full[2,1]"], codafit_stan_draws[s, "sigma[2,1]"]))
  control_probs <- convert_logit_2_prob(control_logit_params)
  
  for(n in 1:N){
    u <- quantile_vec[n]
    if(u < control_probs[1]){ quantile_value_output_control[s,n] <-  0 }
    if(u > control_probs[1]){ quantile_value_output_control[s,n] <-  pareto_quantile_function(((u-control_probs[1])/control_probs[2]), control_shape, data_location) }

  } # close forloop indexed by N 
  
}
revenues_pareto_quantile_value_output_control_median <- apply(quantile_value_output_control,2,median)
revenues_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)


# for expenditures
codafit_stan_draws <- codafit_stan_draws_expenditures
S <- dim(codafit_stan_draws)[1]
M <- 2 # mixture components
component_probs <- matrix(NA,nrow=S,ncol=M)
quantile_vec <- seq(0.05, 0.95,0.1)
N <- length(quantile_vec)
quantile_value_output_control <- matrix(NA, nrow=S, ncol=N)
data_location <- min(data_split_expenditures$expenditures)



  for(s in 1:S){
    
    
    control_shape <- exp(rnorm(1, codafit_stan_draws[s,"control_shape"],codafit_stan_draws[s,"control_sigma"]) )
    control_logit_params <- c(rnorm(1,codafit_stan_draws[s,"beta_full[1,1]"], codafit_stan_draws[s, "sigma[1,1]"]), rnorm(1,codafit_stan_draws[s,"beta_full[2,1]"], codafit_stan_draws[s, "sigma[2,1]"]))
    control_probs <- convert_logit_2_prob(control_logit_params)
    
    for(n in 1:N){
      u <- quantile_vec[n]
      if(u < control_probs[1]){ quantile_value_output_control[s,n] <-  0 }
      if(u > control_probs[1]){ quantile_value_output_control[s,n] <-  pareto_quantile_function(((u-control_probs[1])/control_probs[2]), control_shape, data_location) }
      
    } # close forloop indexed by N 
    
  }

expenditures_pareto_quantile_value_output_control_median <- apply(quantile_value_output_control,2,median)
expenditures_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)


### FOR PROFIT


codafit_stan_draws <- codafit_stan_draws_profit
S <- dim(codafit_stan_draws)[1]
M <- 3 # mixture components
component_probs <- matrix(NA,nrow=S,ncol=M)
quantile_vec <- seq(0.05, 0.95,0.1)
N <- length(quantile_vec)
quantile_value_output_control <- matrix(NA, nrow=S, ncol=N)
data_location_pos <- min(data_split_profit[[2]]$profit)
data_location_neg <- min(-data_split_profit[[1]]$profit)


for(s in 1:S){

  control_shape_1 <- exp(rnorm(1,codafit_stan_draws[s,"control_shape[1]"], codafit_stan_draws[s,"control_sigma[1]"]))
    control_shape_2 <- exp(rnorm(1,codafit_stan_draws[s,"control_shape[2]"], codafit_stan_draws[s,"control_sigma[2]"]))
  control_logit_params <- rnorm(3, codafit_stan_draws[s,c("beta_full[1,1]", "beta_full[2,1]", "beta_full[3,1]")], codafit_stan_draws[s, c( "sigma[1,1]", "sigma[2,1]", "sigma[3,1]") ])
  control_probs <- convert_logit_2_prob(control_logit_params)

  
  
  for(n in 1:N){
    u <- quantile_vec[n]
    if(u < control_probs[1]){ quantile_value_output_control[s,n]  <-  -pareto_quantile_function((u/control_probs[1]), control_shape_1, data_location_neg) }
    if(control_probs[1] < u & u < (control_probs[1]+control_probs[2])){ quantile_value_output_control[s,n]  <-  0 }
    if(u > (control_probs[1]+control_probs[2])){ quantile_value_output_control[s,n]  <-  pareto_quantile_function(((u-(control_probs[1]+control_probs[2]))/(1-control_probs[1]-control_probs[2])), control_shape_2, data_location_pos) }
    
 
  } # close forloop indexed by n 
  
} # closes forloop indexed by s 


profit_pareto_quantile_value_output_control_median <- apply(quantile_value_output_control,2,median)
profit_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)



### NOW MAKE THE TABLE AS A DATAFRAME

data <- data.frame(profit, revenues, expenditures, treatment)
data <- data[complete.cases(data),]

posterior_predicted_data <- data.frame(quantile(data$revenues[data$treatment==0], seq(0.05, 0.95, 0.1), na.rm = TRUE),
                                         revenues_lognormal_quantile_value_output_control_median,
                                       revenues_pareto_quantile_value_output_control_median,
                                      
                                       quantile(data$expenditures[data$treatment==0], seq(0.05, 0.95, 0.1), na.rm = TRUE),
                                       expenditures_lognormal_quantile_value_output_control_median,
                                       expenditures_pareto_quantile_value_output_control_median,
                                    
                                       quantile(data$profit[data$treatment==0], seq(0.05, 0.95, 0.1), na.rm = TRUE),
                                       profit_lognormal_quantile_value_output_control_median,
                                       profit_pareto_quantile_value_output_control_median
                           
                                       )
colnames(posterior_predicted_data) <- c("Actual Revenues Data", "Revenues predicted by Lognormal", "Revenues pedicted by Pareto",
                                        "Actual Expenditures Data", "Expenditures predicted by Lognormal", "Expenditures pedicted by Pareto",
                                        "Actual Profit Data", "Profit predicted by Lognormal", "Profit pedicted by Pareto")


posterior_predicted_data <- round(posterior_predicted_data , 2)


posterior_predicted_data <- t(posterior_predicted_data)

stargazer(posterior_predicted_data, summary = FALSE, digits = 0)



# and now, graphics! 
x_axis <- rep(seq(0.05, 0.95, 0.1),3)
revenues_quantiles <- c(quantile(data$revenues[data$treatment==0], seq(0.05, 0.95, 0.1), na.rm = TRUE),
               revenues_lognormal_quantile_value_output_control_median,
               revenues_pareto_quantile_value_output_control_median)
type <- c(rep("Data",10), rep("LogNormal", 10), rep("Pareto",10))
rev_dataframe <- data.frame(x_axis, revenues_quantiles, type)

p_rev <- ggplot(rev_dataframe, aes(x=x_axis, y=revenues_quantiles, group=type)) +
  geom_line(aes(color=type, fill=type), alpha = 1)+
  geom_point(aes(color=type)) +
   coord_cartesian(ylim=c(0, 1000)) + xlab("Quantiles") +ylab("") +
  ggtitle("Posterior Predictive Distributions of Revenues Data") +
  theme(plot.title = element_text(size = fig_scale*20)) + theme(legend.title = element_blank()) + theme(legend.text=element_text(size=fig_scale*14)) + 
  theme(axis.text = element_text(size=fig_scale*14)) +  theme(axis.title.y = element_text(size = fig_scale*14)) +  theme(axis.title.x = element_text(size = fig_scale*14))
p_rev
filename_as_string <- "output/posterior_predictive_revenues"
myfile <- paste0(filename_as_string,".pdf")
dev.copy(pdf,myfile, width = fig_scale*10, height = fig_scale*6)
dev.off()

expenditures_quantiles <- c(quantile(data$expenditures[data$treatment==0], seq(0.05, 0.95, 0.1), na.rm = TRUE),
                        expenditures_lognormal_quantile_value_output_control_median,
                        expenditures_pareto_quantile_value_output_control_median)
type <- c(rep("Data",10), rep("LogNormal", 10), rep("Pareto",10))
exp_dataframe <- data.frame(x_axis, expenditures_quantiles, type)

p_exp <- ggplot(rev_dataframe, aes(x=x_axis, y=expenditures_quantiles, group=type)) +
  geom_line(aes(color=type, fill=type), alpha = 1)+
  geom_point(aes(color=type)) +
  coord_cartesian(ylim=c(0, 1000)) + xlab("Quantiles") +ylab("") +
  ggtitle("Posterior Predictive Distributions of Expenditures Data") +
  theme(plot.title = element_text(size = fig_scale*20))  + theme(legend.title = element_blank()) + theme(legend.text=element_text(size=fig_scale*14)) + 
  theme(axis.text = element_text(size=fig_scale*14)) +  theme(axis.title.y = element_text(size = fig_scale*14)) +  theme(axis.title.x = element_text(size = fig_scale*14))
p_exp
filename_as_string <- "output/posterior_predictive_expenditures"
myfile <- paste0(filename_as_string,".pdf")
dev.copy(pdf,myfile, width = fig_scale*10, height = fig_scale*6)
dev.off()


profit_quantiles <- c(quantile(data$profit[data$treatment==0], seq(0.05, 0.95, 0.1), na.rm = TRUE),
                        profit_lognormal_quantile_value_output_control_median,
                        profit_pareto_quantile_value_output_control_median)
type <- c(rep("Data",10), rep("LogNormal", 10), rep("Pareto",10))
profit_dataframe <- data.frame(x_axis, profit_quantiles, type)

p_profit <- ggplot(profit_dataframe, aes(x=x_axis, y=profit_quantiles, group=type)) +
  geom_line(aes(color=type, fill=type), alpha = 1)+
  geom_point(aes(color=type)) +
  coord_cartesian(ylim=c(-200, 1000)) + xlab("Quantiles") +ylab("") +
  ggtitle("Posterior Predictive Distributions of Profit Data") +
  theme(plot.title = element_text(size = fig_scale*20))  + theme(legend.title = element_blank()) + theme(legend.text=element_text(size=fig_scale*14)) + 
  theme(axis.text = element_text(size=fig_scale*14)) +  theme(axis.title.y = element_text(size = fig_scale*14)) +  theme(axis.title.x = element_text(size = fig_scale*14))
p_profit
filename_as_string <- "output/posterior_predictive_profit"
myfile <- paste0(filename_as_string,".pdf")
dev.copy(pdf,myfile, width = fig_scale*10, height = fig_scale*6)
dev.off()

