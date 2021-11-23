# graphics for all parameters from all site results from no pooling and partial pooling models together
# Rachael Meager
# First Version: nov 2015, This version: August 2021 


### Notes ###

# Analysis needs to be done before you can run this, obviously

### Preliminaries and Data Intake ###

# clear the workspace to avoid gremlins and past globals from past irresponsible scripts
# but we can't do this if the masterfile is being used to run the script, so we check that first:
if(exists("masterfile_run") == "FALSE"){
  rm(list = ls())
}

installation_needed  <- FALSE
loading_needed <- TRUE
package_list <- c('ggplot2', 'rstan','reshape','reshape2','coda','xtable', 'dplyr', 'Runuran', 'testthat',
                  "MCMCpack", "geoR", "gtools", 'gPdtest', 'fBasics',"PtProcess", "VGAM", "MASS","quantreg",
                  "boot", "gridExtra")
if(installation_needed){install.packages(package_list, repos='http://cran.us.r-project.org')}
if(loading_needed){lapply(package_list, require, character.only = TRUE)}



# Load data
load("output/microcredit_profit_lognormal_tailored_no_pooling_pdf_output_5000_iters.RData")
stan_fit_table_profit_no_pooling <- stan_fit_table
stan_fit_profit_no_pooling <- stan_fit
codafit_stan_draws_profit_no_pooling <- as.matrix(stan2coda(stan_fit_profit_no_pooling))

load("output/microcredit_profit_lognormal_tailored_hierarchical_pdf_output_5000_iters.RData")
stan_fit_table_profit <- stan_fit_table
stan_fit_profit <- stan_fit
codafit_stan_draws_profit <- as.matrix(stan2coda(stan_fit_profit))
nodata_split_profit <- data_split
data_profit <- data


### FUNCTIONS I WILL NEED

convert_logit_2_prob <- function(vector_betas){
  K <- length(vector_betas)
  Z <- sum(exp((vector_betas[1:(K)])))
  vector_betas <- vector_betas - log(Z)
  probs <- rep(NA,K)
  denominator <- sum(exp((vector_betas[1:(K)])))
  probs[K] <- 1/denominator
  for(k in 1:(K)){
    probs[k] <- exp(vector_betas[k])/denominator
  }
  
  return(probs)
}


get_posterior_intervals_function <- function(vector_draws){
  posterior_interval <- quantile(vector_draws, c(0.025, 0.25, 0.5, 0.75,0.975))
  return(posterior_interval)
  print(posterior_interval)
}

# Profit no pooling tabulation 

codafit_stan_draws <- codafit_stan_draws_profit_no_pooling
S <- dim(codafit_stan_draws)[1]
M <- 3 # mixture components
K <- as.numeric(length(unique(data$site)))
component_probs <- matrix(NA,nrow=S,ncol=M)
quantile_vec <- seq(0.05, 0.95,0.1)
N <- length(quantile_vec)
quantile_value_output_control <- array( rep(NA, N*K*S), dim= c(K,S,N))
quantile_value_output_treatment <- array( rep(NA, N*K*S), dim= c(K,S,N))
quantile_value_output_difference <- array( rep(NA, N*K*S), dim= c(K,S,N))
mu_k_name_vec <- grep("mu_k", colnames(codafit_stan_draws_profit), value = TRUE)
tau_k_name_vec <- grep("tau_k", colnames(codafit_stan_draws_profit), value = TRUE)
sigma_control_k_name_vec <- grep("sigma_control_k", colnames(codafit_stan_draws_profit), value = TRUE)
sigma_TE_k_name_vec <- grep("sigma_TE_k", colnames(codafit_stan_draws_profit), value = TRUE)
beta_control_name_vec <- matrix(c("beta_k[1,1,1]", "beta_k[1,2,1]","beta_k[1,3,1]",
                                  "beta_k[2,1,1]", "beta_k[2,2,1]","beta_k[2,3,1]",
                                  "beta_k[3,1,1]", "beta_k[3,2,1]","beta_k[3,3,1]",
                                  "beta_k[4,1,1]", "beta_k[4,2,1]","beta_k[4,3,1]",
                                  "beta_k[5,1,1]", "beta_k[5,2,1]", "beta_k[5,3,1]",
                                  "beta_k[6,1,1]", "beta_k[6,2,1]","beta_k[6,3,1]",
                                  "beta_k[7,1,1]", "beta_k[7,2,1]","beta_k[7,3,1]"), nrow = M, ncol= K)
beta_treatment_name_vec <- matrix(c("beta_k[1,1,2]", "beta_k[1,2,2]", "beta_k[1,3,2]",
                                    "beta_k[2,1,2]", "beta_k[2,2,2]","beta_k[2,3,2]",
                                    "beta_k[3,1,2]", "beta_k[3,2,2]","beta_k[3,3,2]",
                                    "beta_k[4,1,2]", "beta_k[4,2,2]", "beta_k[4,3,2]",
                                    "beta_k[5,1,2]", "beta_k[5,2,2]","beta_k[5,3,2]",
                                    "beta_k[6,1,2]", "beta_k[6,2,2]","beta_k[6,3,2]",
                                    "beta_k[7,1,2]", "beta_k[7,2,2]", "beta_k[7,3,2]"), nrow = M, ncol= K)

for (k in 1:K){
  for(s in 1:S){

    control_logmean_1 <- (codafit_stan_draws[s,mu_k_name_vec[k]])
    control_logmean_2 <- (codafit_stan_draws[s,mu_k_name_vec[k+7]])
    treatment_logmean_1 <- (codafit_stan_draws[s,mu_k_name_vec[k]] + codafit_stan_draws[s,tau_k_name_vec[k]])
    treatment_logmean_2 <- (codafit_stan_draws[s,mu_k_name_vec[k+7]] + codafit_stan_draws[s,tau_k_name_vec[k+7]])
    
    control_logsd_1 <- exp(codafit_stan_draws[s,sigma_control_k_name_vec[k]])
    treatment_logsd_1 <- exp(codafit_stan_draws[s,sigma_control_k_name_vec[k]] + codafit_stan_draws[s,sigma_TE_k_name_vec[k]])
    control_logsd_2 <- exp(codafit_stan_draws[s,sigma_control_k_name_vec[k+7]])
    treatment_logsd_2 <- exp(codafit_stan_draws[s,sigma_control_k_name_vec[k+7]] + codafit_stan_draws[s,sigma_TE_k_name_vec[k+7]])
    
    control_logit_params <- codafit_stan_draws[s,beta_control_name_vec[,k]]
    control_probs <- convert_logit_2_prob(control_logit_params)
    treatment_logit_params <- codafit_stan_draws[s,beta_treatment_name_vec[,k]]
    treatment_probs <- convert_logit_2_prob(control_logit_params + treatment_logit_params)

    # now get the N loop from the old graphics 
    for(n in 1:N){
 
      u <- quantile_vec[n]
      if(u < control_probs[1]){ quantile_value_output_control[k,s,n]  <-  - qlnorm(1-(u/control_probs[1]), meanlog = control_logmean_1, sdlog = control_logsd_1, lower.tail = TRUE, log.p = FALSE) }
      if(control_probs[1] < u & u < (control_probs[1]+control_probs[2])){ quantile_value_output_control[k,s,n]  <-  0 }
      if(u > (control_probs[1]+control_probs[2])){ quantile_value_output_control[k,s,n]  <-  qlnorm(((u-(control_probs[1]+control_probs[2]))/(1-control_probs[1]-control_probs[2])),  meanlog = control_logmean_2, sdlog = control_logsd_2, lower.tail = TRUE, log.p = FALSE)  }
      
      if(u < treatment_probs[1]){ quantile_value_output_treatment[k,s,n]  <-  -qlnorm(1-(u/treatment_probs[1]), meanlog = treatment_logmean_1, sdlog = treatment_logsd_1, lower.tail = TRUE, log.p = FALSE) }
      if(treatment_probs[1] < u & u < (treatment_probs[1]+treatment_probs[2])){ quantile_value_output_treatment[k,s,n]  <-  0 }
      if(u > (treatment_probs[1]+treatment_probs[2])){ quantile_value_output_treatment[k,s,n]  <-  qlnorm(((u-(treatment_probs[1]+treatment_probs[2]))/(1-treatment_probs[1]-treatment_probs[2])), meanlog = treatment_logmean_2, sdlog = treatment_logsd_2, lower.tail = TRUE, log.p = FALSE)  }
      
      
    } # close forloop indexed by n 
    
    
    
  } # closes forloop indexed by s 
  
  
} # closes forloop indexed by K 



quantile_value_output_difference <-  quantile_value_output_treatment -  quantile_value_output_control

quantile_value_output_control_mean <- matrix(NA, nrow=K, ncol=N)
quantile_value_output_treatment_mean <- matrix(NA, nrow=K, ncol=N)
profit_quantile_value_output_difference_mean <- matrix(NA, nrow=K, ncol=N)
profit_recovered_quantile_differences <- array( rep(NA, K*N*5), dim = c(5,N,K))  # due to the dim of "get posterior intervals function" output 
profit_recovered_quantile_control  <- array( rep(NA, K*N*5), dim = c(5,N,K)) 
profit_recovered_quantile_treatment  <- array( rep(NA, K*N*5), dim =c(5,N,K)) 

for(k in 1:K){
  quantile_value_output_control_mean[k,] <- apply(quantile_value_output_control[k,,],2,mean)
  quantile_value_output_treatment_mean[k,] <- apply(quantile_value_output_treatment[k,,],2,mean)
  profit_quantile_value_output_difference_mean[k,] <- apply(quantile_value_output_difference[k,,],2,mean)
  profit_recovered_quantile_differences[,,k] <- apply(quantile_value_output_difference[k,,],2,get_posterior_intervals_function)
  profit_recovered_quantile_control[,,k] <- apply(quantile_value_output_control[k,,],2,get_posterior_intervals_function)
  profit_recovered_quantile_treatment[,,k] <- apply(quantile_value_output_treatment[k,,],2,get_posterior_intervals_function)
  
}


# tabular bit

median_quantile_effects_profit_no_pooling <- round(t(profit_recovered_quantile_differences[3,,]),1)
lower_ci_quantile_effects_profit_no_pooling <- round(t(profit_recovered_quantile_differences[1,,]),1)
upper_ci_quantile_effects_profit_no_pooling <- round(t(profit_recovered_quantile_differences[5,,]),1) 
ci_quantile_effects_profit_no_pooling <- paste0("(", lower_ci_quantile_effects_profit_no_pooling[,], ",", upper_ci_quantile_effects_profit_no_pooling[,], ")" )
ci_quantile_effects_profit_no_pooling <- matrix(ci_quantile_effects_profit_no_pooling, nrow=7, ncol=10)
ci_quantile_effects_profit_no_pooling <- data.frame(country,ci_quantile_effects_profit_no_pooling )
median_quantile_effects_profit_no_pooling <- data.frame(country, median_quantile_effects_profit_no_pooling)
rbindlist(list(median_quantile_effects_profit_no_pooling, ci_quantile_effects_profit_no_pooling))[order(country)]
no_pooling_results_profit <- rbindlist(list(median_quantile_effects_profit_no_pooling, ci_quantile_effects_profit_no_pooling))[order(country)]
colnames(no_pooling_results_profit) <- c("Country", "5th", "15th", "25th", "35th", "45th", "55th", "65th", "75th", "85th", "95th")
saveRDS(no_pooling_results_profit, "output/no_pooling_table_profit.rds")
testy <- readRDS("output/no_pooling_table_profit.rds")

# we need to get the partial for comparison

codafit_stan_draws <- codafit_stan_draws_profit
S <- dim(codafit_stan_draws)[1]
M <- 3 # mixture components
K <- as.numeric(length(unique(data$site)))
component_probs <- matrix(NA,nrow=S,ncol=M)
quantile_vec <- seq(0.05, 0.95,0.1)
N <- length(quantile_vec)
quantile_value_output_control <- array( rep(NA, N*K*S), dim= c(K,S,N))
quantile_value_output_treatment <- array( rep(NA, N*K*S), dim= c(K,S,N))
quantile_value_output_difference <- array( rep(NA, N*K*S), dim= c(K,S,N))
mu_k_name_vec <- grep("mu_k", colnames(codafit_stan_draws_profit), value = TRUE)
tau_k_name_vec <- grep("tau_k", colnames(codafit_stan_draws_profit), value = TRUE)
sigma_control_k_name_vec <- grep("sigma_control_k", colnames(codafit_stan_draws_profit), value = TRUE)
sigma_TE_k_name_vec <- grep("sigma_TE_k", colnames(codafit_stan_draws_profit), value = TRUE)
beta_control_name_vec <- matrix(c("beta_k[1,1,1]", "beta_k[1,2,1]","beta_k[1,3,1]",
                                  "beta_k[2,1,1]", "beta_k[2,2,1]","beta_k[2,3,1]",
                                  "beta_k[3,1,1]", "beta_k[3,2,1]","beta_k[3,3,1]",
                                  "beta_k[4,1,1]", "beta_k[4,2,1]","beta_k[4,3,1]",
                                  "beta_k[5,1,1]", "beta_k[5,2,1]", "beta_k[5,3,1]",
                                  "beta_k[6,1,1]", "beta_k[6,2,1]","beta_k[6,3,1]",
                                  "beta_k[7,1,1]", "beta_k[7,2,1]","beta_k[7,3,1]"), nrow = M, ncol= K)
beta_treatment_name_vec <- matrix(c("beta_k[1,1,2]", "beta_k[1,2,2]", "beta_k[1,3,2]",
                                    "beta_k[2,1,2]", "beta_k[2,2,2]","beta_k[2,3,2]",
                                    "beta_k[3,1,2]", "beta_k[3,2,2]","beta_k[3,3,2]",
                                    "beta_k[4,1,2]", "beta_k[4,2,2]", "beta_k[4,3,2]",
                                    "beta_k[5,1,2]", "beta_k[5,2,2]","beta_k[5,3,2]",
                                    "beta_k[6,1,2]", "beta_k[6,2,2]","beta_k[6,3,2]",
                                    "beta_k[7,1,2]", "beta_k[7,2,2]", "beta_k[7,3,2]"), nrow = M, ncol= K)

for (k in 1:K){
  for(s in 1:S){
    
    control_logmean_1 <- (codafit_stan_draws[s,mu_k_name_vec[k]])
    control_logmean_2 <- (codafit_stan_draws[s,mu_k_name_vec[k+7]])
    treatment_logmean_1 <- (codafit_stan_draws[s,mu_k_name_vec[k]] + codafit_stan_draws[s,tau_k_name_vec[k]])
    treatment_logmean_2 <- (codafit_stan_draws[s,mu_k_name_vec[k+7]] + codafit_stan_draws[s,tau_k_name_vec[k+7]])
    
    control_logsd_1 <- exp(codafit_stan_draws[s,sigma_control_k_name_vec[k]])
    treatment_logsd_1 <- exp(codafit_stan_draws[s,sigma_control_k_name_vec[k]] + codafit_stan_draws[s,sigma_TE_k_name_vec[k]])
    control_logsd_2 <- exp(codafit_stan_draws[s,sigma_control_k_name_vec[k+7]])
    treatment_logsd_2 <- exp(codafit_stan_draws[s,sigma_control_k_name_vec[k+7]] + codafit_stan_draws[s,sigma_TE_k_name_vec[k+7]])
    
    control_logit_params <- codafit_stan_draws[s,beta_control_name_vec[,k]]
    control_probs <- convert_logit_2_prob(control_logit_params)
    treatment_logit_params <- codafit_stan_draws[s,beta_treatment_name_vec[,k]]
    treatment_probs <- convert_logit_2_prob(control_logit_params + treatment_logit_params)
    
    # now get the N loop from the old graphics 
    for(n in 1:N){
      
      u <- quantile_vec[n]
      if(u < control_probs[1]){ quantile_value_output_control[k,s,n]  <-  - qlnorm(1-(u/control_probs[1]), meanlog = control_logmean_1, sdlog = control_logsd_1, lower.tail = TRUE, log.p = FALSE) }
      if(control_probs[1] < u & u < (control_probs[1]+control_probs[2])){ quantile_value_output_control[k,s,n]  <-  0 }
      if(u > (control_probs[1]+control_probs[2])){ quantile_value_output_control[k,s,n]  <-  qlnorm(((u-(control_probs[1]+control_probs[2]))/(1-control_probs[1]-control_probs[2])),  meanlog = control_logmean_2, sdlog = control_logsd_2, lower.tail = TRUE, log.p = FALSE)  }
      
      if(u < treatment_probs[1]){ quantile_value_output_treatment[k,s,n]  <-  -qlnorm(1-(u/treatment_probs[1]), meanlog = treatment_logmean_1, sdlog = treatment_logsd_1, lower.tail = TRUE, log.p = FALSE) }
      if(treatment_probs[1] < u & u < (treatment_probs[1]+treatment_probs[2])){ quantile_value_output_treatment[k,s,n]  <-  0 }
      if(u > (treatment_probs[1]+treatment_probs[2])){ quantile_value_output_treatment[k,s,n]  <-  qlnorm(((u-(treatment_probs[1]+treatment_probs[2]))/(1-treatment_probs[1]-treatment_probs[2])), meanlog = treatment_logmean_2, sdlog = treatment_logsd_2, lower.tail = TRUE, log.p = FALSE)  }
      
      
    } # close forloop indexed by n 
    
    
    
  } # closes forloop indexed by s 
  
  
} # closes forloop indexed by K 



quantile_value_output_difference <-  quantile_value_output_treatment -  quantile_value_output_control

quantile_value_output_control_mean <- matrix(NA, nrow=K, ncol=N)
quantile_value_output_treatment_mean <- matrix(NA, nrow=K, ncol=N)
profit_quantile_value_output_difference_mean <- matrix(NA, nrow=K, ncol=N)
profit_recovered_quantile_differences <- array( rep(NA, K*N*5), dim = c(5,N,K))  # due to the dim of "get posterior intervals function" output 
profit_recovered_quantile_control  <- array( rep(NA, K*N*5), dim = c(5,N,K)) 
profit_recovered_quantile_treatment  <- array( rep(NA, K*N*5), dim =c(5,N,K)) 

for(k in 1:K){
  quantile_value_output_control_mean[k,] <- apply(quantile_value_output_control[k,,],2,mean)
  quantile_value_output_treatment_mean[k,] <- apply(quantile_value_output_treatment[k,,],2,mean)
  profit_quantile_value_output_difference_mean[k,] <- apply(quantile_value_output_difference[k,,],2,mean)
  profit_recovered_quantile_differences[,,k] <- apply(quantile_value_output_difference[k,,],2,get_posterior_intervals_function)
  profit_recovered_quantile_control[,,k] <- apply(quantile_value_output_control[k,,],2,get_posterior_intervals_function)
  profit_recovered_quantile_treatment[,,k] <- apply(quantile_value_output_treatment[k,,],2,get_posterior_intervals_function)
  
}


# tabular bit

median_quantile_effects_profit_partial_pooling <- round(t(profit_recovered_quantile_differences[3,,]),1)
lower_ci_quantile_effects_profit_partial_pooling <- round(t(profit_recovered_quantile_differences[1,,]),1)
upper_ci_quantile_effects_profit_partial_pooling <- round(t(profit_recovered_quantile_differences[5,,]),1) 
ci_quantile_effects_profit_partial_pooling <- paste0("(", lower_ci_quantile_effects_profit_partial_pooling[,], ",", upper_ci_quantile_effects_profit_partial_pooling[,], ")" )
ci_quantile_effects_profit_partial_pooling <- matrix(ci_quantile_effects_profit_partial_pooling, nrow=7, ncol=10)
ci_quantile_effects_profit_partial_pooling <- data.frame(country,ci_quantile_effects_profit_partial_pooling )
median_quantile_effects_profit_partial_pooling <- data.frame(country, median_quantile_effects_profit_partial_pooling)
rbindlist(list(median_quantile_effects_profit_partial_pooling, ci_quantile_effects_profit_partial_pooling))[order(country)]
partial_pooling_results_profit <- rbindlist(list(median_quantile_effects_profit_partial_pooling, ci_quantile_effects_profit_partial_pooling))[order(country)]
colnames(partial_pooling_results_profit) <- c("Country", "5th", "15th", "25th", "35th", "45th", "55th", "65th", "75th", "85th", "95th")
saveRDS(partial_pooling_results_profit, "output/partial_pooling_table_profit.rds")
testy <- readRDS("output/partial_pooling_table_profit.rds")




