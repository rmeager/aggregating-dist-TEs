
# big tables for profit and consumption
# Rachael Meager
# First Version: nov 2018 
# This version: August 2020 

### Notes ###

# Analysis AND GRAPHICS needs to be done before you can run this, it exploits the dataframe structure I used in the graphics 

### Preliminaries and Data Intake ###

installation_needed  <- FALSE
loading_needed <- TRUE
package_list <- c('ggplot2', 'rstan','reshape','reshape2','coda','xtable', 'dplyr', 'Runuran', 'testthat',
                  "MCMCpack", "gtools", 'gPdtest', 'fBasics',"PtProcess", "VGAM", "MASS","quantreg",
                  "boot", "gridExtra", "stargazer", "data.table")
if(installation_needed){install.packages(package_list, repos='http://cran.us.r-project.org')}
if(loading_needed){lapply(package_list, require, character.only = TRUE)}


### PROFIT ####

# Load partial and full pooling bits 
partial_pooling_beta_1 <- readRDS("output/partial_pooling_beta_1_table_profit.rds")
full_pooling_beta_1 <- readRDS("output/full_pooling_beta_1_table_profit.rds")
no_pooling_profit <- readRDS("output/no_pooling_table_profit.rds")
partial_pooling_profit <- readRDS("output/partial_pooling_table_profit.rds")

# stitch them all together
full_pooling_beta_1 <- cbind(c("Average", "Average"), full_pooling_beta_1)
colnames(full_pooling_beta_1)[1] <- "Country"
partial_pooling_beta_1 <- cbind(c("Average", "Average"), partial_pooling_beta_1)
colnames(partial_pooling_beta_1)[1] <- "Country"
profit_table <- rbind(no_pooling_profit, partial_pooling_profit,partial_pooling_beta_1,full_pooling_beta_1, fill = TRUE)

sink("output/table_profit_all.txt")
stargazer(profit_table, summary = FALSE)
dev.off


### COSUMPTION ###
load("output/microcredit_consumption_lognormal_tailored_hierarchical_pdf_output_5000_iters.RData")
stan_fit_table_consumption <- stan_fit_table
stan_fit_consumption <- stan_fit
codafit_stan_draws_consumption <- as.matrix(stan2coda(stan_fit_consumption))
data_split_consumption <- data_split
data_consumption <- data

# quantiles_list <- quantile_vec
# posterior_beta_data <- t(rbind(consumption_recovered_quantile_differences, consumption_quantile_value_output_difference_mean  ))
# posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
# colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")

country <- c( "Mexico", "Mongolia", "Bosnia", "India", "Morocco")


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



codafit_stan_draws <- codafit_stan_draws_consumption
S <- dim(codafit_stan_draws)[1]
M <- 2 # mixture components
K <- as.numeric(length(unique(data$site)))
component_probs <- matrix(NA,nrow=S,ncol=M)
quantile_vec <- seq(0.05, 0.95,0.1)
N <- length(quantile_vec)
quantile_value_output_control <- array( rep(NA, N*K*S), dim= c(K,S,N))
quantile_value_output_treatment <- array( rep(NA, N*K*S), dim= c(K,S,N))
quantile_value_output_difference <- array( rep(NA, N*K*S), dim= c(K,S,N))
mu_k_name_vec <- grep("mu_k", colnames(codafit_stan_draws_consumption), value = TRUE)
tau_k_name_vec <- grep("tau_k", colnames(codafit_stan_draws_consumption), value = TRUE)
sigma_control_k_name_vec <- grep("sigma_control_k", colnames(codafit_stan_draws_consumption), value = TRUE)
sigma_TE_k_name_vec <- grep("sigma_TE_k", colnames(codafit_stan_draws_consumption), value = TRUE)
beta_control_name_vec <- matrix(c("beta_k[1,1,1]", "beta_k[1,2,1]",
                                  "beta_k[2,1,1]", "beta_k[2,2,1]",
                                  "beta_k[3,1,1]", "beta_k[3,2,1]",
                                  "beta_k[4,1,1]", "beta_k[4,2,1]",
                                  "beta_k[5,1,1]", "beta_k[5,2,1]"), nrow = M, ncol= K)
beta_treatment_name_vec <- matrix(c("beta_k[1,1,2]", "beta_k[1,2,2]",
                                    "beta_k[2,1,2]", "beta_k[2,2,2]",
                                    "beta_k[3,1,2]", "beta_k[3,2,2]",
                                    "beta_k[4,1,2]", "beta_k[4,2,2]",
                                    "beta_k[5,1,2]", "beta_k[5,2,2]"), nrow = M, ncol= K)

for (k in 1:K){
  for(s in 1:S){
    
    control_logmean <- (codafit_stan_draws[s,mu_k_name_vec[k]])
    treatment_logmean <- (codafit_stan_draws[s,mu_k_name_vec[k]] + codafit_stan_draws[s,tau_k_name_vec[k]])
    
    control_logsd <- exp(codafit_stan_draws[s,sigma_control_k_name_vec[k]])
    treatment_logsd <- exp(codafit_stan_draws[s,sigma_control_k_name_vec[k]] + codafit_stan_draws[s,sigma_TE_k_name_vec[k]])

    control_logit_params <- codafit_stan_draws[s,beta_control_name_vec[,k]]
    control_probs <- convert_logit_2_prob(control_logit_params)
    treatment_logit_params <- codafit_stan_draws[s,beta_treatment_name_vec[,k]]
    treatment_probs <- convert_logit_2_prob(control_logit_params + treatment_logit_params)
    
    # a question for me: does this structure REALLY propagate like the Gaussian structure does? 
    # One might ask "Should I have been doing everything on the log scale?" but this will not save you in a mixture model. 
    
    # now get the N loop from the old graphics 
    for(n in 1:N){
      
      u <- quantile_vec[n]
      if(u < control_probs[1]){ quantile_value_output_control[k,s,n] <-  0 }
      if(u > control_probs[1]){ quantile_value_output_control[k,s,n] <-  qlnorm(((u-control_probs[1])/control_probs[2]), meanlog = control_logmean, sdlog = control_logsd, lower.tail = TRUE, log.p = FALSE) }
      if(u < treatment_probs[1]){ quantile_value_output_treatment[k,s,n] <-  0 }
      if(u > treatment_probs[1]){ quantile_value_output_treatment[k,s,n] <-  qlnorm(((u-treatment_probs[1])/treatment_probs[2]), meanlog = treatment_logmean, sdlog = treatment_logsd, lower.tail = TRUE, log.p = FALSE)  }

    } # close forloop indexed by n 
    
  } # closes forloop indexed by s 
  
  
} # closes forloop indexed by K 



quantile_value_output_difference <-  quantile_value_output_treatment -  quantile_value_output_control

quantile_value_output_control_mean <- matrix(NA, nrow=K, ncol=N)
quantile_value_output_treatment_mean <- matrix(NA, nrow=K, ncol=N)
consumption_quantile_value_output_difference_mean <- matrix(NA, nrow=K, ncol=N)
consumption_recovered_quantile_differences <- array( rep(NA, K*N*5), dim = c(5,N,K))  # due to the dim of "get posterior intervals function" output 
consumption_recovered_quantile_control  <- array( rep(NA, K*N*5), dim = c(5,N,K)) 
consumption_recovered_quantile_treatment  <- array( rep(NA, K*N*5), dim =c(5,N,K)) 

for(k in 1:K){
  quantile_value_output_control_mean[k,] <- apply(quantile_value_output_control[k,,],2,mean)
  quantile_value_output_treatment_mean[k,] <- apply(quantile_value_output_treatment[k,,],2,mean)
  consumption_quantile_value_output_difference_mean[k,] <- apply(quantile_value_output_difference[k,,],2,mean)
  consumption_recovered_quantile_differences[,,k] <- apply(quantile_value_output_difference[k,,],2,get_posterior_intervals_function)
  consumption_recovered_quantile_control[,,k] <- apply(quantile_value_output_control[k,,],2,get_posterior_intervals_function)
  consumption_recovered_quantile_treatment[,,k] <- apply(quantile_value_output_treatment[k,,],2,get_posterior_intervals_function)
  
}



consumption_recovered_quantile_differences_flat <- t(consumption_recovered_quantile_differences[,,1])
for(k in 2:K){
  consumption_recovered_quantile_differences_flat <- rbind(consumption_recovered_quantile_differences_flat, t(consumption_recovered_quantile_differences[,,k])) }
posterior_k_beta_data <- data.frame(consumption_recovered_quantile_differences_flat)
quantiles_repeated <- rep(c(0.05, .15, .25, .35,.45,.55,.65,.75,.85,.95),K)
posterior_k_beta_data <- cbind(posterior_k_beta_data, quantiles_repeated )
countries_repeated <- c(rep("Mexico", 10),
                        rep("Mongolia", 10),
                        rep("Bosnia", 10),
                        rep("India", 10),
                        rep("Morocco", 10))
countries_repeated <- rep(seq(1,K,1), each=10)

posterior_k_beta_data <- cbind(posterior_k_beta_data, countries_repeated)

posterior_k_beta_data$countries_repeated <- factor(posterior_k_beta_data$countries_repeated, labels = country)

colnames(posterior_k_beta_data) <- c( "2.5%", "25%", "50%", "75%", "97.5%", "quantiles_repeated", "countries_repeated")


# now for a tabular bit
median_quantile_effects_consumption_partial_pooling <- round(t(consumption_recovered_quantile_differences[3,,]),1)
lower_ci_quantile_effects_consumption_partial_pooling <- round(t(consumption_recovered_quantile_differences[1,,]),1)
upper_ci_quantile_effects_consumption_partial_pooling <- round(t(consumption_recovered_quantile_differences[5,,]),1) 
ci_quantile_effects_consumption_partial_pooling <- paste0("(", lower_ci_quantile_effects_consumption_partial_pooling[,], ",", upper_ci_quantile_effects_consumption_partial_pooling[,], ")" )
ci_quantile_effects_consumption_partial_pooling <- matrix(ci_quantile_effects_consumption_partial_pooling, nrow=K, ncol=10)
ci_quantile_effects_consumption_partial_pooling <- data.frame(country,ci_quantile_effects_consumption_partial_pooling )
median_quantile_effects_consumption_partial_pooling <- data.frame(country, median_quantile_effects_consumption_partial_pooling)
partial_pooling_results_consumption <- rbindlist(list(median_quantile_effects_consumption_partial_pooling, ci_quantile_effects_consumption_partial_pooling))[order(country)]
colnames(partial_pooling_results_consumption) <- c("Country", "5th", "15th", "25th", "35th", "45th", "55th", "65th", "75th", "85th", "95th")



partial_pooling_consumption_avg <- readRDS("output/partial_pooling_consumption_avg.rds")
partial_pooling_consumption_avg <- cbind(c("Average","."), partial_pooling_consumption_avg)
colnames(partial_pooling_consumption_avg) <- c("Country", "5th", "15th", "25th", "35th", "45th", "55th", "65th", "75th", "85th", "95th")

partial_pooling_results_consumption <- rbind(partial_pooling_results_consumption, partial_pooling_consumption_avg)



# now we will do the no pooling and full pooling for consumption

# denote quantiles of interest
quantiles_list <- seq(0.05, 0.95,0.10)
N <- as.numeric(length(quantiles_list))
K <- as.numeric(length(unique(site)))

# now quantile reg!! 
base_quantile_picker <- function(quantile_reg){
  quantile_reg$coef[1,1]
}

base_quantile_se_picker <- function(quantile_reg){
  quantile_reg$coef[1,2]
}
te_picker <- function(quantile_reg){
  quantile_reg$coef[2,1]
}
se_picker <- function(quantile_reg){
  quantile_reg$coef[2,2]
}  

quantile_reg <- list()
quantile_summary <- list()
y_0 <- matrix(NA,K,N)
y_0_se <- matrix(NA,K,N)
quantile_effects <- matrix(NA,K,N)
quantile_effects_ses <- matrix(NA,K,N)

for(k in 1:K){
  quantile_reg[[k]] <- rq(data$consumption[(data$site==k)&(data$treatment==0)] ~ 1, tau = quantiles_list)
  
  quantile_summary[[k]] <- summary(quantile_reg[[k]],se = "iid")
  
  y_0[k,] <- unlist(lapply(quantile_summary[[k]],base_quantile_picker))
  y_0_se[k,] <- unlist(lapply(quantile_summary[[k]], base_quantile_se_picker))
  
} # closes forloop indexed by k

quantile_reg <- list()
quantile_summary <- list()
y_1 <- matrix(NA,K,N)
y_1_se <- matrix(NA,K,N)

for(k in 1:K){
  quantile_reg[[k]] <- rq(data$consumption[(data$site==k)&(data$treatment==1)] ~ 1, tau = quantiles_list)
  
  quantile_summary[[k]] <- summary(quantile_reg[[k]],se = "iid")
  
  y_1[k,] <- unlist(lapply(quantile_summary[[k]],base_quantile_picker))
  y_1_se[k,] <- unlist(lapply(quantile_summary[[k]], base_quantile_se_picker))
  
} # closes forloop indexed by k

quantile_reg <- list()
quantile_summary <- list()
beta_1 <- matrix(NA,K,N)
beta_1_se <- matrix(NA,K,N)

for(k in 1:K){
  quantile_reg[[k]] <- rq(data$consumption[(data$site==k)] ~ data$treatment[(data$site==k)], tau = quantiles_list)
  
  quantile_summary[[k]] <- summary(quantile_reg[[k]],se = "iid")
  
  beta_1[k,] <- unlist(lapply(quantile_summary[[k]],te_picker))
  beta_1_se[k,] <- unlist(lapply(quantile_summary[[k]], se_picker))
  
} # closes forloop indexed by k


# testing for my own edification
quantile_reg <- list()
quantile_summary <- list()


for(k in 1:K){
  quantile_reg[[k]] <- rq(data$consumption[(data$site==k)] ~ data$treatment[(data$site==k)] , tau = quantiles_list)
  
  quantile_summary[[k]] <- summary(quantile_reg[[k]],se = "boot")
  quantile_effects[k,] <- unlist(lapply(quantile_summary[[k]],te_picker))
  quantile_effects_ses[k,] <- unlist(lapply(quantile_summary[[k]], se_picker))
  
  
} # closes forloop indexed by k



print(y_0)
print(y_1)
print(y_0_se)
print(y_1_se)
print(y_1-y_0)


# now print out the "no pooling" results


quantile_effects_consumption_no_pooling <- round(quantile_effects,1)
lower_ci_quantile_effects_consumption_no_pooling <- round(quantile_effects - 1.96*quantile_effects_ses,1)
upper_ci_quantile_effects_consumption_no_pooling <- round(quantile_effects + 1.96*quantile_effects_ses,1) 
ci_quantile_effects_consumption_no_pooling <- paste0("(", lower_ci_quantile_effects_consumption_no_pooling[,], ",", upper_ci_quantile_effects_consumption_no_pooling[,], ")" )
ci_quantile_effects_consumption_no_pooling <- matrix(ci_quantile_effects_consumption_no_pooling, nrow=K, ncol=10)
ci_quantile_effects_consumption_no_pooling <- data.frame(country,ci_quantile_effects_consumption_no_pooling )
quantile_effects_consumption_no_pooling <- data.frame(country,quantile_effects_consumption_no_pooling)
no_pooling_results_consumption <- rbindlist(list(quantile_effects_consumption_no_pooling, ci_quantile_effects_consumption_no_pooling))[order(country)]
colnames(no_pooling_results_consumption) <- c("Country", "5th", "15th", "25th", "35th", "45th", "55th", "65th", "75th", "85th", "95th")


# now try the full pooling


full_pooling <- summary(rq(data$consumption ~ factor(data$site) + data$treatment , tau = quantiles_list))
full_pooling_te <- matrix(NA,1,N)
full_pooling_se <- matrix(NA,1,N)
for(i in 1:N){
  full_pooling_te[,i] <- full_pooling[[i]]$coefficients["data$treatment","Value"]
  full_pooling_se[,i] <- full_pooling[[i]]$coefficients["data$treatment","Std. Error"]
}


quantile_effects_avg_consumption_full_pooling <- round(full_pooling_te,1)
lower_ci_quantile_effects_consumption_full_pooling <- round(full_pooling_te - 1.96*full_pooling_se,1)
upper_ci_quantile_effects_consumption_full_pooling <- round(full_pooling_te + 1.96*full_pooling_se,1)
ci_quantile_effects_consumption_full_pooling <- paste0("(", lower_ci_quantile_effects_consumption_full_pooling[,], ",", upper_ci_quantile_effects_consumption_full_pooling[,], ")" )
full_pooling_results_consumption <- rbind(quantile_effects_avg_consumption_full_pooling, ci_quantile_effects_consumption_full_pooling)
full_pooling_results_consumption <- cbind(c("Average", "Average"), full_pooling_results_consumption)
colnames(full_pooling_results_consumption) <- c( "Country", "5th", "15th", "25th", "35th", "45th", "55th", "65th", "75th", "85th", "95th")


consumption_table <- rbind(no_pooling_results_consumption, partial_pooling_results_consumption, full_pooling_results_consumption)
sink("output/table_consumption_all.txt")
stargazer(consumption_table, summary = FALSE)
dev.off


