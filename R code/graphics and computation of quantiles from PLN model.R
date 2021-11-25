# graphics and quantile computation for the lognormals 
# Rachael Meager
# First Version: April 2018 

### Notes ###

# Analysis needs to be done before you can run this, obviously
# I sincerely apologise for the atrocious indenting here. I would like to especially apologise to Jonathan Huggins who taught me better.

### Preliminaries and Data Intake ###

# clear the workspace to avoid gremlins and past globals from past irresponsible scripts
# but we can't do this if the masterfile is being used to run the script, so we check that first:
if(exists("masterfile_run") == "FALSE"){
  rm(list = ls())
}


# you may want to use this, but then again, you may not: setwd("/Users/rachaelmeager/Dropbox/research work/MIT IMPRINT/Research work/aggregating distributional effects/")


# Load data

library(xtable)
library(MCMCpack)
install.packages("NormalLaplace")
library(NormalLaplace)

codafit_stan_draws_profit <- readRDS("output/tailored_hierarchical_pdf_microcredit_output_PLN_4000_iters_codafit.RDS")


### FUNCTIONS I WILL NEED

qPLN_sim <- function(u, alpha, nu, tau, S){
  n <- S
  E1 <- rexp(n, rate = 1)
  Z <- rnorm(n, mean = 0, sd = 1)
  Y <- nu + tau*Z + E1/alpha 
  x <- exp(Y) 
  return(quantile(x,u))
  
}


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

# make draws

codafit_stan_draws <- codafit_stan_draws_profit
S <- dim(codafit_stan_draws)[1]
M <- 3 # mixture components
component_probs <- matrix(NA,nrow=S,ncol=M)
quantile_vec <- seq(0.05, 0.95,0.1)
N <- length(quantile_vec)
quantile_value_output_control <- matrix(NA, nrow=S, ncol=N)
quantile_value_output_treatment <- matrix(NA, nrow=S, ncol=N)
quantile_value_output_difference <- matrix(NA, nrow=S, ncol=N)


for(s in 1:S){

  control_alpha_1 <- exp(codafit_stan_draws[s,"log_alpha[1]"])
  treatment_alpha_1 <- exp(codafit_stan_draws[s,"log_alpha[1]"] + codafit_stan_draws[s,"log_alpha_TE[1]"])
  control_nu_1 <- (codafit_stan_draws[s,"nu[1]"])
  treatment_nu_1 <- (codafit_stan_draws[s,"nu[1]"] + codafit_stan_draws[s,"nu_TE[1]"])
  control_tau_1 <- exp(codafit_stan_draws[s,"log_tau[1]"])
  treatment_tau_1 <- exp(codafit_stan_draws[s,"log_tau[1]"] + codafit_stan_draws[s,"log_tau_TE[1]"])
  
  control_alpha_2 <- exp(codafit_stan_draws[s,"log_alpha[2]"])
  treatment_alpha_2 <- exp(codafit_stan_draws[s,"log_alpha[2]"] + codafit_stan_draws[s,"log_alpha_TE[2]"])
  control_nu_2 <- (codafit_stan_draws[s,"nu[2]"])
  treatment_nu_2 <- (codafit_stan_draws[s,"nu[2]"] + codafit_stan_draws[s,"nu_TE[2]"])
  control_tau_2 <- exp(codafit_stan_draws[s,"log_tau[2]"])
  treatment_tau_2 <- exp(codafit_stan_draws[s,"log_tau[2]"] + codafit_stan_draws[s,"log_tau_TE[2]"])
  
  
  
  control_logit_params <- codafit_stan_draws[s,c("beta_full[1,1]", "beta_full[2,1]", "beta_full[3,1]")]
  control_probs <- convert_logit_2_prob(control_logit_params)
  treatment_logit_params <- codafit_stan_draws[s,c("beta_full[1,2]", "beta_full[2,2]", "beta_full[3,2]")]
  treatment_probs <- convert_logit_2_prob(control_logit_params + treatment_logit_params)
  
  
  for(n in 1:N){

    u <- quantile_vec[n]
    if(u < control_probs[1]){ quantile_value_output_control[s,n]  <-  - qPLN_sim(u = 1-(u/control_probs[1]), alpha = control_alpha_1, nu = control_nu_1, tau = control_tau_1, S=10000) }
    if(control_probs[1] < u & u < (control_probs[1]+control_probs[2])){ quantile_value_output_control[s,n]  <-  0 }
    if(u > (control_probs[1]+control_probs[2])){ quantile_value_output_control[s,n]  <-  qPLN_sim(u = ((u-(control_probs[1]+control_probs[2]))/(1-control_probs[1]-control_probs[2])), alpha = control_alpha_2, nu = control_nu_2, tau = control_tau_2, S=10000)  }
    
    if(u < treatment_probs[1]){ quantile_value_output_treatment[s,n]  <-  - qPLN_sim(u = 1-(u/treatment_probs[1]), alpha = treatment_alpha_1, nu = treatment_nu_1, tau = treatment_tau_1, S=10000) }
    if(treatment_probs[1] < u & u < (treatment_probs[1]+treatment_probs[2])){ quantile_value_output_treatment[s,n]  <-  0 }
    if(u > (treatment_probs[1]+treatment_probs[2])){ quantile_value_output_treatment[s,n]  <-  qPLN_sim(u = ((u-(treatment_probs[1]+treatment_probs[2]))/(1-treatment_probs[1]-treatment_probs[2])), alpha = treatment_alpha_2, nu = treatment_nu_2, tau = treatment_tau_2, S=10000) }
    
    
  } # close forloop indexed by n 
  
} # closes forloop indexed by s 


quantile_value_output_difference <-  quantile_value_output_treatment -  quantile_value_output_control


quantile_value_output_control_mean <- apply(quantile_value_output_control,2,mean)
quantile_value_output_treatment_mean <- apply(quantile_value_output_treatment,2,mean)
profit_quantile_value_output_difference_mean <- apply(quantile_value_output_difference,2,mean)

quantile_value_output_difference <- quantile_value_output_difference[is.finite(rowSums(quantile_value_output_difference)),]
quantile_value_output_difference <- quantile_value_output_difference[complete.cases(quantile_value_output_difference),]
profit_recovered_quantile_differences <- apply(quantile_value_output_difference,2,get_posterior_intervals_function)
profit_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)
profit_recovered_quantile_treatment <- apply(quantile_value_output_treatment,2,get_posterior_intervals_function)



### NOW THE GRAPHICS ###

# PROFIT # 
# make ribbon plot for posterior 
quantiles_list <- quantile_vec
posterior_beta_data <- t(rbind(profit_recovered_quantile_differences, profit_quantile_value_output_difference_mean  ))
posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")


library(ggplot2)
fig_scale = 0.4
posterior_beta_data <- sign(posterior_beta_data)*(abs(posterior_beta_data))^.1
posterior_beta_data_plot <- ggplot(posterior_beta_data, aes(posterior_beta_data$quantiles_list))
pdf("output/posterior_parent_quantile_TEs_profit_PLN.pdf", width=fig_scale*6.5, height=fig_scale*6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "red", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "red", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"50%"]), color = "maroon", size = fig_scale*1.5) +
  ggtitle("Posterior quantile effects on profit from the PLN model") +
  theme(plot.title = element_text(size = fig_scale*16)) + #xlim(0.05,0.95) +  ylim(-100,(6^64)) +
  xlab("Quantiles") + ylab("Quantile treatment effect")+
  theme(axis.text = element_text(size=fig_scale*14)) +  theme(axis.title.y = element_text(size = fig_scale*14)) +  theme(axis.title.x = element_text(size = fig_scale*14))
dev.off()


