# graphics and quantile computation for the lognormals 
# Rachael Meager
# First Version: April 2018 
# This Version: August 2021 

### Notes ###

# Analysis needs to be done before you can run this, obviously
# I sincerely apologise for the atrocious indenting here. I would like to especially apologise to Jonathan Huggins who taught me better.

### Preliminaries and Data Intake ###


install.packages(rstan, coda, ggplot2, stargazer, data.table)
library("rstan", "coda", "ggplot2", "stargazer", "data.table")

# It is on you to set the working directory to the correct location 

# Load data
load("output/microcredit_profit_lognormal_tailored_hierarchical_pdf_output_5000_iters.RData")
stan_fit_table_profit <- stan_fit_table
stan_fit_profit <- stan_fit
codafit_stan_draws_profit <- as.matrix(stan2coda(stan_fit_profit))
nodata_split_profit <- data_split
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

load("output/microcredit_profit_lognormal_tailored_full_pooling_pdf_output_4000_iters.RData")
stan_fit_table_profit_full_pooling <- stan_fit_table
stan_fit_profit_full_pooling <- stan_fit
codafit_stan_draws_profit_full_pooling <- as.matrix(stan2coda(stan_fit_profit_full_pooling))

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

## COMPUTE PARENT QUANTILES (RESULTS FOR FIGURE 1 OF AER VERSION)

# first FOR THE 1 TAIL OUTCOMES

## expenditures

codafit_stan_draws <- codafit_stan_draws_expenditures
S <- dim(codafit_stan_draws)[1]
M <- 2 # mixture components
component_probs <- matrix(NA,nrow=S,ncol=M)
quantile_vec <- seq(0.05, 0.95,0.1)
N <- length(quantile_vec)
quantile_value_output_control <- matrix(NA, nrow=S, ncol=N)
quantile_value_output_treatment <- matrix(NA, nrow=S, ncol=N)
quantile_value_output_difference <- matrix(NA, nrow=S, ncol=N)


for(s in 1:S){

  control_logmean <- (codafit_stan_draws[s,"mu[1]"])
  treatment_logmean <- (codafit_stan_draws[s,"mu[1]"] + codafit_stan_draws[s,"tau[1]"])
  control_logsd <- exp(codafit_stan_draws[s,"sigma_control[1]"])
  treatment_logsd <- exp(codafit_stan_draws[s,"sigma_control[1]"] + codafit_stan_draws[s,"sigma_TE[1]"])
  control_logit_params <- codafit_stan_draws[s,c("beta_full[1,1]", "beta_full[2,1]")]
  control_probs <- convert_logit_2_prob(control_logit_params)
  treatment_logit_params <- codafit_stan_draws[s,c("beta_full[1,2]", "beta_full[2,2]")]
  treatment_probs <- convert_logit_2_prob(control_logit_params + treatment_logit_params)
  

  
  for(n in 1:N){

  u <- quantile_vec[n]
  if(u < control_probs[1]){ quantile_value_output_control[s,n] <-  0 }
  if(u > control_probs[1]){ quantile_value_output_control[s,n] <-  qlnorm(((u-control_probs[1])/control_probs[2]), meanlog = control_logmean, sdlog = control_logsd, lower.tail = TRUE, log.p = FALSE) }
  if(u < treatment_probs[1]){ quantile_value_output_treatment[s,n] <-  0 }
  if(u > treatment_probs[1]){ quantile_value_output_treatment[s,n] <-  qlnorm(((u-treatment_probs[1])/treatment_probs[2]), meanlog = treatment_logmean, sdlog = treatment_logsd, lower.tail = TRUE, log.p = FALSE)  }
  
  } # close forloop indexed by N 

}
quantile_value_output_difference <-  quantile_value_output_treatment -  quantile_value_output_control


quantile_value_output_control_mean <- apply(quantile_value_output_control,2,mean)
quantile_value_output_treatment_mean <- apply(quantile_value_output_treatment,2,mean)
expenditures_quantile_value_output_difference_mean <- apply(quantile_value_output_difference,2,mean)


 
expenditures_recovered_quantile_differences <- apply(quantile_value_output_difference,2,get_posterior_intervals_function)
expenditures_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)
expenditures_recovered_quantile_treatment <- apply(quantile_value_output_treatment,2,get_posterior_intervals_function)

## revenues! 

codafit_stan_draws <- codafit_stan_draws_revenues
S <- dim(codafit_stan_draws)[1]
M <- 2 # mixture components
component_probs <- matrix(NA,nrow=S,ncol=M)
quantile_vec <- seq(0.05, 0.95,0.1)
N <- length(quantile_vec)
quantile_value_output_control <- matrix(NA, nrow=S, ncol=N)
quantile_value_output_treatment <- matrix(NA, nrow=S, ncol=N)
quantile_value_output_difference <- matrix(NA, nrow=S, ncol=N)


for(s in 1:S){
  
  control_logmean <- (codafit_stan_draws[s,"mu[1]"])
  treatment_logmean <- (codafit_stan_draws[s,"mu[1]"] + codafit_stan_draws[s,"tau[1]"])
  control_logsd <- exp(codafit_stan_draws[s,"sigma_control[1]"])
  treatment_logsd <- exp(codafit_stan_draws[s,"sigma_control[1]"] + codafit_stan_draws[s,"sigma_TE[1]"])
  control_logit_params <- codafit_stan_draws[s,c("beta_full[1,1]", "beta_full[2,1]")]
  control_probs <- convert_logit_2_prob(control_logit_params)
  treatment_logit_params <- codafit_stan_draws[s,c("beta_full[1,2]", "beta_full[2,2]")]
  treatment_probs <- convert_logit_2_prob(control_logit_params + treatment_logit_params)
  

  
  for(n in 1:N){
    
    u <- quantile_vec[n]
    if(u < control_probs[1]){ quantile_value_output_control[s,n] <-  0 }
    if(u > control_probs[1]){ quantile_value_output_control[s,n] <-  qlnorm(((u-control_probs[1])/control_probs[2]), meanlog = control_logmean, sdlog = control_logsd, lower.tail = TRUE, log.p = FALSE) }
    if(u < treatment_probs[1]){ quantile_value_output_treatment[s,n] <-  0 }
    if(u > treatment_probs[1]){ quantile_value_output_treatment[s,n] <-  qlnorm(((u-treatment_probs[1])/treatment_probs[2]), meanlog = treatment_logmean, sdlog = treatment_logsd, lower.tail = TRUE, log.p = FALSE)  }
    
  } # close forloop indexed by N 
  
}
quantile_value_output_difference <-  quantile_value_output_treatment -  quantile_value_output_control


quantile_value_output_control_mean <- apply(quantile_value_output_control,2,mean)
quantile_value_output_treatment_mean <- apply(quantile_value_output_treatment,2,mean)
revenues_quantile_value_output_difference_mean <- apply(quantile_value_output_difference,2,mean)

revenues_recovered_quantile_differences <- apply(quantile_value_output_difference,2,get_posterior_intervals_function)
revenues_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)
revenues_recovered_quantile_treatment <- apply(quantile_value_output_treatment,2,get_posterior_intervals_function)


### FOR PROFIT

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

  control_logmean_1 <- (codafit_stan_draws[s,"mu[1]"])
  treatment_logmean_1 <- (codafit_stan_draws[s,"mu[1]"] + codafit_stan_draws[s,"tau[1]"])
  control_logsd_1 <- exp(codafit_stan_draws[s,"sigma_control[1]"])
  treatment_logsd_1 <- exp(codafit_stan_draws[s,"sigma_control[1]"] + codafit_stan_draws[s,"sigma_TE[1]"])
  
  control_logmean_2 <- (codafit_stan_draws[s,"mu[2]"])
  treatment_logmean_2 <- (codafit_stan_draws[s,"mu[2]"] + codafit_stan_draws[s,"tau[2]"])
  control_logsd_2 <- exp(codafit_stan_draws[s,"sigma_control[2]"])
  treatment_logsd_2 <- exp(codafit_stan_draws[s,"sigma_control[2]"] + codafit_stan_draws[s,"sigma_TE[2]"])
  
  control_logit_params <- codafit_stan_draws[s,c("beta_full[1,1]", "beta_full[2,1]", "beta_full[3,1]")]
  control_probs <- convert_logit_2_prob(control_logit_params)
  treatment_logit_params <- codafit_stan_draws[s,c("beta_full[1,2]", "beta_full[2,2]", "beta_full[3,2]")]
  treatment_probs <- convert_logit_2_prob(control_logit_params + treatment_logit_params)
  
  
  for(n in 1:N){
    u <- quantile_vec[n]
    if(u < control_probs[1]){ quantile_value_output_control[s,n]  <-  - qlnorm(1-(u/control_probs[1]), meanlog = control_logmean_1, sdlog = control_logsd_1, lower.tail = TRUE, log.p = FALSE) }
    if(control_probs[1] < u & u < (control_probs[1]+control_probs[2])){ quantile_value_output_control[s,n]  <-  0 }
    if(u > (control_probs[1]+control_probs[2])){ quantile_value_output_control[s,n]  <-  qlnorm(((u-(control_probs[1]+control_probs[2]))/(1-control_probs[1]-control_probs[2])),  meanlog = control_logmean_2, sdlog = control_logsd_2, lower.tail = TRUE, log.p = FALSE)  }
    
    if(u < treatment_probs[1]){ quantile_value_output_treatment[s,n]  <-  -qlnorm(1-(u/treatment_probs[1]), meanlog = treatment_logmean_1, sdlog = treatment_logsd_1, lower.tail = TRUE, log.p = FALSE) }
    if(treatment_probs[1] < u & u < (treatment_probs[1]+treatment_probs[2])){ quantile_value_output_treatment[s,n]  <-  0 }
    if(u > (treatment_probs[1]+treatment_probs[2])){ quantile_value_output_treatment[s,n]  <-  qlnorm(((u-(treatment_probs[1]+treatment_probs[2]))/(1-treatment_probs[1]-treatment_probs[2])), meanlog = treatment_logmean_2, sdlog = treatment_logsd_2, lower.tail = TRUE, log.p = FALSE)  }
    
    
  } # close forloop indexed by n 
  
} # closes forloop indexed by s 


quantile_value_output_difference <-  quantile_value_output_treatment -  quantile_value_output_control


quantile_value_output_control_mean <- apply(quantile_value_output_control,2,mean)
quantile_value_output_treatment_mean <- apply(quantile_value_output_treatment,2,mean)
profit_quantile_value_output_difference_mean <- apply(quantile_value_output_difference,2,mean)


profit_recovered_quantile_differences <- apply(quantile_value_output_difference,2,get_posterior_intervals_function)
profit_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)
profit_recovered_quantile_treatment <- apply(quantile_value_output_treatment,2,get_posterior_intervals_function)



### NOW THE GRAPHICS  (RESULTS FOR FIGURE 1 OF AER VERSION) ###

# PROFIT # 
# make ribbon plot for posterior 
quantiles_list <- quantile_vec
posterior_beta_data <- t(rbind(profit_recovered_quantile_differences, profit_quantile_value_output_difference_mean  ))
posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")

# tabular! 

median_beta_1_profit_partial_pooling <- round(t(profit_recovered_quantile_differences[3,]),1)
lower_ci_beta_1_profit_partial_pooling <- round(t(profit_recovered_quantile_differences[1,]),1)
upper_ci_beta_1_profit_partial_pooling <- round(t(profit_recovered_quantile_differences[5,]),1) 
ci_beta_1_profit_partial_pooling <- paste0("(", lower_ci_beta_1_profit_partial_pooling[,], ",", upper_ci_beta_1_profit_partial_pooling[,], ")" )
ci_beta_1_profit_partial_pooling <- matrix(ci_beta_1_profit_partial_pooling, nrow=1, ncol=10)
ci_beta_1_profit_partial_pooling <- data.frame(ci_beta_1_profit_partial_pooling )
median_beta_1_profit_partial_pooling <- data.frame(median_beta_1_profit_partial_pooling)
partial_pooling_beta_1_results_profit <- rbindlist(list(median_beta_1_profit_partial_pooling, ci_beta_1_profit_partial_pooling))
colnames(partial_pooling_beta_1_results_profit) <- c( "5th", "15th", "25th", "35th", "45th", "55th", "65th", "75th", "85th", "95th")

saveRDS(partial_pooling_beta_1_results_profit, "output/partial_pooling_beta_1_table_profit.rds")

fig_scale <- 0.4 # because the AER template will not let you scale it in tex 

posterior_beta_data_plot <- ggplot(posterior_beta_data, aes(posterior_beta_data$quantiles_list))
pdf("output/posterior_parent_quantile_TEs_profit_lognormal.pdf", width=fig_scale*6.5, height=fig_scale*6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "red", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "red", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"50%"]), color = "dark red", size = 1.5*fig_scale) +
  ggtitle("Posterior quantile effects on profit") +
  theme(plot.title = element_text(size = fig_scale*16)) + xlim(0.05,0.95) + ylim(-100,400) +
  xlab("Quantiles") + ylab("Quantile treatment effect")+
  theme(axis.text = element_text(size=fig_scale*14)) +  theme(axis.title.y = element_text(size = fig_scale*14)) +  theme(axis.title.x = element_text(size = fig_scale*14))
dev.off()

# expenditures # 
# make ribbon plot for posterior 
quantiles_list <- quantile_vec
posterior_beta_data <- t(rbind(expenditures_recovered_quantile_differences, expenditures_quantile_value_output_difference_mean  ))
posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")


posterior_beta_data_plot <- ggplot(posterior_beta_data, aes(posterior_beta_data$quantiles_list))
pdf("output/posterior_parent_quantile_TEs_expenditures_lognormal.pdf", width=fig_scale*6.5, height=fig_scale*6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "red", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "red", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"50%"]), color = "dark red", size = fig_scale*1.5) +
  ggtitle("Posterior quantile effects on expenditures") +
  theme(plot.title = element_text(size = fig_scale*16)) + xlim(0.05,0.95) + ylim(-100,400) +
  xlab("Quantiles") + ylab("Quantile treatment effect")+
  theme(axis.text = element_text(size=fig_scale*14)) +  theme(axis.title.y = element_text(size = fig_scale*14)) +  theme(axis.title.x = element_text(size = fig_scale*14))
dev.off()


# revenues # 
# make ribbon plot for posterior 
quantiles_list <- quantile_vec
posterior_beta_data <- t(rbind(revenues_recovered_quantile_differences, revenues_quantile_value_output_difference_mean  ))
posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")


posterior_beta_data_plot <- ggplot(posterior_beta_data, aes(posterior_beta_data$quantiles_list))
pdf("output/posterior_parent_quantile_TEs_revenues_lognormal.pdf", width=fig_scale*6.5, height=fig_scale*6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "red", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "red", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"50%"]), color = "dark red", size = fig_scale*1.5) +
  ggtitle("Posterior quantile effects on revenues") +
  theme(plot.title = element_text(size = fig_scale*16)) + xlim(0.05,0.95) + ylim(-100,400) +
  xlab("Quantiles") + ylab("Quantile treatment effect")+
  theme(axis.text = element_text(size=fig_scale*14)) +  theme(axis.title.y = element_text(size = fig_scale*14)) +  theme(axis.title.x = element_text(size = fig_scale*14))
dev.off()


### NOW THE FULL POOLING RESULTS FOR PROFIT 


codafit_stan_draws <- codafit_stan_draws_profit_full_pooling
S <- dim(codafit_stan_draws)[1]
M <- 3 # mixture components
component_probs <- matrix(NA,nrow=S,ncol=M)
quantile_vec <- seq(0.05, 0.95,0.1)
N <- length(quantile_vec)
quantile_value_output_control <- matrix(NA, nrow=S, ncol=N)
quantile_value_output_treatment <- matrix(NA, nrow=S, ncol=N)
quantile_value_output_difference <- matrix(NA, nrow=S, ncol=N)



for(s in 1:S){
  
  control_logmean_1 <- (codafit_stan_draws[s,"mu[1]"])
  treatment_logmean_1 <- (codafit_stan_draws[s,"mu[1]"] + codafit_stan_draws[s,"tau[1]"])
  control_logsd_1 <- exp(codafit_stan_draws[s,"sigma_control[1]"])
  treatment_logsd_1 <- exp(codafit_stan_draws[s,"sigma_control[1]"] + codafit_stan_draws[s,"sigma_TE[1]"])
  
  control_logmean_2 <- (codafit_stan_draws[s,"mu[2]"])
  treatment_logmean_2 <- (codafit_stan_draws[s,"mu[2]"] + codafit_stan_draws[s,"tau[2]"])
  control_logsd_2 <- exp(codafit_stan_draws[s,"sigma_control[2]"])
  treatment_logsd_2 <- exp(codafit_stan_draws[s,"sigma_control[2]"] + codafit_stan_draws[s,"sigma_TE[2]"])
  
  control_logit_params <- codafit_stan_draws[s,c("beta_full[1,1]", "beta_full[2,1]", "beta_full[3,1]")]
  control_probs <- convert_logit_2_prob(control_logit_params)
  treatment_logit_params <- codafit_stan_draws[s,c("beta_full[1,2]", "beta_full[2,2]", "beta_full[3,2]")]
  treatment_probs <- convert_logit_2_prob(control_logit_params + treatment_logit_params)
  
  
  for(n in 1:N){
    u <- quantile_vec[n]
    if(u < control_probs[1]){ quantile_value_output_control[s,n]  <-  - qlnorm(1-(u/control_probs[1]), meanlog = control_logmean_1, sdlog = control_logsd_1, lower.tail = TRUE, log.p = FALSE) }
    if(control_probs[1] < u & u < (control_probs[1]+control_probs[2])){ quantile_value_output_control[s,n]  <-  0 }
    if(u > (control_probs[1]+control_probs[2])){ quantile_value_output_control[s,n]  <-  qlnorm(((u-(control_probs[1]+control_probs[2]))/(1-control_probs[1]-control_probs[2])),  meanlog = control_logmean_2, sdlog = control_logsd_2, lower.tail = TRUE, log.p = FALSE)  }
    
    if(u < treatment_probs[1]){ quantile_value_output_treatment[s,n]  <-  -qlnorm(1-(u/treatment_probs[1]), meanlog = treatment_logmean_1, sdlog = treatment_logsd_1, lower.tail = TRUE, log.p = FALSE) }
    if(treatment_probs[1] < u & u < (treatment_probs[1]+treatment_probs[2])){ quantile_value_output_treatment[s,n]  <-  0 }
    if(u > (treatment_probs[1]+treatment_probs[2])){ quantile_value_output_treatment[s,n]  <-  qlnorm(((u-(treatment_probs[1]+treatment_probs[2]))/(1-treatment_probs[1]-treatment_probs[2])), meanlog = treatment_logmean_2, sdlog = treatment_logsd_2, lower.tail = TRUE, log.p = FALSE)  }
    
    
  } # close forloop indexed by n 
  
} # closes forloop indexed by s 


quantile_value_output_difference <-  quantile_value_output_treatment -  quantile_value_output_control


quantile_value_output_control_mean <- apply(quantile_value_output_control,2,mean)
quantile_value_output_treatment_mean <- apply(quantile_value_output_treatment,2,mean)
profit_quantile_value_output_difference_mean <- apply(quantile_value_output_difference,2,mean)


profit_recovered_quantile_differences <- apply(quantile_value_output_difference,2,get_posterior_intervals_function)
profit_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)
profit_recovered_quantile_treatment <- apply(quantile_value_output_treatment,2,get_posterior_intervals_function)



### NOW THE TABLE INPUT ###

# PROFIT # 
# make ribbon plot for posterior 
quantiles_list <- quantile_vec
posterior_beta_data <- t(rbind(profit_recovered_quantile_differences, profit_quantile_value_output_difference_mean  ))
posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")

country <- c("Mexico", "Mongolia", "Bosnia", "India", "Morocco", "Philippines", "Ethiopia")
median_beta_1_profit_full_pooling <- round(t(profit_recovered_quantile_differences[3,]),1)
lower_ci_beta_1_profit_full_pooling <- round(t(profit_recovered_quantile_differences[1,]),1)
upper_ci_beta_1_profit_full_pooling <- round(t(profit_recovered_quantile_differences[5,]),1) 
ci_beta_1_profit_full_pooling <- paste0("(", lower_ci_beta_1_profit_full_pooling[,], ",", upper_ci_beta_1_profit_full_pooling[,], ")" )
ci_beta_1_profit_full_pooling <- matrix(ci_beta_1_profit_full_pooling, nrow=1, ncol=10)
ci_beta_1_profit_full_pooling <- data.frame(ci_beta_1_profit_full_pooling )
median_beta_1_profit_full_pooling <- data.frame(median_beta_1_profit_full_pooling)
full_pooling_beta_1_results_profit <- rbindlist(list(median_beta_1_profit_full_pooling, ci_beta_1_profit_full_pooling))
colnames(full_pooling_beta_1_results_profit) <- c( "5th", "15th", "25th", "35th", "45th", "55th", "65th", "75th", "85th", "95th")

saveRDS(full_pooling_beta_1_results_profit, "output/full_pooling_beta_1_table_profit.rds")


##### POSTERIOR PREDICTIVE MODEL QUANTILE COMPUTATION 

# DRAW POSTERIOR PREDICTED PARAMETERS 

posterior_predicted_control_logmean_draws_profit_negative <- rnorm(length(codafit_stan_draws_profit[,"mu[1]"]),codafit_stan_draws_profit[,"mu[1]"],codafit_stan_draws_profit[,"sd_mu[1]"]  )
posterior_predicted_control_logmean_draws_profit_positive <- rnorm(length(codafit_stan_draws_profit[,"mu[2]"]),codafit_stan_draws_profit[,"mu[2]"],codafit_stan_draws_profit[,"sd_mu[2]"]  )
posterior_predicted_control_logmean_draws_expenditures <- rnorm(length(codafit_stan_draws_expenditures[,"mu[1]"]),codafit_stan_draws_expenditures[,"mu[1]"],codafit_stan_draws_expenditures[,"sd_mu[1]"]  )
posterior_predicted_control_logmean_draws_revenues <- rnorm(length(codafit_stan_draws_revenues[,"mu[1]"]),codafit_stan_draws_revenues[,"mu[1]"],codafit_stan_draws_revenues[,"sd_mu[1]"]  )

posterior_predicted_tau_logmean_draws_profit_negative <- rnorm(length(codafit_stan_draws_profit[,"tau[1]"]),codafit_stan_draws_profit[,"tau[1]"],codafit_stan_draws_profit[,"sd_tau[1]"]  )
posterior_predicted_tau_logmean_draws_profit_positive <- rnorm(length(codafit_stan_draws_profit[,"tau[2]"]),codafit_stan_draws_profit[,"tau[2]"],codafit_stan_draws_profit[,"sd_tau[2]"]  )
posterior_predicted_tau_logmean_draws_expenditures <- rnorm(length(codafit_stan_draws_expenditures[,"tau[1]"]),codafit_stan_draws_expenditures[,"tau[1]"],codafit_stan_draws_expenditures[,"sd_tau[1]"]  )
posterior_predicted_tau_logmean_draws_revenues <- rnorm(length(codafit_stan_draws_revenues[,"tau[1]"]),codafit_stan_draws_revenues[,"tau[1]"],codafit_stan_draws_revenues[,"sd_tau[1]"]  )


posterior_predicted_control_sigma_draws_profit_negative <- rnorm(length(codafit_stan_draws_profit[,"sigma_control[1]"]),codafit_stan_draws_profit[,"sigma_control[1]"],codafit_stan_draws_profit[,"sd_sigma_control[1]"]  )
posterior_predicted_control_sigma_draws_profit_positive <- rnorm(length(codafit_stan_draws_profit[,"sigma_control[2]"]),codafit_stan_draws_profit[,"sigma_control[2]"],codafit_stan_draws_profit[,"sd_sigma_control[2]"]  )
posterior_predicted_control_sigma_draws_expenditures <- rnorm(length(codafit_stan_draws_expenditures[,"sigma_control[1]"]),codafit_stan_draws_expenditures[,"sigma_control[1]"],codafit_stan_draws_expenditures[,"sd_sigma_control[1]"]  )
posterior_predicted_control_sigma_draws_revenues <- rnorm(length(codafit_stan_draws_revenues[,"sigma_control[1]"]),codafit_stan_draws_revenues[,"sigma_control[1]"],codafit_stan_draws_revenues[,"sd_sigma_control[1]"]  )

posterior_predicted_sigma_TE_draws_profit_negative <- rnorm(length(codafit_stan_draws_profit[,"sigma_TE[1]"]),codafit_stan_draws_profit[,"sigma_TE[1]"],codafit_stan_draws_profit[,"sd_sigma_TE[1]"]  )
posterior_predicted_sigma_TE_draws_profit_positive <- rnorm(length(codafit_stan_draws_profit[,"sigma_TE[2]"]),codafit_stan_draws_profit[,"sigma_TE[2]"],codafit_stan_draws_profit[,"sd_sigma_TE[2]"]  )
posterior_predicted_sigma_TE_draws_expenditures <- rnorm(length(codafit_stan_draws_expenditures[,"sigma_TE[1]"]),codafit_stan_draws_expenditures[,"sigma_TE[1]"],codafit_stan_draws_expenditures[,"sd_sigma_TE[1]"]  )
posterior_predicted_sigma_TE_draws_revenues <- rnorm(length(codafit_stan_draws_revenues[,"sigma_TE[1]"]),codafit_stan_draws_revenues[,"sigma_TE[1]"],codafit_stan_draws_revenues[,"sd_sigma_TE[1]"]  )

posterior_predicted_logit_control_draws_profit_negative <- rnorm(length(codafit_stan_draws_profit[,"beta_full[1,1]"]),codafit_stan_draws_profit[,"beta_full[1,1]"],codafit_stan_draws_profit[,"sigma[1,1]"]  )
posterior_predicted_logit_control_draws_profit_positive <- rnorm(length(codafit_stan_draws_profit[,"beta_full[2,1]"]),codafit_stan_draws_profit[,"beta_full[2,1]"],codafit_stan_draws_profit[,"sigma[2,1]"]  )
posterior_predicted_logit_control_draws_expenditures <- rnorm(length(codafit_stan_draws_expenditures[,"beta_full[1,1]"]),codafit_stan_draws_expenditures[,"beta_full[1,1]"],codafit_stan_draws_expenditures[,"sigma[1,1]"]  )
posterior_predicted_logit_control_draws_revenues <- rnorm(length(codafit_stan_draws_revenues[,"beta_full[1,1]"]),codafit_stan_draws_revenues[,"beta_full[1,1]"],codafit_stan_draws_revenues[,"sigma[1,1]"]  )


posterior_predicted_logit_te_draws_profit_negative <- rnorm(length(codafit_stan_draws_profit[,"beta_full[1,2]"]),codafit_stan_draws_profit[,"beta_full[1,2]"],codafit_stan_draws_profit[,"sigma[1,2]"]  )
posterior_predicted_logit_te_draws_profit_positive <- rnorm(length(codafit_stan_draws_profit[,"beta_full[2,2]"]),codafit_stan_draws_profit[,"beta_full[2,2]"],codafit_stan_draws_profit[,"sigma[2,2]"]  )
posterior_predicted_logit_te_draws_expenditures <- rnorm(length(codafit_stan_draws_expenditures[,"beta_full[1,2]"]),codafit_stan_draws_expenditures[,"beta_full[1,2]"],codafit_stan_draws_expenditures[,"sigma[1,2]"]  )
posterior_predicted_logit_te_draws_revenues <- rnorm(length(codafit_stan_draws_revenues[,"beta_full[1,2]"]),codafit_stan_draws_revenues[,"beta_full[1,2]"],codafit_stan_draws_revenues[,"sigma[1,2]"]  )

# NOW THE POSTERIOR PREDICTED QUANTILES

## FOR THE 1 TAIL OUTCOMES

## revenues! 

codafit_stan_draws <- codafit_stan_draws_revenues
S <- dim(codafit_stan_draws)[1]
M <- 2 # mixture components
component_probs <- matrix(NA,nrow=S,ncol=M)
quantile_vec <- seq(0.05, 0.95,0.1)
N <- length(quantile_vec)
quantile_value_output_control <- matrix(NA, nrow=S, ncol=N)
quantile_value_output_treatment <- matrix(NA, nrow=S, ncol=N)
quantile_value_output_difference <- matrix(NA, nrow=S, ncol=N)


for(s in 1:S){
  
  control_logmean <- posterior_predicted_control_logmean_draws_revenues[s]
  treatment_logmean <- posterior_predicted_control_logmean_draws_revenues[s] + posterior_predicted_tau_logmean_draws_revenues[s]
  control_logsd <- exp(posterior_predicted_control_sigma_draws_revenues[s])
  treatment_logsd <- exp(posterior_predicted_control_sigma_draws_revenues[s] + posterior_predicted_sigma_TE_draws_revenues[s])
  
  control_logit_params <- c(posterior_predicted_logit_control_draws_revenues[s],0)
  control_probs <- convert_logit_2_prob(control_logit_params)
  treatment_logit_params <- c(posterior_predicted_logit_te_draws_revenues[s],0)
  treatment_probs <- convert_logit_2_prob(control_logit_params + treatment_logit_params)
  
  
  for(n in 1:N){
    
    u <- quantile_vec[n]
    if(u < control_probs[1]){ quantile_value_output_control[s,n] <-  0 }
    if(u > control_probs[1]){ quantile_value_output_control[s,n] <-  qlnorm(((u-control_probs[1])/control_probs[2]), meanlog = control_logmean, sdlog = control_logsd, lower.tail = TRUE, log.p = FALSE)}
    if(u < treatment_probs[1]){ quantile_value_output_treatment[s,n] <-  0 }
    if(u > treatment_probs[1]){ quantile_value_output_treatment[s,n] <-  qlnorm(((u-treatment_probs[1])/treatment_probs[2]), meanlog = treatment_logmean, sdlog = treatment_logsd, lower.tail = TRUE, log.p = FALSE) }
    
  } # close forloop indexed by N 
  
}
quantile_value_output_difference <-  quantile_value_output_treatment -  quantile_value_output_control


quantile_value_output_control_mean <- apply(quantile_value_output_control,2,mean)
quantile_value_output_treatment_mean <- apply(quantile_value_output_treatment,2,mean)
revenues_quantile_value_output_difference_mean <- apply(quantile_value_output_difference,2,mean)

revenues_recovered_quantile_differences <- apply(quantile_value_output_difference,2,get_posterior_intervals_function)
revenues_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)
revenues_recovered_quantile_treatment <- apply(quantile_value_output_treatment,2,get_posterior_intervals_function)


## expenditures! 

codafit_stan_draws <- codafit_stan_draws_expenditures
S <- dim(codafit_stan_draws)[1]
M <- 2 # mixture components
component_probs <- matrix(NA,nrow=S,ncol=M)
quantile_vec <- seq(0.05, 0.95,0.1)
N <- length(quantile_vec)
quantile_value_output_control <- matrix(NA, nrow=S, ncol=N)
quantile_value_output_treatment <- matrix(NA, nrow=S, ncol=N)
quantile_value_output_difference <- matrix(NA, nrow=S, ncol=N)


for(s in 1:S){

  control_logmean <- posterior_predicted_control_logmean_draws_expenditures[s]
  treatment_logmean <- posterior_predicted_control_logmean_draws_expenditures[s] + posterior_predicted_tau_logmean_draws_expenditures[s]
  control_logsd <- exp(posterior_predicted_control_sigma_draws_expenditures[s])
  treatment_logsd <- exp(posterior_predicted_control_sigma_draws_expenditures[s] + posterior_predicted_sigma_TE_draws_expenditures[s])
  
  control_logit_params <- c(posterior_predicted_logit_control_draws_expenditures[s],0)
  control_probs <- convert_logit_2_prob(control_logit_params)
  treatment_logit_params <- c(posterior_predicted_logit_te_draws_expenditures[s],0)
  treatment_probs <- convert_logit_2_prob(control_logit_params + treatment_logit_params)
  
  
  for(n in 1:N){

    u <- quantile_vec[n]
    if(u < control_probs[1]){ quantile_value_output_control[s,n] <-  0 }
    if(u > control_probs[1]){ quantile_value_output_control[s,n] <-  qlnorm(((u-control_probs[1])/control_probs[2]), meanlog = control_logmean, sdlog = control_logsd, lower.tail = TRUE, log.p = FALSE)}
    if(u < treatment_probs[1]){ quantile_value_output_treatment[s,n] <-  0 }
    if(u > treatment_probs[1]){ quantile_value_output_treatment[s,n] <-  qlnorm(((u-treatment_probs[1])/treatment_probs[2]), meanlog = treatment_logmean, sdlog = treatment_logsd, lower.tail = TRUE, log.p = FALSE) }
    
  } # close forloop indexed by N 
  
}
quantile_value_output_difference <-  quantile_value_output_treatment -  quantile_value_output_control


quantile_value_output_control_mean <- apply(quantile_value_output_control,2,mean)
quantile_value_output_treatment_mean <- apply(quantile_value_output_treatment,2,mean)
expenditures_quantile_value_output_difference_mean <- apply(quantile_value_output_difference,2,mean)

expenditures_recovered_quantile_differences <- apply(quantile_value_output_difference,2,get_posterior_intervals_function)
expenditures_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)
expenditures_recovered_quantile_treatment <- apply(quantile_value_output_treatment,2,get_posterior_intervals_function)

# PROFIT

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
  

  control_logmean_1 <- posterior_predicted_control_logmean_draws_profit_negative[s]
  treatment_logmean_1 <- posterior_predicted_control_logmean_draws_profit_negative[s] + posterior_predicted_tau_logmean_draws_profit_negative[s]
  control_logsd_1 <- exp(posterior_predicted_control_sigma_draws_profit_negative[s])
  treatment_logsd_1 <- exp(posterior_predicted_control_sigma_draws_profit_negative[s] + posterior_predicted_sigma_TE_draws_profit_negative[s])
  
  control_logmean_2 <- posterior_predicted_control_logmean_draws_profit_positive[s]
  treatment_logmean_2 <- posterior_predicted_control_logmean_draws_profit_positive[s] + posterior_predicted_tau_logmean_draws_profit_positive[s]
  control_logsd_2 <- exp(posterior_predicted_control_sigma_draws_profit_positive[s])
  treatment_logsd_2 <- exp(posterior_predicted_control_sigma_draws_profit_positive[s] + posterior_predicted_sigma_TE_draws_profit_positive[s])
  
  
  control_logit_params <- c(posterior_predicted_logit_control_draws_profit_negative[s], posterior_predicted_logit_control_draws_profit_positive[s],0)
  control_probs <- convert_logit_2_prob(control_logit_params)
  treatment_logit_params <- c(posterior_predicted_logit_te_draws_profit_negative[s], posterior_predicted_logit_te_draws_profit_positive[s],0)
  treatment_probs <- convert_logit_2_prob(control_logit_params + treatment_logit_params)
  
  
  for(n in 1:N){
    
    u <- quantile_vec[n]
    if(u < control_probs[1]){ quantile_value_output_control[s,n]  <-  - qlnorm(1-(u/control_probs[1]), meanlog = control_logmean_1, sdlog = control_logsd_1, lower.tail = TRUE, log.p = FALSE) }
    if(control_probs[1] < u & u < (control_probs[1]+control_probs[2])){ quantile_value_output_control[s,n]  <-  0 }
    if(u > (control_probs[1]+control_probs[2])){ quantile_value_output_control[s,n]  <-  qlnorm(((u-(control_probs[1]+control_probs[2]))/(1-control_probs[1]-control_probs[2])),  meanlog = control_logmean_2, sdlog = control_logsd_2, lower.tail = TRUE, log.p = FALSE)  }
    
    if(u < treatment_probs[1]){ quantile_value_output_treatment[s,n]  <-  -qlnorm(1-(u/treatment_probs[1]), meanlog = treatment_logmean_1, sdlog = treatment_logsd_1, lower.tail = TRUE, log.p = FALSE) }
    if(treatment_probs[1] < u & u < (treatment_probs[1]+treatment_probs[2])){ quantile_value_output_treatment[s,n]  <-  0 }
    if(u > (treatment_probs[1]+treatment_probs[2])){ quantile_value_output_treatment[s,n]  <-  qlnorm(((u-(treatment_probs[1]+treatment_probs[2]))/(1-treatment_probs[1]-treatment_probs[2])), meanlog = treatment_logmean_2, sdlog = treatment_logsd_2, lower.tail = TRUE, log.p = FALSE)  }
    
    
  } # close forloop indexed by n 
  
} # closes forloop indexed by s 


quantile_value_output_difference <-  quantile_value_output_treatment -  quantile_value_output_control


quantile_value_output_control_mean <- apply(quantile_value_output_control,2,mean)
quantile_value_output_treatment_mean <- apply(quantile_value_output_treatment,2,mean)
profit_quantile_value_output_difference_mean <- apply(quantile_value_output_difference,2,mean)


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


posterior_beta_data_plot <- ggplot(posterior_beta_data, aes(posterior_beta_data$quantiles_list))
pdf("output/posterior_predicted_quantile_TEs_profit_lognormal.pdf", width=fig_scale*6.5, height=fig_scale*6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "red", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "red", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"50%"]), color = "dark red", size = fig_scale*1.5) +
  ggtitle("Posterior predicted quantile effects on profit") +
  theme(plot.title = element_text(size = 16)) + xlim(0.05,0.95) + ylim(-1500,6000) +
  xlab("Quantiles") + ylab("Quantile treatment effect")+
  theme(axis.text = element_text(size=fig_scale*14)) +  theme(axis.title.y = element_text(size = fig_scale*14)) +  theme(axis.title.x = element_text(size = fig_scale*14))
dev.off()

# expenditures # 
# make ribbon plot for posterior 
quantiles_list <- quantile_vec
posterior_beta_data <- t(rbind(expenditures_recovered_quantile_differences, expenditures_quantile_value_output_difference_mean  ))
posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")


posterior_beta_data_plot <- ggplot(posterior_beta_data, aes(posterior_beta_data$quantiles_list))
pdf("output/posterior_predicted_quantile_TEs_expenditures_lognormal.pdf", width=fig_scale*6.5, height=fig_scale*6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "red", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "red", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"50%"]), color = "dark red", size = fig_scale*1.5) +
  ggtitle("Posterior predicted quantile effects on expenditures") +
  theme(plot.title = element_text(size = fig_scale*14)) + xlim(0.05,0.95) + ylim(-2500,10600) +
  xlab("Quantiles") +ylab("Quantile treatment effect")+
  theme(axis.text = element_text(size=fig_scale*14)) +  theme(axis.title.y = element_text(size = fig_scale*14)) +  theme(axis.title.x = element_text(size = fig_scale*14))
dev.off()


# revenues # 
# make ribbon plot for posterior 
quantiles_list <- quantile_vec
posterior_beta_data <- t(rbind(revenues_recovered_quantile_differences, revenues_quantile_value_output_difference_mean  ))
posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")


posterior_beta_data_plot <- ggplot(posterior_beta_data, aes(posterior_beta_data$quantiles_list))
pdf("output/posterior_predicted_quantile_TEs_revenues_lognormal.pdf", width=fig_scale*6.5, height=fig_scale*6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "red", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "red", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"50%"]), color = "dark red", size = fig_scale*1.5) +
  ggtitle("Posterior predicted quantile effects on revenues") +
  theme(plot.title = element_text(size = fig_scale*14)) + xlim(0.05,0.95)  + ylim(-2525,10100) +
  xlab("Quantiles") + ylab("Quantile treatment effect")+
  theme(axis.text = element_text(size=fig_scale*14)) +  theme(axis.title.y = element_text(size = fig_scale*14)) +  theme(axis.title.x = element_text(size = fig_scale*14))
dev.off()

