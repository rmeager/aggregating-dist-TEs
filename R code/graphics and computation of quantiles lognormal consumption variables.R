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

# install and load packages


installation_needed  <- FALSE
loading_needed <- TRUE
package_list <- c('ggplot2', 'rstan','reshape','reshape2','coda','xtable', 'dplyr', 'Runuran', 'testthat',
                  "MCMCpack", "geoR", "gtools", 'gPdtest', 'fBasics',"PtProcess", "VGAM", "MASS","quantreg",
                  "boot")
if(installation_needed){install.packages(package_list, repos='http://cran.us.r-project.org')}
if(loading_needed){lapply(package_list, require, character.only = TRUE)}


# Load data
load("output/microcredit_consumerdurables_lognormal_tailored_hierarchical_pdf_output_5000_iters.RData")
stan_fit_table_consumerdurables <- stan_fit_table
stan_fit_consumerdurables <- stan_fit
codafit_stan_draws_consumerdurables <- as.matrix(stan2coda(stan_fit_consumerdurables))
data_split_consumerdurables <- data_split
data_consumerdurables <- data

load("output/microcredit_consumption_lognormal_tailored_hierarchical_pdf_output_5000_iters.RData")
stan_fit_table_consumption <- stan_fit_table
stan_fit_consumption <- stan_fit
codafit_stan_draws_consumption <- as.matrix(stan2coda(stan_fit_consumption))
data_split_consumption <- data_split
data_consumption <- data

load("output/microcredit_temptation_lognormal_tailored_hierarchical_pdf_output_5000_iters.RData")
stan_fit_table_temptation <- stan_fit_table
stan_fit_temptation <- stan_fit
codafit_stan_draws_temptation <- as.matrix(stan2coda(stan_fit_temptation))
data_split_temptation <- data_split
data_temptation <- data


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

fig_scale <- 0.4 # because the AER template will not let you scale it in tex 

## CONSUMPTION



codafit_stan_draws <- codafit_stan_draws_consumption
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
consumption_quantile_value_output_difference_mean <- apply(quantile_value_output_difference,2,mean)

consumption_recovered_quantile_differences <- apply(quantile_value_output_difference,2,get_posterior_intervals_function)
consumption_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)
consumption_recovered_quantile_treatment <- apply(quantile_value_output_treatment,2,get_posterior_intervals_function)




## CONSUMERDURABLES
codafit_stan_draws <- codafit_stan_draws_consumerdurables
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
consumerdurables_quantile_value_output_difference_mean <- apply(quantile_value_output_difference,2,mean)

consumerdurables_recovered_quantile_differences <- apply(quantile_value_output_difference,2,get_posterior_intervals_function)
consumerdurables_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)
consumerdurables_recovered_quantile_treatment <- apply(quantile_value_output_treatment,2,get_posterior_intervals_function)

## temptation
codafit_stan_draws <- codafit_stan_draws_temptation
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
temptation_quantile_value_output_difference_mean <- apply(quantile_value_output_difference,2,mean)

temptation_recovered_quantile_differences <- apply(quantile_value_output_difference,2,get_posterior_intervals_function)
temptation_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)
temptation_recovered_quantile_treatment <- apply(quantile_value_output_treatment,2,get_posterior_intervals_function)



# temptation # 
# make ribbon plot for posterior 
quantiles_list <- quantile_vec
posterior_beta_data <- t(rbind(temptation_recovered_quantile_differences, temptation_quantile_value_output_difference_mean  ))
posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")


posterior_beta_data_plot <- ggplot(posterior_beta_data, aes(posterior_beta_data$quantiles_list))
pdf("output/posterior_parent_quantile_TEs_temptation_lognormal.pdf", width=fig_scale*6.5, height=fig_scale*6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "red", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "red", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"mean"]), color = "dark red", size = fig_scale*1.5) +
  ggtitle("Posterior quantile effects on temptation") +
  theme(plot.title = element_text(size = fig_scale*16)) + xlim(0.05,0.95) + ylim(-15,60) +
  xlab("Quantiles") + ylab("Quantile treatment effect")+
  theme(axis.text = element_text(size=fig_scale*14)) +  theme(axis.title.y = element_text(size = fig_scale*14)) +  theme(axis.title.x = element_text(size = fig_scale*14))
dev.off()


# consumerdurables # 
# make ribbon plot for posterior 
quantiles_list <- quantile_vec
posterior_beta_data <- t(rbind(consumerdurables_recovered_quantile_differences, consumerdurables_quantile_value_output_difference_mean  ))
posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")


posterior_beta_data_plot <- ggplot(posterior_beta_data, aes(posterior_beta_data$quantiles_list))
pdf("output/posterior_parent_quantile_TEs_consumerdurables_lognormal.pdf", width=fig_scale*6.5, height=fig_scale*6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "red", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "red", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"50%"]), color = "dark red", size = fig_scale*1.5) +
  ggtitle("Posterior quantile effects on consumerdurables") +
  theme(plot.title = element_text(size = fig_scale*16)) + xlim(0.05,0.95) + ylim(-275,1100) +
  xlab("Quantiles") + ylab("Quantile treatment effect")+
  theme(axis.text = element_text(size=fig_scale*14)) +  theme(axis.title.y = element_text(size = fig_scale*14)) +  theme(axis.title.x = element_text(size = fig_scale*14))
dev.off()


# consumption # 
# make ribbon plot for posterior 
quantiles_list <- quantile_vec
posterior_beta_data <- t(rbind(consumption_recovered_quantile_differences, consumption_quantile_value_output_difference_mean  ))
posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")
posterior_beta_data

posterior_beta_data_plot <- ggplot(posterior_beta_data, aes(posterior_beta_data$quantiles_list))
pdf("output/posterior_parent_quantile_TEs_consumption_lognormal.pdf", width=fig_scale*6.5, height=fig_scale*6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "red", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "red", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"mean"]), color = "dark red", size = 1.5) +
  ggtitle("Posterior quantile effects on consumption") +
  theme(plot.title = element_text(size = fig_scale*16)) + xlim(0.05,0.95) + ylim(-100,400) +
  xlab("Quantiles") + ylab("Quantile treatment effect")+
  theme(axis.text = element_text(size=fig_scale*14)) +  theme(axis.title.y = element_text(size = fig_scale*14)) +  theme(axis.title.x = element_text(size = fig_scale*14))
dev.off()


# add the average, again this must be computed earlier
quantiles_list <- quantile_vec
posterior_beta_data <- t(rbind(consumption_recovered_quantile_differences, consumption_quantile_value_output_difference_mean  ))
posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")

# tabular!

median_beta_1_consumption_partial_pooling <- round(t(consumption_recovered_quantile_differences[3,]),1)
lower_ci_beta_1_consumption_partial_pooling <- round(t(consumption_recovered_quantile_differences[1,]),1)
upper_ci_beta_1_consumption_partial_pooling <- round(t(consumption_recovered_quantile_differences[5,]),1) 
ci_beta_1_consumption_partial_pooling <- paste0("(", lower_ci_beta_1_consumption_partial_pooling[,], ",", upper_ci_beta_1_consumption_partial_pooling[,], ")" )
ci_beta_1_consumption_partial_pooling <- matrix(ci_beta_1_consumption_partial_pooling, nrow=1, ncol=10)
ci_beta_1_consumption_partial_pooling <- data.frame(ci_beta_1_consumption_partial_pooling )
median_beta_1_consumption_partial_pooling <- data.frame(median_beta_1_consumption_partial_pooling)
partial_pooling_beta_1_results_consumption <- rbindlist(list(median_beta_1_consumption_partial_pooling, ci_beta_1_consumption_partial_pooling))
colnames(partial_pooling_beta_1_results_consumption) <- c( "5th", "15th", "25th", "35th", "45th", "55th", "65th", "75th", "85th", "95th")
saveRDS(partial_pooling_beta_1_results_consumption, "output/partial_pooling_consumption_avg.rds")


##### {POSTERIOR PREDICTIVE QUANTILE DIFFERENCES}

# DRAW POSTERIOR PREDICTED PARAMETERS 

posterior_predicted_control_logmean_draws_consumerdurables <- rnorm(length(codafit_stan_draws_consumerdurables[,"mu[1]"]),codafit_stan_draws_consumerdurables[,"mu[1]"],codafit_stan_draws_consumerdurables[,"sd_mu[1]"]  )
posterior_predicted_control_logmean_draws_consumption <- rnorm(length(codafit_stan_draws_consumption[,"mu[1]"]),codafit_stan_draws_consumption[,"mu[1]"],codafit_stan_draws_consumption[,"sd_mu[1]"]  )
posterior_predicted_control_logmean_draws_temptation <- rnorm(length(codafit_stan_draws_temptation[,"mu[1]"]),codafit_stan_draws_temptation[,"mu[1]"],codafit_stan_draws_temptation[,"sd_mu[1]"]  )

posterior_predicted_tau_logmean_draws_consumerdurables <- rnorm(length(codafit_stan_draws_consumerdurables[,"tau[1]"]),codafit_stan_draws_consumerdurables[,"tau[1]"],codafit_stan_draws_consumerdurables[,"sd_tau[1]"]  )
posterior_predicted_tau_logmean_draws_consumption <- rnorm(length(codafit_stan_draws_consumption[,"tau[1]"]),codafit_stan_draws_consumption[,"tau[1]"],codafit_stan_draws_consumption[,"sd_tau[1]"]  )
posterior_predicted_tau_logmean_draws_temptation <- rnorm(length(codafit_stan_draws_temptation[,"tau[1]"]),codafit_stan_draws_temptation[,"tau[1]"],codafit_stan_draws_temptation[,"sd_tau[1]"]  )


posterior_predicted_control_sigma_draws_consumerdurables <- rnorm(length(codafit_stan_draws_consumerdurables[,"sigma_control[1]"]),codafit_stan_draws_consumerdurables[,"sigma_control[1]"],codafit_stan_draws_consumerdurables[,"sd_sigma_control[1]"]  )
posterior_predicted_control_sigma_draws_consumption <- rnorm(length(codafit_stan_draws_consumption[,"sigma_control[1]"]),codafit_stan_draws_consumption[,"sigma_control[1]"],codafit_stan_draws_consumption[,"sd_sigma_control[1]"]  )
posterior_predicted_control_sigma_draws_temptation <- rnorm(length(codafit_stan_draws_temptation[,"sigma_control[1]"]),codafit_stan_draws_temptation[,"sigma_control[1]"],codafit_stan_draws_temptation[,"sd_sigma_control[1]"]  )

posterior_predicted_sigma_TE_draws_consumerdurables <- rnorm(length(codafit_stan_draws_consumerdurables[,"sigma_TE[1]"]),codafit_stan_draws_consumerdurables[,"sigma_TE[1]"],codafit_stan_draws_consumerdurables[,"sd_sigma_TE[1]"]  )
posterior_predicted_sigma_TE_draws_consumption <- rnorm(length(codafit_stan_draws_consumption[,"sigma_TE[1]"]),codafit_stan_draws_consumption[,"sigma_TE[1]"],codafit_stan_draws_consumption[,"sd_sigma_TE[1]"]  )
posterior_predicted_sigma_TE_draws_temptation <- rnorm(length(codafit_stan_draws_temptation[,"sigma_TE[1]"]),codafit_stan_draws_temptation[,"sigma_TE[1]"],codafit_stan_draws_temptation[,"sd_sigma_TE[1]"]  )

posterior_predicted_logit_control_draws_consumerdurables <- rnorm(length(codafit_stan_draws_consumerdurables[,"beta_full[1,1]"]),codafit_stan_draws_consumerdurables[,"beta_full[1,1]"],codafit_stan_draws_consumerdurables[,"sigma[1,1]"]  )
posterior_predicted_logit_control_draws_consumption <- rnorm(length(codafit_stan_draws_consumption[,"beta_full[1,1]"]),codafit_stan_draws_consumption[,"beta_full[1,1]"],codafit_stan_draws_consumption[,"sigma[1,1]"]  )
posterior_predicted_logit_control_draws_temptation <- rnorm(length(codafit_stan_draws_temptation[,"beta_full[1,1]"]),codafit_stan_draws_temptation[,"beta_full[1,1]"],codafit_stan_draws_temptation[,"sigma[1,1]"]  )


posterior_predicted_logit_te_draws_consumerdurables <- rnorm(length(codafit_stan_draws_consumerdurables[,"beta_full[1,2]"]),codafit_stan_draws_consumerdurables[,"beta_full[1,2]"],codafit_stan_draws_consumerdurables[,"sigma[1,2]"]  )
posterior_predicted_logit_te_draws_consumption <- rnorm(length(codafit_stan_draws_consumption[,"beta_full[1,2]"]),codafit_stan_draws_consumption[,"beta_full[1,2]"],codafit_stan_draws_consumption[,"sigma[1,2]"]  )
posterior_predicted_logit_te_draws_temptation <- rnorm(length(codafit_stan_draws_temptation[,"beta_full[1,2]"]),codafit_stan_draws_temptation[,"beta_full[1,2]"],codafit_stan_draws_temptation[,"sigma[1,2]"]  )

# NOW THE POSTERIOR PREDICTED QUANTILES

## FOR THE 1 TAIL OUTCOMES

## consumption! 

codafit_stan_draws <- codafit_stan_draws_consumption
S <- dim(codafit_stan_draws)[1]
M <- 2 # mixture components
component_probs <- matrix(NA,nrow=S,ncol=M)
quantile_vec <- seq(0.05, 0.95,0.1)
N <- length(quantile_vec)
quantile_value_output_control <- matrix(NA, nrow=S, ncol=N)
quantile_value_output_treatment <- matrix(NA, nrow=S, ncol=N)
quantile_value_output_difference <- matrix(NA, nrow=S, ncol=N)


for(s in 1:S){
  
  control_logmean <- posterior_predicted_control_logmean_draws_consumption[s]
  treatment_logmean <- posterior_predicted_control_logmean_draws_consumption[s] + posterior_predicted_tau_logmean_draws_consumption[s]
  control_logsd <- exp(posterior_predicted_control_sigma_draws_consumption[s])
  treatment_logsd <- exp(posterior_predicted_control_sigma_draws_consumption[s] + posterior_predicted_sigma_TE_draws_consumption[s])
  
  control_logit_params <- c(posterior_predicted_logit_control_draws_consumption[s],0)
  control_probs <- convert_logit_2_prob(control_logit_params)
  treatment_logit_params <- c(posterior_predicted_logit_te_draws_consumption[s],0)
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
quantile_value_output_difference <- quantile_value_output_difference[complete.cases(quantile_value_output_difference),]


quantile_value_output_control_mean <- apply(quantile_value_output_control,2,mean)
quantile_value_output_treatment_mean <- apply(quantile_value_output_treatment,2,mean)
consumption_quantile_value_output_difference_mean <- apply(quantile_value_output_difference,2,mean)

consumption_recovered_quantile_differences <- apply(quantile_value_output_difference,2,get_posterior_intervals_function)
consumption_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)
consumption_recovered_quantile_treatment <- apply(quantile_value_output_treatment,2,get_posterior_intervals_function)
consumption_recovered_quantile_treatment - consumption_recovered_quantile_control # this probably is NOT a valid propagation 

## consumerdurables! 

codafit_stan_draws <- codafit_stan_draws_consumerdurables
S <- dim(codafit_stan_draws)[1]
M <- 2 # mixture components
component_probs <- matrix(NA,nrow=S,ncol=M)
quantile_vec <- seq(0.05, 0.95,0.1)
N <- length(quantile_vec)
quantile_value_output_control <- matrix(NA, nrow=S, ncol=N)
quantile_value_output_treatment <- matrix(NA, nrow=S, ncol=N)
quantile_value_output_difference <- matrix(NA, nrow=S, ncol=N)


for(s in 1:S){

  control_logmean <- posterior_predicted_control_logmean_draws_consumerdurables[s]
  treatment_logmean <- posterior_predicted_control_logmean_draws_consumerdurables[s] + posterior_predicted_tau_logmean_draws_consumerdurables[s]
  control_logsd <- exp(posterior_predicted_control_sigma_draws_consumerdurables[s])
  treatment_logsd <- exp(posterior_predicted_control_sigma_draws_consumerdurables[s] + posterior_predicted_sigma_TE_draws_consumerdurables[s])
  
  control_logit_params <- c(posterior_predicted_logit_control_draws_consumerdurables[s],0)
  control_probs <- convert_logit_2_prob(control_logit_params)
  treatment_logit_params <- c(posterior_predicted_logit_te_draws_consumerdurables[s],0)
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
consumerdurables_quantile_value_output_difference_mean <- apply(quantile_value_output_difference,2,mean)
quantile_value_output_difference <- quantile_value_output_difference[complete.cases(quantile_value_output_difference),]


consumerdurables_recovered_quantile_differences <- apply(quantile_value_output_difference,2,get_posterior_intervals_function)
consumerdurables_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)
consumerdurables_recovered_quantile_treatment <- apply(quantile_value_output_treatment,2,get_posterior_intervals_function)

## temptation! 

codafit_stan_draws <- codafit_stan_draws_temptation
S <- dim(codafit_stan_draws)[1]
M <- 2 # mixture components
component_probs <- matrix(NA,nrow=S,ncol=M)
quantile_vec <- seq(0.05, 0.95,0.1)
N <- length(quantile_vec)
quantile_value_output_control <- matrix(NA, nrow=S, ncol=N)
quantile_value_output_treatment <- matrix(NA, nrow=S, ncol=N)
quantile_value_output_difference <- matrix(NA, nrow=S, ncol=N)


for(s in 1:S){
  
  control_logmean <- posterior_predicted_control_logmean_draws_temptation[s]
  treatment_logmean <- posterior_predicted_control_logmean_draws_temptation[s] + posterior_predicted_tau_logmean_draws_temptation[s]
  control_logsd <- exp(posterior_predicted_control_sigma_draws_temptation[s])
  treatment_logsd <- exp(posterior_predicted_control_sigma_draws_temptation[s] + posterior_predicted_sigma_TE_draws_temptation[s])
  
  control_logit_params <- c(posterior_predicted_logit_control_draws_temptation[s],0)
  control_probs <- convert_logit_2_prob(control_logit_params)
  treatment_logit_params <- c(posterior_predicted_logit_te_draws_temptation[s],0)
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
temptation_quantile_value_output_difference_mean <- apply(quantile_value_output_difference,2,mean)

temptation_recovered_quantile_differences <- apply(quantile_value_output_difference,2,get_posterior_intervals_function)
temptation_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)
temptation_recovered_quantile_treatment <- apply(quantile_value_output_treatment,2,get_posterior_intervals_function)




### NOW THE GRAPHICS ###

# temptation # 
# make ribbon plot for posterior 
quantiles_list <- quantile_vec
posterior_beta_data <- t(rbind(temptation_recovered_quantile_differences, temptation_quantile_value_output_difference_mean  ))
posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")


posterior_beta_data_plot <- ggplot(posterior_beta_data, aes(posterior_beta_data$quantiles_list))
pdf("output/posterior_predicted_quantile_TEs_temptation_lognormal.pdf", width=fig_scale*6.5, height=fig_scale*6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "red", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "red", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"50%"]), color = "dark red", size = fig_scale*1.5) +
  ggtitle("Posterior predicted quantile effects on temptation") +
  theme(plot.title = element_text(size = fig_scale*14)) + xlim(0.05,0.95) + ylim(-100,400) +
  xlab("Quantiles") +ylab("Quantile treatment effect")+
  theme(axis.text = element_text(size=fig_scale*14)) +  theme(axis.title.y = element_text(size = fig_scale*14)) +  theme(axis.title.x = element_text(size = fig_scale*14))
dev.off()


# consumerdurables # 
# make ribbon plot for posterior 
quantiles_list <- quantile_vec
posterior_beta_data <- t(rbind(consumerdurables_recovered_quantile_differences, consumerdurables_quantile_value_output_difference_mean  ))
posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")


posterior_beta_data_plot <- ggplot(posterior_beta_data, aes(posterior_beta_data$quantiles_list))
pdf("output/posterior_predicted_quantile_TEs_consumerdurables_lognormal.pdf", width=fig_scale*6.5, height=fig_scale*6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "red", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "red", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"50%"]), color = "dark red", size = fig_scale*1.5) +
  ggtitle("Posterior predicted quantile effects on consumerdurables") +
  theme(plot.title = element_text(size = fig_scale*14)) + xlim(0.05,0.95) +  ylim(-60000,240000) +
  xlab("Quantiles") +ylab("Quantile treatment effect")+
  theme(axis.text = element_text(size=fig_scale*14)) +  theme(axis.title.y = element_text(size = fig_scale*14)) +  theme(axis.title.x = element_text(size = fig_scale*14))
dev.off()


# consumption # 
# make ribbon plot for posterior 
quantiles_list <- quantile_vec
posterior_beta_data <- t(rbind(consumption_recovered_quantile_differences, consumption_quantile_value_output_difference_mean  ))
posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")


posterior_beta_data_plot <- ggplot(posterior_beta_data, aes(posterior_beta_data$quantiles_list))
pdf("output/posterior_predicted_quantile_TEs_consumption_lognormal.pdf", width=fig_scale*6.5, height=fig_scale*6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "red", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "red", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"50%"]), color = "dark red", size = fig_scale*1.5) +
  ggtitle("Posterior predicted quantile effects on consumption") +
  theme(plot.title = element_text(size = fig_scale*14)) + xlim(0.05,0.95)  + ylim(-1500,6000) +
  xlab("Quantiles") + ylab("Quantile treatment effect")+
  theme(axis.text = element_text(size=fig_scale*14)) +  theme(axis.title.y = element_text(size = fig_scale*14)) +  theme(axis.title.x = element_text(size = fig_scale*14))
dev.off()



