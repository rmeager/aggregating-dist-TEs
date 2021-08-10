# graphics and quantile computation for the composite tail
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
load("output/microcredit_profit_composite_tails_tailored_hierarchical_pdf_output.RData")

stan_fit_table_profit <- stan_fit_table
stan_fit_profit <- stan_fit
codafit_stan_draws_profit <- as.matrix(stan2coda(stan_fit_profit))
data_split_profit <- data_split
data_profit <- data

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
  posterior_interval <- quantile(vector_draws, c(0.025, 0.25, 0.5, 0.75,0.975), na.rm = TRUE)
  return(posterior_interval)
  print(posterior_interval)
}

library(extraDistr)

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
  
  control_pareto_scale_1 <- (codafit_stan_draws[s,"pareto_control_shape[1]"])
  control_pareto_scale_2 <- (codafit_stan_draws[s,"pareto_control_shape[2]"])
  treatment_pareto_scale_1 <- exp( (codafit_stan_draws[s,"pareto_control_shape[1]"]) + (codafit_stan_draws[s,"pareto_tau_shape[1]"])       )
  treatment_pareto_scale_2 <- exp(  (codafit_stan_draws[s,"pareto_control_shape[2]"]) + (codafit_stan_draws[s,"pareto_tau_shape[2]"])          )
  
  control_logit_params <- codafit_stan_draws[s,c("beta_full[1,1]", "beta_full[2,1]", "beta_full[3,1]")]
  control_probs <- convert_logit_2_prob(control_logit_params)
  treatment_logit_params <- codafit_stan_draws[s,c("beta_full[1,2]", "beta_full[2,2]", "beta_full[3,2]")]
  treatment_probs <- convert_logit_2_prob(control_logit_params + treatment_logit_params)
  

  for(n in 1:N){
    u <- quantile_vec[n]
    if(u < 0.2*control_probs[1]){ # we need to compute the location parameter by understanding where to cut the lognormal off
      control_pareto_location_1 <- - qlnorm(1-(u/(0.2*control_probs[1])), meanlog = control_logmean_1, sdlog = control_logsd_1, lower.tail = TRUE, log.p = FALSE)
      quantile_value_output_control[s,n] <- - qpareto(1-(u/(0.2*control_probs[1])), a = control_pareto_scale_1, b = -control_pareto_location_1, lower.tail = TRUE, log.p = FALSE)
      } # i am worried there is a bug here 
    if(0.2*control_probs[1] <= u & u < control_probs[1]){ quantile_value_output_control[s,n]  <-  - qlnorm(1-(u/control_probs[1]), meanlog = control_logmean_1, sdlog = control_logsd_1, lower.tail = TRUE, log.p = FALSE)  }
    if(control_probs[1] < u & u < (control_probs[1]+control_probs[2])){ quantile_value_output_control[s,n]  <-  0 }
    if(u > (control_probs[1]+control_probs[2])){ quantile_value_output_control[s,n]  <-  qlnorm(((u-(control_probs[1]+control_probs[2]))/(1-control_probs[1]-control_probs[2])),  meanlog = control_logmean_2, sdlog = control_logsd_2, lower.tail = TRUE, log.p = FALSE)  }
    if( u >  1 - 0.2*(1-(control_probs[1]+control_probs[2]) )){ # we need to compute the location parameter
      control_pareto_location_2 <- qlnorm(  1 - 0.2*(1-(control_probs[1]+control_probs[2]) )   ,  meanlog = control_logmean_2, sdlog = control_logsd_2, lower.tail = TRUE, log.p = FALSE)    
      quantile_value_output_control[s,n]  <-   qpareto(1 - 0.2*(1-(control_probs[1]+control_probs[2]) ), a = control_pareto_scale_2, b = control_pareto_location_2, lower.tail = TRUE, log.p = FALSE)            
      }
    
    if(u < 0.2*treatment_probs[1]){ # we need to compute the location parameter by understanding where to cut the lognormal off
      treatment_pareto_location_1 <- - qlnorm(1-(u/(0.2*treatment_probs[1])), meanlog = treatment_logmean_1, sdlog = treatment_logsd_1, lower.tail = TRUE, log.p = FALSE)
      quantile_value_output_treatment[s,n] <- - qpareto(1-(u/(0.2*treatment_probs[1])), a = treatment_pareto_scale_1, b = -treatment_pareto_location_1, lower.tail = TRUE, log.p = FALSE)
    } # i am worried there is a bug here
    
    if(0.2*treatment_probs[1] <= u & u < treatment_probs[1]){ quantile_value_output_treatment[s,n]  <-  -qlnorm(1-(u/treatment_probs[1]), meanlog = treatment_logmean_1, sdlog = treatment_logsd_1, lower.tail = TRUE, log.p = FALSE) }
    if(treatment_probs[1] < u & u < (treatment_probs[1]+treatment_probs[2])){ quantile_value_output_treatment[s,n]  <-  0 }
    if(u > (treatment_probs[1]+treatment_probs[2])){ quantile_value_output_treatment[s,n]  <-  qlnorm(((u-(treatment_probs[1]+treatment_probs[2]))/(1-treatment_probs[1]-treatment_probs[2])), meanlog = treatment_logmean_2, sdlog = treatment_logsd_2, lower.tail = TRUE, log.p = FALSE)  }
    if( u >  1 - 0.2*(1-(treatment_probs[1]+treatment_probs[2]) )){ # we need to compute the location parameter
      treatment_pareto_location_2 <- qlnorm(  1 - 0.2*(1-(treatment_probs[1]+treatment_probs[2]) )   ,  meanlog = treatment_logmean_2, sdlog = treatment_logsd_2, lower.tail = TRUE, log.p = FALSE)    
      quantile_value_output_treatment[s,n]  <-   qpareto(1 - 0.2*(1-(treatment_probs[1]+treatment_probs[2]) ), a = treatment_pareto_scale_2, b = treatment_pareto_location_2, lower.tail = TRUE, log.p = FALSE)            
    }
    
  } # close forloop indexed by n 
  
} # closes forloop indexed by s 


quantile_value_output_difference <-  quantile_value_output_treatment -  quantile_value_output_control
quantiles_output <- data.frame(quantile_value_output_control, quantile_value_output_treatment, quantile_value_output_difference)
quantiles_output <- quantiles_output[complete.cases(quantiles_output),]
quantiles_output <- quantiles_output[is.finite(rowSums(quantiles_output)),]

clean_mean <- function(x){ mean(x, na.rm= TRUE)}
quantile_value_output_control_mean <- apply(quantiles_output[,1:10],2,clean_mean)
quantile_value_output_treatment_mean <- apply(quantiles_output[,11:20],2,clean_mean)
profit_quantile_value_output_difference_mean <- apply(quantiles_output[21:30],2,clean_mean)

quantile_value_output_difference <- quantile_value_output_difference[is.finite(rowSums(quantile_value_output_difference)),]
quantile_value_output_difference <- quantile_value_output_difference[complete.cases(quantile_value_output_difference),]
profit_recovered_quantile_differences <- apply(quantile_value_output_control,2,get_posterior_intervals_function)
profit_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)
profit_recovered_quantile_treatment <- apply(quantile_value_output_treatment,2,get_posterior_intervals_function)



### NOW THE GRAPHICS ###

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
install.packages("data.table")
library(data.table)
partial_pooling_beta_1_results_profit <- rbindlist(list(median_beta_1_profit_partial_pooling, ci_beta_1_profit_partial_pooling))
colnames(partial_pooling_beta_1_results_profit) <- c( "5th", "15th", "25th", "35th", "45th", "55th", "65th", "75th", "85th", "95th")

saveRDS(partial_pooling_beta_1_results_profit, "output/partial_pooling_beta_1_table_profit_composite_tail.rds")



posterior_beta_data <- sign(posterior_beta_data)*sqrt(abs(posterior_beta_data))
posterior_beta_data_plot <- ggplot(posterior_beta_data, aes(posterior_beta_data$quantiles_list))
pdf("output/posterior_parent_quantile_TEs_profit_composite.pdf", width=6.5, height=6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "red", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "red", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"50%"]), color = "maroon", size = 1.5) +
  ggtitle("Posterior quantile effects on profit from the Composite Tail") +
  theme(plot.title = element_text(size = 16)) + # xlim(0.05,0.95) + ylim(-100,100000) +
  xlab("Quantiles") + ylab("Quantile treatment effect")+
  theme(axis.text = element_text(size=14)) +  theme(axis.title.y = element_text(size = 14)) +  theme(axis.title.x = element_text(size = 14))
dev.off()
