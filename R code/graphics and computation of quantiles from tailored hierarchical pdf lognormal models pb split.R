# graphics and quantile computation for the pb split lognormal
# Rachael Meager
# First Version: Sept 2015

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
                  "MCMCpack", "geoR", "gtools", 'gPdtest', 'fBasics',"PtProcess", "VGAM", "MASS","quantreg",
                  "boot")
if(installation_needed){install.packages(package_list, repos='http://cran.us.r-project.org')}
if(loading_needed){lapply(package_list, require, character.only = TRUE)}


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

data_splitter_profit <- function(data){
  # check that the data is the right format
  if(!"site" %in% colnames(data) | !"treatment" %in% colnames(data) | !"profit" %in% colnames(data)  ){stop("data must be a dataframe with columns including site, treatment and profit")}
  
  N <- length(data$profit) # number of draws from the distribution / dimensionality of data vector
  M <- 3 # mixture components
  K <- length(unique(data$site)) # sites
  P <- 2 # dimensions of X 
  X <- cbind(rep(1,N), data$treatment)
  
  # create categorical allocations
  cat <- rep(NA,N) # storage
  for(i in 1:N){
    if(data$profit[i] < 0){ cat[i] <- 1  }
    else if(identical(data$profit[i],0)){cat[i] <- 2}
    else{cat[i] <- 3}
  }
  
  data <- data.frame(data, cat)
  
  data_split <- list( data[cat==1,],data[cat==3,])
  return(data_split)} # end data splitter profit function

data_splitter_expenditures <- function(data){
  # check that the data is the right format
  if(!"site" %in% colnames(data) | !"treatment" %in% colnames(data) | !"expenditures" %in% colnames(data)  ){stop("data must be a dataframe with columns including site, treatment and  expenditures")}
  

N <- length(data$expenditures) # number of draws from the distribution / dimensionality of data vector
M <- 2 # mixture components
K <- length(unique(data$site)) # sites
P <- 2 # dimensions of X 
X <- cbind(rep(1,N), data$treatment)

# create categorical allocations
cat <- rep(NA,N) # storage
for(i in 1:N){
  if(data$expenditures[i] < 0)stop("negative expenditures entry detected")
  if(identical(data$expenditures[i],0)){cat[i] <- 1}
  else{cat[i] <- 2}
}

data <- data.frame(data, cat)

data_split <- data[cat==2,] 
return(data_split)}

data_splitter_revenues <- function(data){
  # check that the data is the right format
  if(!"site" %in% colnames(data) | !"treatment" %in% colnames(data) | !"revenues" %in% colnames(data)  ){stop("data must be a dataframe with columns including site, treatment and revenues")}
  
N <- length(data$revenues) # number of draws from the distribution / dimensionality of data vector
M <- 2 # mixture components
K <- length(unique(data$site)) # sites
P <- 2 # dimensions of X 
X <- cbind(rep(1,N), data$treatment)

# create categorical allocations
cat <- rep(NA,N) # storage
for(i in 1:N){
  if(data$revenues[i] < 0)stop("negative revenues entry detected")
  if(identical(data$revenues[i],0)){cat[i] <- 1}
  else{cat[i] <- 2}
}

data <- data.frame(data, cat)

data_split <- data[cat==2,]
return(data_split)}

### Load data ###

load("output/microcredit_profit_tailored_hierarchical_pdf_output_pb_split_lognormal.RData")
stan_fit_table_0 <- xtable(stan_fit_summary_0$summary)
stan_fit_table_1 <- xtable(stan_fit_summary_1$summary)
stan_fit_table_profit_0 <- stan_fit_table_0
stan_fit_profit_0 <- stan_fit_pb_0
codafit_stan_draws_profit_0 <- as.matrix(stan2coda(stan_fit_profit_0))
data_split_profit_0 <- data_splitter_profit(data_0)
data_profit_0 <- data_0
stan_fit_table_profit_1 <- stan_fit_table_1
stan_fit_profit_1 <- stan_fit_pb_1
codafit_stan_draws_profit_1 <- as.matrix(stan2coda(stan_fit_profit_1))
data_split_profit_1 <- data_splitter_profit(data_1)
data_profit_1 <- data_1

load("output/microcredit_expenditures_tailored_hierarchical_pdf_output_pb_split_lognormal.RData")
stan_fit_table_0 <- xtable(stan_fit_summary_0$summary)
stan_fit_table_1 <- xtable(stan_fit_summary_1$summary)
stan_fit_table_expenditures_0 <- stan_fit_table_0
stan_fit_expenditures_0 <- stan_fit_pb_0
codafit_stan_draws_expenditures_0 <- as.matrix(stan2coda(stan_fit_expenditures_0))
data_split_expenditures_0 <- data_splitter_expenditures(data_0)
data_expenditures_0 <- data_0
stan_fit_table_expenditures_1 <- stan_fit_table_1
stan_fit_expenditures_1 <- stan_fit_pb_1
codafit_stan_draws_expenditures_1 <- as.matrix(stan2coda(stan_fit_expenditures_1))
data_split_expenditures_1 <- data_splitter_expenditures(data_1)
data_expenditures_1 <- data_1


load("output/microcredit_revenues_tailored_hierarchical_pdf_output_pb_split_lognormal.RData")
stan_fit_table_0 <- xtable(stan_fit_summary_0$summary)
stan_fit_table_1 <- xtable(stan_fit_summary_1$summary)
stan_fit_table_revenues_0 <- stan_fit_table_0
stan_fit_revenues_0 <- stan_fit_pb_0
codafit_stan_draws_revenues_0 <- as.matrix(stan2coda(stan_fit_revenues_0))
data_split_revenues_0 <- data_splitter_revenues(data_0)
data_revenues_0 <- data_0
stan_fit_table_revenues_1 <- stan_fit_table_1
stan_fit_revenues_1 <- stan_fit_pb_1
codafit_stan_draws_revenues_1 <- as.matrix(stan2coda(stan_fit_revenues_1))
data_split_revenues_1 <- data_splitter_revenues(data_1)
data_revenues_1 <- data_1

### FOR PB = 1 ###

## FOR THE 1 TAIL OUTCOMES

## revenues! 

codafit_stan_draws <- codafit_stan_draws_revenues_1

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
  
  # a question for me: does this structure REALLY propagate like the Gaussian structure does? 
  # One might ask "Should I have been doing everything on the log scale?" but this will not save you in a mixture model. 
  
  
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


# for expenditures
codafit_stan_draws <- codafit_stan_draws_expenditures_1
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
  
  # a question for me: does this structure REALLY propagate like the Gaussian structure does? 
  # One might ask "Should I have been doing everything on the log scale?" but this will not save you in a mixture model. 
  
  
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


### FOR PROFIT


codafit_stan_draws <- codafit_stan_draws_profit_1
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
    if(u < control_probs[1]){ quantile_value_output_control[s,n]  <-  - qlnorm((u/control_probs[1]), meanlog = control_logmean_1, sdlog = control_logsd_1, lower.tail = TRUE, log.p = FALSE) }
    if(control_probs[1] < u & u < (control_probs[1]+control_probs[2])){ quantile_value_output_control[s,n]  <-  0 }
    if(u > (control_probs[1]+control_probs[2])){ quantile_value_output_control[s,n]  <-  qlnorm(((u-(control_probs[1]+control_probs[2]))/(1-control_probs[1]-control_probs[2])),  meanlog = control_logmean_2, sdlog = control_logsd_2, lower.tail = TRUE, log.p = FALSE)  }
    
    if(u < treatment_probs[1]){ quantile_value_output_treatment[s,n]  <-  -qlnorm((u/treatment_probs[1]), meanlog = treatment_logmean_1, sdlog = treatment_logsd_1, lower.tail = TRUE, log.p = FALSE) }
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
pdf("output/posterior_parent_quantile_TEs_profit_pb_1_lognormal.pdf", width=6.5, height=6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "orange", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "orange", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"50%"]), color = "dark orange", size = 1.5) +
  ggtitle("Posterior quantile effects on profit PB=1 ") +
  theme(plot.title = element_text(size = 16)) + xlim(0.05,0.95) + ylim(-350,350) +
  xlab("Quantiles") + ylab("Quantile treatment effect ")+
  theme(axis.text = element_text(size=14)) +  theme(axis.title.y = element_text(size = 14)) +  theme(axis.title.x = element_text(size = 14))
dev.off()

# expenditures # 
# make ribbon plot for posterior 
quantiles_list <- quantile_vec
posterior_beta_data <- t(rbind(expenditures_recovered_quantile_differences, expenditures_quantile_value_output_difference_mean  ))
posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")


posterior_beta_data_plot <- ggplot(posterior_beta_data, aes(posterior_beta_data$quantiles_list))
pdf("output/posterior_parent_quantile_TEs_expenditures_pb_1_lognormal.pdf", width=6.5, height=6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "orange", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "orange", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"50%"]), color = "dark orange", size = 1.5) +
  ggtitle("Posterior quantile effects on expenditures PB=1") +
  theme(plot.title = element_text(size = 16)) + xlim(0.05,0.95) +  ylim(-100,500) +
  xlab("Quantiles") + ylab("Quantile treatment effect ")+
  theme(axis.text = element_text(size=14)) +  theme(axis.title.y = element_text(size = 14)) +  theme(axis.title.x = element_text(size = 14))
dev.off()


# revenues # 
# make ribbon plot for posterior 
quantiles_list <- quantile_vec
posterior_beta_data <- t(rbind(revenues_recovered_quantile_differences, revenues_quantile_value_output_difference_mean  ))
posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")


posterior_beta_data_plot <- ggplot(posterior_beta_data, aes(posterior_beta_data$quantiles_list))
pdf("output/posterior_parent_quantile_TEs_revenues_pb_1_lognormal.pdf", width=6.5, height=6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "orange", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "orange", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"50%"]), color = "dark orange", size = 1.5) +
  ggtitle("Posterior quantile effects on revenues PB=1") +
  theme(plot.title = element_text(size = 16)) + xlim(0.05,0.95) +  ylim(-250,600) +
  xlab("Quantiles") + ylab("Quantile treatment effect ")+
  theme(axis.text = element_text(size=14)) +  theme(axis.title.y = element_text(size = 14)) +  theme(axis.title.x = element_text(size = 14))
dev.off()



### FOR PB = 0 ###

## FOR THE 1 TAIL OUTCOMES

## revenues! 

codafit_stan_draws <- codafit_stan_draws_revenues_0

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
  
  # a question for me: does this structure REALLY propagate like the Gaussian structure does? 
  # One might ask "Should I have been doing everything on the log scale?" but this will not save you in a mixture model. 
  
  
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


# for expenditures
codafit_stan_draws <- codafit_stan_draws_expenditures_0
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
  
  # a question for me: does this structure REALLY propagate like the Gaussian structure does? 
  # One might ask "Should I have been doing everything on the log scale?" but this will not save you in a mixture model. 
  
  
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


### FOR PROFIT


codafit_stan_draws <- codafit_stan_draws_profit_0
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
    if(u < control_probs[1]){ quantile_value_output_control[s,n]  <-  - qlnorm((u/control_probs[1]), meanlog = control_logmean_1, sdlog = control_logsd_1, lower.tail = TRUE, log.p = FALSE) }
    if(control_probs[1] < u & u < (control_probs[1]+control_probs[2])){ quantile_value_output_control[s,n]  <-  0 }
    if(u > (control_probs[1]+control_probs[2])){ quantile_value_output_control[s,n]  <-  qlnorm(((u-(control_probs[1]+control_probs[2]))/(1-control_probs[1]-control_probs[2])),  meanlog = control_logmean_2, sdlog = control_logsd_2, lower.tail = TRUE, log.p = FALSE)  }
    
    if(u < treatment_probs[1]){ quantile_value_output_treatment[s,n]  <-  -qlnorm((u/treatment_probs[1]), meanlog = treatment_logmean_1, sdlog = treatment_logsd_1, lower.tail = TRUE, log.p = FALSE) }
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
pdf("output/posterior_parent_quantile_TEs_profit_pb_0_lognormal.pdf", width=6.5, height=6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "blue", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "blue", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"50%"]), color = "dark blue", size = 1.5) +
  ggtitle("Posterior quantile effects on profit PB=0 (cubed root) ") +
  theme(plot.title = element_text(size = 16)) + xlim(0.05,0.95) +  ylim(-350,350) +
  xlab("Quantiles") + ylab("Quantile treatment effect ")+
  theme(axis.text = element_text(size=14)) +  theme(axis.title.y = element_text(size = 14)) +  theme(axis.title.x = element_text(size = 14))
dev.off()

# expenditures # 
# make ribbon plot for posterior 
quantiles_list <- quantile_vec
posterior_beta_data <- t(rbind(expenditures_recovered_quantile_differences, expenditures_quantile_value_output_difference_mean  ))
posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")


posterior_beta_data_plot <- ggplot(posterior_beta_data, aes(posterior_beta_data$quantiles_list))
pdf("output/posterior_parent_quantile_TEs_expenditures_pb_0_lognormal.pdf", width=6.5, height=6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "blue", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "blue", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"50%"]), color = "dark blue", size = 1.5) +
  ggtitle("Posterior quantile effects on expenditures PB=0 ") +
  theme(plot.title = element_text(size = 16)) + xlim(0.05,0.95) +  ylim(-100,500) +
  xlab("Quantiles") + ylab("Quantile treatment effect ")+
  theme(axis.text = element_text(size=14)) +  theme(axis.title.y = element_text(size = 14)) +  theme(axis.title.x = element_text(size = 14))
dev.off()


# revenues # 
# make ribbon plot for posterior 
quantiles_list <- quantile_vec
posterior_beta_data <- t(rbind(revenues_recovered_quantile_differences, revenues_quantile_value_output_difference_mean  ))
posterior_beta_data <- data.frame(posterior_beta_data, quantiles_list)
colnames(posterior_beta_data) <- c("2.5%","25%","50%", "75%", "97.5%", "mean", "quantiles_list")


posterior_beta_data_plot <- ggplot(posterior_beta_data, aes(posterior_beta_data$quantiles_list))
pdf("output/posterior_parent_quantile_TEs_revenues_pb_0_lognormal.pdf", width=6.5, height=6)
posterior_beta_data_plot +
  geom_ribbon(aes(ymin = posterior_beta_data[,"2.5%"], ymax = posterior_beta_data[,"97.5%"]), fill = "blue", alpha=0.3) +
  geom_ribbon(aes(ymin = posterior_beta_data[,"25%"], ymax = posterior_beta_data[,"75%"]), fill = "blue", alpha=0.6) +
  geom_line(aes(y = posterior_beta_data[,"50%"]), color = "dark blue", size = 1.5) +
  ggtitle("Posterior quantile effects on revenues PB=0") +
  theme(plot.title = element_text(size = 16)) + xlim(0.05,0.95) + ylim(-250,600) +
  xlab("Quantiles") + ylab("Quantile treatment effect ")+
  theme(axis.text = element_text(size=14)) +  theme(axis.title.y = element_text(size = 14)) +  theme(axis.title.x = element_text(size = 14))
dev.off()


