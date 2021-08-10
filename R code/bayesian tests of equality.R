



# Rachael Meager
# First Version: Feb 2020


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

# Load data
load("output/microcredit_profit_lognormal_tailored_hierarchical_pdf_output_5000_iters.RData")
stan_fit_table_profit <- stan_fit_table
stan_fit_profit <- stan_fit
codafit_stan_draws_profit <- as.matrix(stan2coda(stan_fit_profit))
data_split_profit <- data_split
data_profit <- data

### PROFIT AT THE PARENT LEVEL



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
dim(quantile_value_output_difference)

# now for the test
tol <- c(0.5,1,3,5)
T <- length(tol)
post_odds <- rep(NA,T)
for(t in 1:T){
  tol_run <- tol[t]
close_enough_indicator <- rep(NA, S)
for(s in 1:S){
  x <- quantile_value_output_difference[s,]
  close_enough_indicator[s] <- (abs(max(x) - min(x)) < tol_run)
}
tab_out <- table(close_enough_indicator)
post_odds[t] <- tab_out[2]/tab_out[1]
}



# now let's do it for the business ownership split
# independent models means I have got myself stuck with independent posterior draws. oh well. here we go!
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
quantile_value_output_difference_PB1 <- quantile_value_output_difference

quantile_value_output_control_mean <- apply(quantile_value_output_control,2,mean)
quantile_value_output_treatment_mean <- apply(quantile_value_output_treatment,2,mean)
profit_quantile_value_output_difference_mean <- apply(quantile_value_output_difference,2,mean)


profit_recovered_quantile_differences <- apply(quantile_value_output_difference,2,get_posterior_intervals_function)
profit_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)
profit_recovered_quantile_treatment <- apply(quantile_value_output_treatment,2,get_posterior_intervals_function)



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
quantile_value_output_difference_PB0 <- quantile_value_output_difference

quantile_value_output_control_mean <- apply(quantile_value_output_control,2,mean)
quantile_value_output_treatment_mean <- apply(quantile_value_output_treatment,2,mean)
profit_quantile_value_output_difference_mean <- apply(quantile_value_output_difference,2,mean)


profit_recovered_quantile_differences <- apply(quantile_value_output_difference,2,get_posterior_intervals_function)
profit_recovered_quantile_control <- apply(quantile_value_output_control,2,get_posterior_intervals_function)
profit_recovered_quantile_treatment <- apply(quantile_value_output_treatment,2,get_posterior_intervals_function)



quantiles_values_both_groups <- cbind(quantile_value_output_difference_PB0, quantile_value_output_difference_PB1)


# now for the test
# quantile by quantile at 4 tolerances
tol <- c(0.5,1,3,5)
T <- length(tol)
table_out <- matrix(NA, 10, T)
for(t in 1:T){
close_enough_indicator <- rep(NA, S)
close_enough_qbyq <- matrix(NA, S, 10)
for(s in 1:S){
  tol_run <- tol[t]
  x <- quantile_value_output_difference_PB0[s,]
  y <- quantile_value_output_difference_PB1[s,]
  for(i in 1:10){
  close_enough_qbyq[s,i] <- (abs(x[i] - y[i]) < tol_run)
  }
}


for(i in 1:10){
table_out[i,t] <- (mean(close_enough_qbyq[,i]))}

}
table_out <- as.data.frame(table_out)
table_out <- data.frame(quantiles, table_out)
colnames(table_out) <- c("Quantile", "Tolerance = 0.5 USD PPP", "1 USD PPP", "3 USD PPP", "5 USD PPP")
stargazer(table_out, summary = FALSE, title = "Posterior odds that effects are equal for groups PB = 0 and 1")

