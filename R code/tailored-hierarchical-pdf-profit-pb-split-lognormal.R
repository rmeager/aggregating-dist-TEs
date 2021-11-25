# microcredit tailored hierarchical pdf for profit lognormal model
# Rachael Meager
# April 2016 

### Notes

### Preliminaries 

# clear the workspace to avoid gremlins and past globals from past irresponsible scripts
if(exists("masterfile_run") == "FALSE"){
  rm(list = ls())
}


# install and load packages

installPackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {
    print(paste('installing', p))
    install.packages(p, dep = TRUE)
  } else {
    print(paste(p, 'is already installed'))
  }
  require(p, character.only = TRUE) # test if it loads
}
installation_needed  <- FALSE
loading_needed <- TRUE
package_list <- c('ggplot2', 'rstan','reshape','reshape2','coda','xtable', 'dplyr', 'Runuran', 'testthat',
                  "MCMCpack", "geoR", "gtools", 'gPdtest', 'fBasics',"PtProcess", "VGAM")
if(installation_needed){install.packages(package_list, repos='http://cran.us.r-project.org')}
if(loading_needed){lapply(package_list, require, character.only = TRUE)}

# configure rstan

rstan_options(auto_write = TRUE)

cores <- as.numeric(Sys.getenv('NSLOTS'))

options(mc.cores = cores)

# load data

load("data/microcredit_project_data.RData")

if( exists("USD_convert_to_2009_dollars")!=TRUE){stop("OLD DATA FILE LOADED! RELOAD CORRECT FILE!")}
# This preps data so that the model can be written efficiently - it's not necessary but it's the best way I know to code it in STAN
site <- c(angelucci_indicator, attanasio_indicator, augsberg_indicator,banerjee_indicator, crepon_indicator, karlan_indicator, tarozzi_indicator)
profit <- c(angelucci_profit, attanasio_profit, augsberg_profit,banerjee_profit, crepon_profit, karlan_profit, tarozzi_profit)
treatment <- c(angelucci_treatment, attanasio_treatment, augsberg_treatment,banerjee_treatment, crepon_treatment, karlan_treatment, tarozzi_treatment)
priorbiz <- c(angelucci_existingbusiness, attanasio_existingbusiness, augsberg_existingbusiness,banerjee_existingbusiness, crepon_existingbusiness, karlan_existingbusiness, tarozzi_existingbusiness)

# Now we have to standardise any variables which are in local currency units to USD PPP per fortnight in 2009 dollars


expanded_standardiser_USD_PPP_per_fortnight <- c( rep(the_profit_standardiser_USD_PPP_per_fortnight[1],length(angelucci_indicator)),
                                                  rep(the_profit_standardiser_USD_PPP_per_fortnight[2],length(attanasio_indicator)),
                                                  rep(the_profit_standardiser_USD_PPP_per_fortnight[3],length(augsberg_indicator)),
                                                  rep(the_profit_standardiser_USD_PPP_per_fortnight[4],length(banerjee_indicator)),
                                                  rep(the_profit_standardiser_USD_PPP_per_fortnight[5],length(crepon_indicator)),
                                                  rep(the_profit_standardiser_USD_PPP_per_fortnight[6],length(karlan_indicator)),
                                                  rep(the_profit_standardiser_USD_PPP_per_fortnight[7],length(tarozzi_indicator)))

profit <- profit*expanded_standardiser_USD_PPP_per_fortnight

# bind everything into a data frame
data <- data.frame(site, profit, treatment, priorbiz)

# We gotta remove the NA values
data <- data[complete.cases(data),]



tailored_hierarchical_pdf_aggregator <- function(data){

  # check that the data is the right format
  if(!"site" %in% colnames(data) | !"treatment" %in% colnames(data) | !"profit" %in% colnames(data)  ){stop("data must be a dataframe with columns including site, treatment and profit")}
  
#stan needs these declared as data
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

# now STAN MODEL! 


fit_data <- list(M = 3,
                 N = N,
                 K = K,
                 P = 2, 
                 Q = 10,
                 cat = cat,
                 site = data$site,
                 x = X,
                 treatment_neg = data_split[[1]]$treatment,
                 y_neg = -data_split[[1]]$profit,
                 treatment_pos = data_split[[2]]$treatment,
                 y_pos = data_split[[2]]$profit,
                 N_neg = length(data_split[[1]]$profit),
                 N_pos = length(data_split[[2]]$profit),
                 site_neg = data_split[[1]]$site,
                 site_pos = data_split[[2]]$site,
                 quantile_vec = seq(0.05, 0.95, .1))


stan_fit <- stan("stan-code/tailored-hierarchical-pdf-log-normal.stan", iter = 5000, cores = 8, data = fit_data)

print(stan_fit)

return(stan_fit)
} # close function

data_0 <- data[data$priorbiz==0,]
data_0$site[data_0$site==7] <- 6 #ensures no attempt to call element 7 of the resulting 6 dimensional vector
stan_fit_pb_0 <- tailored_hierarchical_pdf_aggregator(data_0)

data_1 <- data[data$priorbiz==1,]
stan_fit_pb_1 <- tailored_hierarchical_pdf_aggregator(data_1)

stan_fit_summary_0 <- summary(stan_fit_pb_0)
stan_fit_table_0 <- xtable(stan_fit_summary_0$summary)

saveRDS(stan_fit_summary_0, file = "output/tailored_hierarchical_pdf_microcredit_output_pb_split_0_lognormal.RDS", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

stan_fit_summary_1 <- summary(stan_fit_pb_1)
stan_fit_table_1 <- xtable(stan_fit_summary_1$summary)

saveRDS(stan_fit_summary_1, file = "output/tailored_hierarchical_pdf_microcredit_output_pb_split_1_lognormal.RDS", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

# now, because how I fit this, the stanfit object is massive
# to be able to upload this to github, i have to abandon it and work only with the reduced codafit object

codafit_stan_draws_profit_pb_1 <- as.matrix(stan2coda(stan_fit_pb_1))
saveRDS(codafit_stan_draws_profit_pb_1, file = "output/tailored_hierarchical_pdf_microcredit_profit_pb_split_1.RDS", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

codafit_stan_draws_profit_pb_0 <- as.matrix(stan2coda(stan_fit_pb_0))
saveRDS(codafit_stan_draws_profit_pb_1, file = "output/tailored_hierarchical_pdf_microcredit_profit_pb_split_0.RDS", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

save.image("output/microcredit_profit_tailored_hierarchical_pdf_output_pb_split_lognormal.RData")



