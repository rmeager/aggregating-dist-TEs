# full tailored hierarchical pdf model for temptation 
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
site <- c(angelucci_indicator, attanasio_indicator, augsberg_indicator,banerjee_indicator, crepon_indicator)
temptation <- c(angelucci_temptation, attanasio_temptation, augsberg_temptation,banerjee_temptation, crepon_temptation)
treatment <- c(angelucci_treatment, attanasio_treatment, augsberg_treatment,banerjee_treatment, crepon_treatment)

priorbiz <- c(angelucci_existingbusiness, attanasio_existingbusiness, augsberg_existingbusiness,banerjee_existingbusiness, crepon_existingbusiness)

# Now we have to standardise any variables which are in local currency units to USD PPP per fortnight in 2009 dollars


expanded_standardiser_USD_PPP_per_fortnight <- c( rep(the_temptation_standardiser_USD_PPP_per_fortnight[1],length(angelucci_indicator)),
                                                  rep(the_temptation_standardiser_USD_PPP_per_fortnight[2],length(attanasio_indicator)),
                                                  rep(the_temptation_standardiser_USD_PPP_per_fortnight[3],length(augsberg_indicator)),
                                                  rep(the_temptation_standardiser_USD_PPP_per_fortnight[4],length(banerjee_indicator)),
                                                  rep(the_temptation_standardiser_USD_PPP_per_fortnight[5],length(crepon_indicator)))

temptation <- temptation*expanded_standardiser_USD_PPP_per_fortnight

# bind everything into a data frame
data <- data.frame(site, temptation, treatment, priorbiz)

# We gotta remove the NA values
data <- data[complete.cases(data),]

tailored_hierarchical_pdf_aggregator <- function(data){
  
  # check that the data is the right format
  if(!"site" %in% colnames(data) | !"treatment" %in% colnames(data) | !"temptation" %in% colnames(data)  ){stop("data must be a dataframe with columns including site, treatment and temptation")}
  

#stan needs these declared as data
N <- length(data$temptation) # number of draws from the distribution / dimensionality of data vector
M <- 2 # mixture components
K <- length(unique(data$site)) # sites
P <- 2 # dimensions of X 
X <- cbind(rep(1,N), data$treatment)

# create categorical allocations
cat <- rep(NA,N) # storage
for(i in 1:N){
  if(data$temptation[i] < 0)stop("negative temptation entry detected")
  if(identical(data$temptation[i],0)){cat[i] <- 1}
  else{cat[i] <- 2}
}

data <- data.frame(data, cat)

data_split <- data[cat==2,]


# now STAN MODEL! 

fit_data <- list(M = M,
                 N = N,
                 K = K,
                 P = P, 
                 Q = 10,
                 cat = cat,
                 site = data$site,
                 x = X,
                 treatment_pos = data_split$treatment,
                 y_pos = data_split$temptation,
                 N_pos = length(data_split$temptation),
                 site_pos = data_split$site,
                 quantile_vec = seq(0.05, 0.95,0.1))


stan_fit <- stan("stan-code/tailored-hierarchical-pdf-log-normal-1-tail.stan", iter = 5000, cores = 8, data = fit_data)

print(stan_fit)
return(stan_fit)
} # close function

data_0 <- data[data$priorbiz==0,]
data_0$site[data_0$site==7] <- 6
stan_fit_pb_0 <- tailored_hierarchical_pdf_aggregator(data_0)

data_1 <- data[data$priorbiz==1,]
stan_fit_pb_1 <- tailored_hierarchical_pdf_aggregator(data_1)

stan_fit_summary_0 <- summary(stan_fit_pb_0)
stan_fit_table_0 <- xtable(stan_fit_summary_0$summary)

save.image("output/microcredit_temptation_tailored_hierarchical_pdf_output_pb_split_lognormal.RData")


saveRDS(stan_fit_summary_0, file = "output/tailored_hierarchical_pdf_temptation_microcredit_output_pb_split_0_lognormal.RDS", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

stan_fit_summary_1 <- summary(stan_fit_pb_1)
stan_fit_table_1 <- xtable(stan_fit_summary_1$summary)

saveRDS(stan_fit_summary_1, file = "output/tailored_hierarchical_pdf_temptation_microcredit_output_pb_split_1_lognormal.RDS", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)


save.image("output/microcredit_temptation_tailored_hierarchical_pdf_output_pb_split_lognormal.RData")



