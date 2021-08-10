# microcredit tailored hierarchical pdf for consumerdurables with lognormal tails
# Rachael Meager
# April 2016 

### Notes

### Preliminaries 



# clear the workspace to avoid gremlins and past globals from past irresponsible scripts
# but we can't do this if the masterfile is being used to run the script, so we check that first:
if(exists("masterfile_run") == "FALSE"){
  rm(list = ls())
}


# install and load packages
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
if(installation_needed){for (p in packagelist) {installPackage(p, repos='http://cran.us.r-project.org')}}
if(loading_needed){lapply(package_list, require, character.only = TRUE)}

# configure rstan
print("loaded packages")
rstan_options(auto_write = TRUE)
print("ran rstan_options command")
cores <- as.numeric(Sys.getenv('NSLOTS'))
print("cores variable assigned")
options(mc.cores = cores)
print("mc.cores = cores")
# load data

# load data

load("data/microcredit_project_data.RData")
print("loaded data")
if( exists("USD_convert_to_2009_dollars")!=TRUE){stop("OLD DATA FILE LOADED! RELOAD CORRECT FILE!")}
# This preps data so that the model can be written efficiently - it's not necessary but it's the best way I know to code it in STAN
site <- c(attanasio_indicator, augsberg_indicator,banerjee_indicator, crepon_indicator)-1
consumerdurables <- c(attanasio_consumerdurables, augsberg_consumerdurables,banerjee_consumerdurables, crepon_consumerdurables)
treatment <- c(attanasio_treatment, augsberg_treatment,banerjee_treatment, crepon_treatment)

# Now we have to standardise any variables which are in local currency units to USD PPP per fortnight in 2009 dollars


expanded_standardiser_USD_PPP_per_fortnight <- c( rep(the_consumerdurables_standardiser_USD_PPP_per_fortnight[1],length(attanasio_indicator)),
                                                  rep(the_consumerdurables_standardiser_USD_PPP_per_fortnight[2],length(augsberg_indicator)),
                                                  rep(the_consumerdurables_standardiser_USD_PPP_per_fortnight[3],length(banerjee_indicator)),
                                                  rep(the_consumerdurables_standardiser_USD_PPP_per_fortnight[4],length(crepon_indicator)))

consumerdurables <- consumerdurables*expanded_standardiser_USD_PPP_per_fortnight

# bind everything into a data frame
data <- data.frame(site, consumerdurables, treatment)
# We gotta remove the NA values
data <- data[complete.cases(data),]
print(dim(data))
#stan needs these declared as data
N <- length(data$consumerdurables) # number of draws from the distribution / dimensionality of data vector
M <- 2 # mixture components
K <- length(unique(data$site)) # sites
P <- 2 # dimensions of X 
X <- cbind(rep(1,N), data$treatment)

# create categorical allocations
cat <- rep(NA,N) # storage
for(i in 1:N){
  if(data$consumerdurables[i] < 0)stop("negative consumerdurables entry detected")
  if(identical(data$consumerdurables[i],0)){cat[i] <- 1}
  else{cat[i] <- 2}
}
print("category assignment for logit complete")
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
                 y_pos = data_split$consumerdurables,
                 N_pos = length(data_split$consumerdurables),
                 site_pos = data_split$site,
                 quantile_vec = seq(0.05, 0.95,0.1))

print("declared stan data")
stan_fit <- stan("stan-code/tailored-hierarchical-pdf-log-normal-1-tail.stan", 
                 iter = 5000, chains = 6, data = fit_data,
                 control = list(adapt_delta = 0.99, max_treedepth = 15))
print("ran stan allegedly")
print(stan_fit)

stan_fit_summary <- summary(stan_fit)
stan_fit_table <- xtable(stan_fit_summary$summary)

saveRDS(stan_fit_summary, file = "tailored_hierarchical_pdf_microcredit_output_1.RDS", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

save.image("output/microcredit_consumerdurables_lognormal_tailored_hierarchical_pdf_output_5000_iters.RData")

sink("output/microcredit_consumerdurables_lognormal_tailored_hierarchical_output_table_5000_iters.txt")
print(stan_fit)
dev.off()







