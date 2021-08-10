# microcredit tailored hierarchical pdf for profit run 
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
data <- data.frame(site, profit, treatment)

# We gotta remove the NA values
data <- data[complete.cases(data),]


# ok i have now generated the data

#stan needs these declared as data
M <- 3 # mixture components
P <- 2 # dimensions of X 
N <- length(data$site)
K <- length(unique(data$site))
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

# now find the quantile cut locations 

N_splits <- matrix(NA, K, M )
for(k in 1:K){
  for(m in 1:M){
    N_splits[k,m] <- as.numeric(length(data$profit[data$site==k & cat==m]))
  }
}

twenty_percent_of_N_splits <- round(0.2*N_splits,0)
cut_locations <- matrix(NA, K,2)
for(k in 1:K){
  cut_locations[k,1] <- quantile(data$profit[data$cat==1 & data$site==k],.2)
  cut_locations[k,2] <- quantile(data$profit[data$cat==3 & data$site==k],.8)
}

#manually override the NA, it doesnt matter as the NA signals no negative tail
cut_locations[3,1] <- -1000


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
                 cut_locations = abs(cut_locations),
                 quantile_vec = seq(0.05, 0.95, .1))

stan_fit <- stan("stan-code/tailored-hierarchical-pdf-composite-tails.stan", iter = 4000, chains = 6, data = fit_data)

print(stan_fit)

stan_fit_summary <- summary(stan_fit)
stan_fit_table <- xtable(stan_fit_summary$summary)

saveRDS(stan_fit_summary, file = "composite_tail_tailored_hierarchical_pdf_microcredit_output.RDS", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

save.image("output/microcredit_profit_composite_tails_tailored_hierarchical_pdf_output.RData")

sink("output/microcredit_profit_composite_tails_tailored_hierarchical_output_table.txt")
print(stan_fit)
dev.off()






