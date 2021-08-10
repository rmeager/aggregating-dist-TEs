



data {
int M; // number of categories
int N; // number of observations
int P; // dimentionality of beta parameter
int K; // number of sites 
int Q; // number of quantiles
int cat[N]; // category indicator
int N_neg; // number of obvs in negative tail
int N_pos; // number of obvs in positive tail
vector[P] x[N]; // covariates
int site[N]; // site indicator 
int treatment_neg[N_neg]; // treatment in neg tail
int treatment_pos[N_pos]; // treatment in pos tail
int site_neg[N_neg]; // site indicator in neg tail
int site_pos[N_pos]; // site indicator in pos tail
real y_neg[N_neg]; // the negative tail
real y_pos[N_pos]; // the positive tail 
matrix[K,2] cut_locations; // the split between pareto and lognormal
vector[Q] quantile_vec ; // vector of quantiles of interest
}

parameters {

real mu[2]; #lognormals control locations
real tau[2]; #lognormals location TEs
real pareto_control_shape[2]; #pareto scale or shape parameter for control groups
real pareto_tau_shape[2]; #pareto scale or shape TE
real<lower=0> sd_mu[2]; #hypervariance of mu
real<lower=0> sd_tau[2];#hypervariance of tau
real<lower=0> sd_pareto_control_shape[2]; #hypervariance of pareto scales
real<lower=0> sd_pareto_tau_shape[2]; #hypervariance of pareto scales
real sigma_control[2]; #lognormals control scales
real sigma_TE[2]; #lognormals scales TEs
real<lower=0> sd_sigma_control[2]; #hypervariance of lognormal control scales
real<lower=0> sd_sigma_TE[2]; #hypervariance of lognormal scale TE
matrix[K,2] mu_k;
matrix[K,2] tau_k;
matrix[K,2] sigma_control_k;
matrix[K,2] sigma_TE_k;
matrix[K,2] pareto_control_shape_k;
matrix[K,2] pareto_tau_shape_k;
matrix[M-1,P] beta; // the parent parameters minus the Mth category
matrix<lower=0>[M,P] sigma; // the set of M*P parent variances (not a covariance matrix)
matrix[M,P] beta_k_raw[K]; // the hierarchical increments 

}
transformed parameters{
matrix[M,P] beta_full;
matrix[M,P] beta_k[K];
beta_full <- append_row(beta,rep_row_vector(0, P));
for (m in 1:M){
 for (k in 1:K){
  for (p in 1:P){
   beta_k[k,m,p] <- beta_full[m,p] + sigma[m,p]*beta_k_raw[k,m,p];
}}}
}
model {
to_vector(beta) ~ normal(0,5);
to_vector(sigma) ~ cauchy(0,2);
for (m in 1:M){
  for (k in 1:K){
    beta_k_raw[k,m] ~ normal(0,1);
  }}
for (n in 1:N)
cat[n] ~ categorical_logit(beta_k[site[n]] * x[n]);

mu ~ normal(0,100);
tau ~ normal(0,100);
sd_mu ~ cauchy(0,2);
sd_tau ~ cauchy(0,2);
sigma_control ~ normal(0,100);
sd_sigma_control ~ cauchy(0,2);
sigma_TE ~ normal(0,100);
sd_sigma_TE ~ cauchy(0,2);
pareto_tau_shape ~ normal(0,5);
sd_pareto_tau_shape ~ cauchy(0,2);
pareto_control_shape ~ normal(2,2);
sd_pareto_control_shape ~ cauchy(0,2);
for (k in 1:K){
  for (i in 1:2){
  mu_k[k,i] ~ normal(mu[i],sd_mu[i]);
  tau_k[k,i] ~ normal(tau[i], sd_tau[i]);
  pareto_tau_shape_k[k,i] ~ normal(pareto_tau_shape[i],sd_pareto_tau_shape[i]);
  pareto_control_shape_k[k,i] ~ normal(pareto_control_shape[i], sd_pareto_control_shape[i]);
  sigma_control_k[k,i] ~ normal(sigma_control[i],sd_sigma_control[i]);
  sigma_TE_k[k,i] ~ normal(sigma_TE[i], sd_sigma_TE[i]);
  }
}
for (n in 1:N_neg){
  if(y_neg[n] < cut_locations[site_neg[n],1]){
    y_neg[n] ~ lognormal(mu_k[site_neg[n],1] + tau_k[site_neg[n],1]*treatment_neg[n], exp(sigma_control_k[site_neg[n],1] + sigma_TE_k[site_neg[n],1]*treatment_neg[n]));
    }
    else{ y_neg[n] ~ pareto((cut_locations[site_neg[n] ,1]),exp(pareto_control_shape_k[site_neg[n],1]+pareto_tau_shape_k[site_neg[n],1]*treatment_neg[n]) );}
    
}
for(n in 1:N_pos){
  if(y_pos[n] < cut_locations[site_pos[n],2]){
    y_pos[n] ~ lognormal(mu_k[site_pos[n],2] + tau_k[site_pos[n],2]*treatment_pos[n], exp(sigma_control_k[site_pos[n],2] + sigma_TE_k[site_pos[n],2]*treatment_pos[n])) ;
    }
    else{y_pos[n] ~ pareto((cut_locations[site_pos[n] ,2]),exp(pareto_control_shape_k[site_pos[n],2]+pareto_tau_shape_k[site_pos[n],2]*treatment_pos[n]));}
}
}
