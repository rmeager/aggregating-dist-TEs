



data {
int M; // number of categories
int N; // number of observations
int P; // dimensionality of beta parameter
int K; // number of sites 
int Q; // number of quantiles
int cat[N]; // category indicator
int N_pos; // number of obvs in positive tail
vector[P] x[N]; // covariates
int site[N]; // site indicator 
int treatment_pos[N_pos]; // treatment in pos tail
int site_pos[N_pos]; // site indicator in pos tail
real y_pos[N_pos]; // the positive tail 
vector[Q] quantile_vec ; // vector of quantiles of interest, this is currently not used 
}

parameters {

real mu[1];
real tau[1];
real<lower=0> sd_mu[1];
real<lower=0> sd_tau[1];
real sigma_control[1];
real sigma_TE[1];
real<lower=0> sd_sigma_control[1];
real<lower=0> sd_sigma_TE[1];
vector[K] mu_k;
vector[K] tau_k;
vector[K] sigma_control_k;
vector[K] sigma_TE_k;
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
for (k in 1:K){
  mu_k[k] ~ normal(mu,sd_mu);
  tau_k[k] ~ normal(tau, sd_tau);
  sigma_control_k[k] ~ normal(sigma_control,sd_sigma_control);
  sigma_TE_k[k] ~ normal(sigma_TE, sd_sigma_TE);

}
for(n in 1:N_pos){
    y_pos[n] ~ lognormal(mu_k[site_pos[n]] + tau_k[site_pos[n]]*treatment_pos[n], exp(sigma_control_k[site_pos[n]] + sigma_TE_k[site_pos[n]]*treatment_pos[n])) ;
}
}
