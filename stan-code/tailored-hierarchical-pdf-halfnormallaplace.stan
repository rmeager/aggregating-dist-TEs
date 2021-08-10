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
vector[Q] quantile_vec ; // vector of quantiles of interest
}

parameters {
real log_alpha[2]; // will be entered within an exponent into the distribution
real nu[2];
real log_tau[2]; //will be entered within an exponent into the distribution
real<lower=0> sd_log_alpha[2];
real<lower=0> sd_nu[2];
real<lower=0> sd_log_tau[2];
matrix[K,2] log_alpha_k;
matrix[K,2] nu_k;
matrix[K,2] log_tau_k;
real log_alpha_TE[2]; // will be entered within an exponent into the distribution
real nu_TE[2];
real log_tau_TE[2]; //will be entered within an exponent into the distribution
real<lower=0> sd_log_alpha_TE[2];
real<lower=0> sd_nu_TE[2];
real<lower=0> sd_log_tau_TE[2];
matrix[K,2] log_alpha_TE_k;
matrix[K,2] nu_TE_k;
matrix[K,2] log_tau_TE_k;
matrix[M-1,P] beta; // the parent parameters minus the Mth category
matrix<lower=0>[M,P] sigma; // the set of M*P parent variances (not a covariance matrix)
matrix[M,P] beta_k_raw[K]; // the hierarchical increments 
real<lower=0> epsilon_neg[N_neg];
real<lower=0> epsilon_pos[N_pos];
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
vector[N_neg] alpha_slot_y_neg;
vector[N_neg] nu_slot_y_neg;
vector[N_neg] tau_slot_y_neg;
vector[N_pos] alpha_slot_y_pos;
vector[N_pos] nu_slot_y_pos;
vector[N_pos] tau_slot_y_pos;
real eps_pos_over_alpha[N_pos];
real eps_neg_over_alpha[N_neg];

to_vector(beta) ~ normal(0,5);
to_vector(sigma) ~ cauchy(0,2);
for (m in 1:M){
  for (k in 1:K){
    beta_k_raw[k,m] ~ normal(0,1);
  }}
for (n in 1:N)
cat[n] ~ categorical_logit(beta_k[site[n]] * x[n]);

// priors
log_alpha  ~ normal(3,2); # we must do this for identification reasons, see here https://discourse.mc-stan.org/t/double-pareto-lognormal-distribution-in-stan/10097/20
nu ~ normal(0,5);
log_tau   ~ normal(0,5);
sd_log_alpha ~ cauchy(0,2);
sd_nu ~ cauchy(0,2);
sd_log_tau ~ cauchy(0,2);
log_alpha_TE  ~ normal(0,5);
nu_TE ~ normal(0,5);
log_tau_TE   ~ normal(0,5);
sd_log_alpha_TE ~ cauchy(0,2);
sd_nu_TE ~ cauchy(0,2);
sd_log_tau_TE ~ cauchy(0,2);
for (k in 1:K){
  for (i in 1:2){
  log_alpha_k[k,i] ~ normal(log_alpha[i],sd_log_alpha[i]);
  log_tau_k[k,i] ~ normal(log_tau[i], sd_log_tau[i]);
  nu_k[k,i] ~ normal(nu[i],sd_nu[i]);
  log_alpha_TE_k[k,i] ~ normal(log_alpha_TE[i],sd_log_alpha_TE[i]);
  log_tau_TE_k[k,i] ~ normal(log_tau_TE[i], sd_log_tau_TE[i]);
  nu_TE_k[k,i] ~ normal(nu_TE[i],sd_nu_TE[i]);
  }
}
epsilon_neg ~ exponential(1);
epsilon_pos ~ exponential(1);

for (n in 1:N_neg){
  alpha_slot_y_neg[n] = exp(log_alpha_k[site_neg[n],1] + log_alpha_TE_k[site_neg[n],1]*treatment_neg[n]);
  eps_neg_over_alpha[n] = epsilon_neg[n]*inv(alpha_slot_y_neg[n]);
  nu_slot_y_neg[n]  = (nu_k[site_neg[n],1] + nu_TE_k[site_neg[n],1]*treatment_neg[n]);
  tau_slot_y_neg[n]  = exp(log_tau_k[site_neg[n],1] + log_tau_TE_k[site_neg[n],1]*treatment_neg[n]);
      y_neg[n] ~  normal(nu_slot_y_neg[n] + eps_neg_over_alpha[n], tau_slot_y_neg[n]);

}
for(n in 1:N_pos){
  alpha_slot_y_pos[n] = exp(log_alpha_k[site_pos[n],1] + log_alpha_TE_k[site_pos[n],1]*treatment_pos[n]);
  eps_pos_over_alpha[n] = epsilon_pos[n]*inv(alpha_slot_y_pos[n]);
  nu_slot_y_pos[n]  = (nu_k[site_pos[n],1] + nu_TE_k[site_pos[n],1]*treatment_pos[n]);
  tau_slot_y_pos[n]  = exp(log_tau_k[site_pos[n],1] + log_tau_TE_k[site_pos[n],1]*treatment_pos[n]);
      y_pos[n] ~ normal(nu_slot_y_pos[n] + eps_pos_over_alpha[n], tau_slot_y_pos[n]);
}


}
