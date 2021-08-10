



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

vector[K] mu_k;
vector[K] tau_k;
vector[K] sigma_control_k;
vector[K] sigma_TE_k;
matrix[M-1,P] beta_k_raw[K]; // the hierarchical increments 

}
transformed parameters{
matrix[M,P] beta_k[K];
  for (k in 1:K){
  beta_k[k] <- append_row(beta_k_raw[k],rep_row_vector(0, P));
  }
}
model {
for (m in 1:(M-1)){
  for (k in 1:K){
    beta_k_raw[k,m] ~ normal(0,5);
  }}
for (n in 1:N)
cat[n] ~ categorical_logit(beta_k[site[n]] * x[n]);


mu_k ~ normal(0,100);
tau_k ~ normal(0,100);
sigma_control_k ~ normal(0,10);
sigma_TE_k ~ normal(0,10);


for(n in 1:N_pos){
    y_pos[n] ~ lognormal(mu_k[site_pos[n]] + tau_k[site_pos[n]]*treatment_pos[n], exp(sigma_control_k[site_pos[n]] + sigma_TE_k[site_pos[n]]*treatment_pos[n])) ;
}
}
