data {
int M; // number of categories
int N; // number of observations
int P; // dimentionality of beta parameter
int K; // number of sites 
int cat[N]; // category indicator
int N_pos; // number of obvs in positive tail
vector[P] x[N]; // covariates
int site[N]; // site indicator 
int treatment_pos[N_pos]; // treatment in pos tail
int site_pos[N_pos]; // site indicator in pos tail
real y_pos[N_pos]; // the positive tail 
}
transformed data{
}
parameters {

real<lower=0> control_mean;
real tau_mean;
vector[K] tau_mean_k;
vector[K] control_mean_k;
real<lower=0> tau_sigma; 
real<lower=0> control_sigma;
real<lower=0> control_sd;
real tau_sd;
vector[K] tau_sd_k;
vector[K] control_sd_k;
real<lower=0> tau_sd_sigma; 
real<lower=0> control_sd_sigma;
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
to_vector(sigma) ~ cauchy(0,2.5);
for (m in 1:M){
  for (k in 1:K){
    beta_k_raw[k,m] ~ normal(0,1);
  }}
for (n in 1:N)
cat[n] ~ categorical_logit(beta_k[site[n]] * x[n]);

tau_sd ~ normal(0,5);
tau_sd_sigma ~ cauchy(0,2.5);
control_sd ~ normal(2,5);
control_sd_sigma ~ cauchy(0,2.5);
for (k in 1:K){
  tau_sd_k[k] ~ normal(tau_sd,tau_sd_sigma);
  control_sd_k[k] ~ normal(control_sd, control_sd_sigma);
  tau_mean_k[k] ~ normal(tau_mean, tau_sigma);
  control_mean_k[k] ~ normal(control_mean, control_sigma);
}

for(n in 1:N_pos){
    y_pos[n] ~ normal(control_mean_k[site_pos[n]] + tau_mean_k[site_pos[n]]*treatment_pos[n], exp(control_sd_k[site_pos[n]]+tau_sd_k[site_pos[n]]*treatment_pos[n]) );
}
}

generated quantities{

  real posterior_predicted_tau_sd_k;
  real posterior_predicted_control_sd_k;
  real posterior_predicted_tau_mean_k;
  real posterior_predicted_control_mean_k;
  real posterior_predicted_cat;
  real posterior_cat;
  real posterior_predicted_y;
  real posterior_y;
  real treatment_sim;
  vector[P] x_sim;
  matrix[M,P] posterior_predicted_beta_k;

  treatment_sim = bernoulli_rng(0.5);
  x_sim[1] = 1;
  x_sim[2] = treatment_sim;
  
  for (m in 1:M){
    for (p in 1:P){
    posterior_predicted_beta_k[m,p] = normal_rng(beta_full[m,p],sigma[m,p]);
  }}
  posterior_cat = categorical_rng(softmax(beta_full * x_sim));
  posterior_predicted_cat = categorical_rng(softmax(posterior_predicted_beta_k * x_sim));

  if (posterior_cat == 1) posterior_y = 0;
  if (posterior_cat == 2) posterior_y = normal_rng(control_mean + tau_mean*treatment_sim, exp(control_sd+tau_sd*treatment_sim) );

  posterior_predicted_tau_sd_k= normal_rng(tau_sd,tau_sd_sigma); 
  posterior_predicted_control_sd_k = normal_rng(control_sd, control_sd_sigma);
  posterior_predicted_tau_mean_k= normal_rng(tau_mean,tau_sigma); 
  posterior_predicted_control_mean_k = normal_rng(control_mean, control_sigma);
  if (posterior_predicted_cat == 1) posterior_predicted_y = 0;
  if (posterior_predicted_cat == 2) posterior_predicted_y = normal_rng(posterior_predicted_control_mean_k + posterior_predicted_tau_mean_k*treatment_sim, exp(posterior_predicted_control_sd_k+posterior_predicted_tau_sd_k*treatment_sim) );
}



