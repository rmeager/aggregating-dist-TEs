



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
transformed data{
real y_min_neg;
real y_min_pos;
y_min_neg <- min(y_neg);
y_min_pos <- min(y_pos);
}
parameters {

real control_shape[2];
real tau_shape[2];
matrix[K,2] tau_shape_k;
matrix[K,2] control_shape_k;
real<lower=0> tau_sigma[2]; 
real<lower=0> control_sigma[2];
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

tau_shape ~ normal(0,5);
tau_sigma ~ cauchy(0,2);
control_shape ~ normal(2,5);
control_sigma ~ cauchy(0,2);
for (k in 1:K){
  for (i in 1:2){
  tau_shape_k[k,i] ~ normal(tau_shape[i],tau_sigma[i]);
  control_shape_k[k,i] ~ normal(control_shape[i], control_sigma[i]);
  }
}
for (n in 1:N_neg){
    y_neg[n] ~ pareto(y_min_neg,exp(control_shape_k[site_neg[n],1]+tau_shape_k[site_neg[n],1]*treatment_neg[n]) );
}
for(n in 1:N_pos){
    y_pos[n] ~ pareto(y_min_pos,exp(control_shape_k[site_pos[n],2]+tau_shape_k[site_pos[n],2]*treatment_pos[n]) );
}
}
