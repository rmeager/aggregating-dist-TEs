functions {
real pareto_quantile(real u, real dist_min, real shape) {
real quantile_value;
quantile_value = pow((shape*pow(dist_min, shape)/u), (1/shape+1));
return quantile_value;
}
}


data {
int M; // number of categories
int N; // number of observations
int P; // dimentionality of beta parameter
int K; // number of sites 
int Q; // number of quantiles
int cat[N]; // category indicator
int N_pos; // number of obvs in positive tail
vector[P] x[N]; // covariates
int site[N]; // site indicator 
int treatment_pos[N_pos]; // treatment in pos tail
int site_pos[N_pos]; // site indicator in pos tail
real y_pos[N_pos]; // the positive tail 
vector[Q] quantile_vec ; // vector of quantiles of interest
}
transformed data{
real y_min_pos;
y_min_pos <- min(y_pos);
}
parameters {

real control_shape;
real tau_shape;
vector[K] tau_shape_k;
vector[K] control_shape_k;
real<lower=0> tau_sigma; 
real<lower=0> control_sigma;
matrix[M-1,P] beta; // the parent parameters minus the Mth category
matrix<lower=0>[M,P] sigma; // the set of M*P parent variances (not a covariance matrix)
matrix[M,P] beta_k_raw[K]; // the hierarchical increments 

}
transformed parameters{
matrix[M,P] beta_full;
matrix[M,P] beta_k[K];
real treatment_shape;
vector[K] treatment_shape_k;

beta_full <- append_row(beta,rep_row_vector(0, P));
for (m in 1:M){
 for (k in 1:K){
  for (p in 1:P){
   beta_k[k,m,p] <- beta_full[m,p] + sigma[m,p]*beta_k_raw[k,m,p];
}}}


treatment_shape <- tau_shape + control_shape;
treatment_shape_k <- tau_shape_k + control_shape_k;


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
  tau_shape_k[k] ~ normal(tau_shape,tau_sigma);
  control_shape_k[k] ~ normal(control_shape, control_sigma);
}

for(n in 1:N_pos){
    y_pos[n] ~ pareto(y_min_pos,exp(control_shape_k[site_pos[n]]+tau_shape_k[site_pos[n]]*treatment_pos[n]) );
}
}





