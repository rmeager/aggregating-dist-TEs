



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


matrix[K,2] tau_shape_k;
matrix[K,2] control_shape_k;
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

for(m in 1:2){
tau_shape_k[m] ~ normal(0,5);
control_shape_k[m] ~ normal(2,5);
}
for (n in 1:N_neg){
    y_neg[n] ~ pareto(y_min_neg,exp(control_shape_k[site_neg[n],1]+tau_shape_k[site_neg[n],1]*treatment_neg[n]) );
}
for(n in 1:N_pos){
    y_pos[n] ~ pareto(y_min_pos,exp(control_shape_k[site_pos[n],2]+tau_shape_k[site_pos[n],2]*treatment_pos[n]) );
}
}
