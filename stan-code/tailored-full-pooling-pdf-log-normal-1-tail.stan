



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

parameters {

real mu[1];
real tau[1];
real sigma_control[1];
real sigma_TE[1];
matrix[M-1,P] beta_raw; // the parent parameters minus the Mth category

}
transformed parameters{
matrix[M,P] beta_full;

  beta_full<- append_row(beta_raw,rep_row_vector(0, P));

}
model {
for (m in 1:(M-1)){
    beta_raw[m] ~ normal(0,5);
}
for (n in 1:N)
cat[n] ~ categorical_logit(beta_full * x[n]);

mu ~ normal(0,100);
tau ~ normal(0,100);
sigma_control ~ normal(0,100);
sigma_TE ~ normal(0,100);


for(n in 1:N_pos){
    y_pos[n] ~ lognormal(mu[1] + tau[1]*treatment_pos[n], exp(sigma_control[1] + sigma_TE[1]*treatment_pos[n]));
}
}
