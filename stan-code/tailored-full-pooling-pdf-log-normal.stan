



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

real mu[2];
real tau[2];
real sigma_control[2];
real sigma_TE[2];
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

for(m in 1:2){
mu[m] ~ normal(0,100);
tau[m] ~ normal(0,100);
sigma_control[m] ~ normal(0,100);
sigma_TE[m] ~ normal(0,100);
}
for (n in 1:N_neg){
    y_neg[n] ~ lognormal(mu[1] + tau[1]*treatment_neg[n], exp(sigma_control[1] + sigma_TE[1]*treatment_neg[n]));
}
for(n in 1:N_pos){
    y_pos[n] ~ lognormal(mu[2] + tau[2]*treatment_pos[n], exp(sigma_control[2] + sigma_TE[2]*treatment_pos[n]));
}
}
