


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


real tau_shape;
real control_shape;
matrix[M-1,P] beta_raw; // the hierarchical increments 

}
transformed parameters{
matrix[M,P] beta_full;

  beta_full <- append_row(beta_raw,rep_row_vector(0, P));

}
model {

for (m in 1:(M-1)){
    beta_raw[m] ~ normal(0,5);
  }
for (n in 1:N)
cat[n] ~ categorical_logit(beta_full * x[n]);

tau_shape ~ normal(0,5);
control_shape ~ normal(2,5);


for(n in 1:N_pos){
    y_pos[n] ~ pareto(y_min_pos,exp(control_shape+tau_shape*treatment_pos[n]) );
}
}

