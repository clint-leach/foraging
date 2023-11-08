data {
  int<lower = 1> nb;
  int nd[nb];
  int<lower = 1> np;
  int<lower = 0> P_local;
  int<lower = 0> P_global;
  int<lower = 0> P_env;
  int y[np, nb];
  matrix[nb, P_local] X_local;
  matrix[nb, P_global] X_global;
  matrix[nb, P_env] X_env;
  matrix[nb, nb] dists;
  matrix[nb, nb] L;
  real<lower = 1> nu_lambda;
  real<lower = 1> nu_tau;
  real<lower = 0> tau0;
}
transformed data{
  real log_nd[nb];
  log_nd = log(nd);
}
parameters {
  
  // psi-level environmental effects
  row_vector[np] alpha0;
  matrix[P_env, np] alpha;

  // lambda-level local effects
  matrix[P_local, np] delta_std;
  row_vector[np] delta0;
  vector<lower = 0>[np] r1_delta;
  vector<lower = 0>[np] r2_delta;
  
  // lambda-level global effects
  matrix[P_global, np] gamma_std;
  vector<lower = 0>[np] r1_gamma;
  vector<lower = 0>[np] r2_gamma;
  
  //lambda-level coefficient shrinkage
  vector<lower = 0>[np] r1_tau_count;
  vector<lower = 0>[np] r2_tau_count;
  
  // real<lower = 0> sigma;
  // matrix[nb, np] omega_raw;
}
transformed parameters {
  
  matrix[nb, np] logit_psi;
  matrix[nb, np] log_lambda;
  vector<lower = 0>[np] tau_count;
  vector<lower = 0>[np] lambda_delta;
  vector<lower = 0>[np] lambda_gamma;
  matrix[P_global, np] gamma;
  matrix[P_local, np] delta;

  // Computing t-distributed horseshoe prior components
  lambda_delta = r1_delta .* sqrt(r2_delta);
  lambda_gamma = r1_gamma .* sqrt(r2_gamma);
  tau_count = r1_tau_count .* sqrt(r2_tau_count);

  // Scaling coefficient matrices
  delta = delta_std * diag_matrix(lambda_delta .* tau_count);
  gamma = gamma_std * diag_matrix(lambda_gamma .* tau_count);
  
  // Latent bout-level species/size occupancy
  logit_psi = rep_matrix(alpha0, nb) + X_env * alpha; #+ sigma * L * omega_raw;
  
  // Counts by species/size
  log_lambda = rep_matrix(delta0, nb) + X_local * delta + X_global * gamma;
  
}

model {
  
  // Occupancy regression coefficients
  alpha0 ~ normal(0, 5);
  to_vector(alpha) ~ normal(0, 2.5);
  
  // Count regression coefficients
  delta0 ~ normal(0, 2.5);
  to_vector(delta_std) ~ std_normal();
  to_vector(gamma_std) ~ std_normal();
  
  // Half-t priors on lambdas
  r1_delta ~ std_normal();
  r1_gamma ~ std_normal();
  r2_delta ~ inv_gamma(0.5 * nu_lambda, 0.5 * nu_lambda);
  r2_gamma ~ inv_gamma(0.5 * nu_lambda, 0.5 * nu_lambda);

  // Half-t prior on tau
  r1_tau_count ~ normal(0, tau0);
  r2_tau_count ~ inv_gamma(0.5 * nu_tau, 0.5 * nu_tau);
  
  // Sd of bout-level spatial random effects
  // sigma ~ normal(0, 1.0);
  
  // Unscaled bout-level random effects
  // for(k in 1:np){
  //   omega_raw[, k] ~ std_normal();
  // }
  
  // Data model 
  for(k in 1:np){
    for(i in 1:nb){
      if(y[k, i] == 0){
        target += log_sum_exp(log1m_inv_logit(logit_psi[i, k]), 
                              log_inv_logit(logit_psi[i, k]) + poisson_log_lpmf(0 | log_nd[i] + log_lambda[i, k]));
      } else{
        target += log_inv_logit(logit_psi[i, k]) + poisson_log_lpmf(y[k, i] | log_nd[i] + log_lambda[i, k]);
      }
    }
  }
}
generated quantities{
  
  array[np, nb] int yhat;

  for(k in 1:np){
    for(i in 1:nb){

      int zhat = bernoulli_logit_rng(logit_psi[i, k]);

      yhat[k][i] = zhat * poisson_log_rng(log_nd[i] + log_lambda[i, k]);
    }
  }
}
