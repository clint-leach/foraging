data {
  int<lower = 1> nb;
  int nd[nb];
  int<lower = 1> np;
  int<lower = 1> nsp;
  int<lower = 0> nr;
  int<lower = 0> P_local;
  int<lower = 0> P_global;
  int y[np, nb];
  matrix[nb, P_local] X_local;
  matrix[nb, P_global] X_global;
  int bout_idx[nb];
  int sp_idx[np];
  matrix[nsp, np] spxsz;
  matrix[nb, nb] dists;
  matrix[nb, nb] L;
  int blocks[2, nr];
  int block_size[nr];
  real<lower = 1> nu_lambda;
  real<lower = 1> nu_tau;
  real<lower = 0> tau0;
}
transformed data{
  real log_nd[nb];
  log_nd = log(nd);
}
parameters {
  // psi-level local effects
  matrix[P_local, np] beta_std;
  row_vector[np] beta0;
  vector<lower = 0>[np] r1_beta;
  vector<lower = 0>[np] r2_beta;
  
  // psi-level global effects
  matrix[P_global, np] alpha_std;
  vector<lower = 0>[np] r1_alpha;
  vector<lower = 0>[np] r2_alpha;
  
  // lambda-level local effects
  matrix[P_global, np] delta_std;
  row_vector[np] delta0;
  vector<lower = 0>[np] r1_delta;
  vector<lower = 0>[np] r2_delta;
  
  // lambda-level global effects
  matrix[P_global, np] gamma_std;
  vector<lower = 0>[np] r1_gamma;
  vector<lower = 0>[np] r2_gamma;
  
  // psi-level coefficient shrinkage
  vector<lower = 0>[np] r1_tau_occ;
  vector<lower = 0>[np] r2_tau_occ;
  
  //lambda-level coefficient shrinkage
  vector<lower = 0>[np] r1_tau_count;
  vector<lower = 0>[np] r2_tau_count;
  
  vector<lower = 0>[np] sigma;
  matrix[nb, np] omega_raw;
}
transformed parameters {
  
  matrix[nb, np] logit_psi;
  matrix[nb, np] omega;
  matrix[nb, np] log_lambda;
  vector<lower = 0>[np] tau_occ;
  vector<lower = 0>[np] tau_count;
  vector<lower = 0>[np] lambda_alpha;
  vector<lower = 0>[np] lambda_beta;
  vector<lower = 0>[np] lambda_delta;
  vector<lower = 0>[np] lambda_gamma;
  matrix[P_global, np] alpha;
  matrix[P_local, np] beta;
  matrix[P_global, np] gamma;
  matrix[P_local, np] delta;

  // Computing t-distributed horseshoe prior components
  lambda_alpha = r1_alpha .* sqrt(r2_alpha);
  lambda_beta = r1_beta .* sqrt(r2_beta);
  lambda_delta = r1_delta .* sqrt(r2_delta);
  lambda_gamma = r1_gamma .* sqrt(r2_gamma);
  tau_count = r1_tau_count .* sqrt(r2_tau_count);
  tau_occ = r1_tau_occ .* sqrt(r2_tau_occ);

  // Scaling coefficient matrices
  alpha = alpha_std * diag_matrix(lambda_alpha .* tau_occ);
  beta = beta_std * diag_matrix(lambda_beta .* tau_occ);
  
  delta = delta_std * diag_matrix(lambda_delta .* tau_count);
  gamma = gamma_std * diag_matrix(lambda_gamma .* tau_count);

  // Occupancy level random effect
  for (k in 1:np){
    omega[, k] = sigma[k] * omega_raw[, k];
  }
  
  // Latent bout-level species/size occupancy
  logit_psi = rep_matrix(beta0, nb) + X_local * beta + X_global * alpha + omega;
  
  // Counts by species/size
  log_lambda = rep_matrix(delta0, nb) + X_local * delta + X_global * gamma;
  
}

model {
  
  // Occupancy regression coefficients
  beta0 ~ normal(0, 5);
  to_vector(beta_std) ~ std_normal();
  to_vector(alpha_std) ~ std_normal();
  
  // Count regression coefficients
  delta0 ~ normal(0, 5);
  to_vector(delta_std) ~ std_normal();
  to_vector(gamma_std) ~ std_normal();
  
  // Half-t priors on lambdas
  r1_alpha ~ std_normal();
  r1_beta ~ std_normal();
  r1_delta ~ std_normal();
  r1_gamma ~ std_normal();
  r2_alpha ~ inv_gamma(0.5 * nu_lambda, 0.5 * nu_lambda);
  r2_beta ~ inv_gamma(0.5 * nu_lambda, 0.5 * nu_lambda);
  r2_delta ~ inv_gamma(0.5 * nu_lambda, 0.5 * nu_lambda);
  r2_gamma ~ inv_gamma(0.5 * nu_lambda, 0.5 * nu_lambda);

  // Half-t prior on tau
  r1_tau_occ ~ normal(0, tau0);
  r2_tau_occ ~ inv_gamma(0.5 * nu_tau, 0.5 * nu_tau);
  r1_tau_count ~ normal(0, tau0);
  r2_tau_count ~ inv_gamma(0.5 * nu_tau, 0.5 * nu_tau);
  
  // Sds of bout-level random effects
  sigma ~ normal(0, 1.0);
  
  // Unscaled bout-level random effects
  for(k in 1:np){
    omega_raw[, k] ~ std_normal();
  }
  
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
  
  // int yhat[np, sum(nd)];
  // 
  // for(k in 1:np){
  //   for(i in 1:nb){
  // 
  //     int zhat = bernoulli_logit_rng(logit_psi[i, k]);
  // 
  //     for(j in 1:nd[i]){
  //       yhat[k][bout_idx[i] + j - 1] = zhat * bernoulli_logit_rng(logit_phi[i, k]) * poisson_log_rng(log_lambda[i, k]);
  //     }
  //   }
  // }
}
