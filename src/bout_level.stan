data {
  int<lower = 1> nb;
  int nd[nb];
  int<lower = 1> np;
  int<lower = 1> nsp;
  int<lower = 0> nr;
  int<lower = 0> P_env;
  int<lower = 0> P_dyn;
  int y[np, nb];
  matrix[nb, P_env] X_env;
  matrix[nb, P_dyn] X_dyn;
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
  matrix[P_env, nsp] beta;
  matrix[P_dyn, np] alpha_std;
  row_vector[np] alpha0;
  matrix<lower = 0>[P_dyn, np] r1_alpha;
  matrix<lower = 0>[P_dyn, np] r2_alpha;
  matrix[P_dyn, np] delta_std;
  row_vector[np] delta0;
  matrix<lower = 0>[P_dyn, np] r1_delta;
  matrix<lower = 0>[P_dyn, np] r2_delta;  
  vector<lower = 0>[np] r1_tau;
  vector<lower = 0>[np] r2_tau;
  vector<lower = 0>[nsp] sigma;
  matrix[nb, nsp] omega_raw;
}
transformed parameters {
  
  matrix[nb, np] logit_psi;
  matrix[nb, nsp] omega;
  matrix[nb, nsp] spsp;
  matrix[nb, np] log_lambda;
  vector<lower = 0>[np] tau;
  matrix<lower = 0>[P_dyn, np] lambda_alpha;
  matrix<lower = 0>[P_dyn, np] lambda_delta;
  matrix[P_dyn, np] alpha;
  matrix[P_dyn, np] delta;

  // Computing t-distributed horseshoe prior components
  lambda_alpha = r1_alpha .* sqrt(r2_alpha);
  lambda_delta = r1_delta .* sqrt(r2_delta);
  tau = r1_tau .* sqrt(r2_tau);
  
  // Scaling coefficient matrices
  alpha = alpha_std .* lambda_alpha * diag_matrix(tau);
  delta = delta_std .* lambda_delta * diag_matrix(tau);

  // Occupancy level random effect
  for (k in 1:nsp){
    omega[, k] = sigma[k] * omega_raw[, k];
  }
  
  // Static species-level occupancy effects
  spsp = X_env * beta + omega;
  
  // Latent bout-level species/size occupancy
  logit_psi = rep_matrix(delta0, nb) + spsp * spxsz + X_dyn * delta;
  
  // Counts by species/size
  log_lambda = rep_matrix(alpha0, nb) + X_dyn * alpha;
  
}

model {
  
  // Static occupancy regression coefficients
  to_vector(beta) ~ normal(0, 2.5);
  
  // Dynamic occupancy regression coefficients
  delta0 ~ normal(0, 5);
  to_vector(delta_std) ~ std_normal();
  
  // Count regression coefficients
  alpha0 ~ normal(0, 5);
  to_vector(alpha_std) ~ std_normal();
  
  // Half-t priors on lambdas
  to_vector(r1_alpha) ~ std_normal();
  to_vector(r1_delta) ~ std_normal();
  to_vector(r2_alpha) ~ inv_gamma(0.5 * nu_lambda, 0.5 * nu_lambda);
  to_vector(r2_delta) ~ inv_gamma(0.5 * nu_lambda, 0.5 * nu_lambda);
  
  // Half-t prior on tau
  r1_tau ~ normal(0, tau0);
  r2_tau ~ inv_gamma(0.5 * nu_tau, 0.5 * nu_tau);

  // Sds of bout-level random effects
  sigma ~ normal(0, 1.0);
  
  // Unscaled bout-level random effects
  for(k in 1:nsp){
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
