data {
  int<lower = 1> nb;
  vector[nb] nd;
  int<lower = 1> np;
  int<lower = 1> nr;
  int<lower = 1> nt;
  int<lower = 1> npred_x;
  int<lower = 1> npred_t;
  int<lower = 0> P_t;
  array[np, nb] int<lower = 0> y;
  vector[nb] otters;
  vector[npred_x] otter_pred;
  matrix[nb, P_t] X_t;
  matrix[npred_t, P_t] X_t_pred;
  matrix[nb, nr] K_r;
  matrix[nb, nt] K_t;
  matrix[nb, nb] dists;
  matrix[nb, nb] L;
}
transformed data{
  vector[nb] log_nd;
  log_nd = log(nd);
}
parameters {
  vector[np] lambda0;
  vector[np] h_raw;
  vector<lower = 0>[np] tau;
  matrix[P_t, np] gamma_raw;
  vector[P_t] mu_gamma;
  real<lower = 0> sigma_gamma;
  vector<lower = 0>[np] phi_inv;
  row_vector[np] mu_psi;
  vector<lower = 0>[np] sigma_eta;
  vector<lower = 0>[np] sigma_eps;
  matrix[nr, np] epsilon_raw;
  matrix[nt, np] eta_raw;
}
transformed parameters {
  
  matrix[nb, np] log_lambda;
  matrix[nb, np] logit_psi;
  matrix[P_t, np] gamma;
  vector[np] h;
  vector[np] phi;
  matrix[nb, np] delta;

  // Overdispersion
  phi = 1.0 ./ phi_inv;
  
  // Occupancy
  logit_psi = rep_matrix(mu_psi, nb) + K_r * diag_post_multiply(epsilon_raw, sigma_eps) + K_t * diag_post_multiply(eta_raw, sigma_eta);  

  // Counts
  gamma = rep_matrix(mu_gamma, np) + sigma_gamma * gamma_raw;
  delta = X_t * gamma;

  h = 2.0 * h_raw;
  
  for(k in 1:np){
    for(i in 1:nb){
      log_lambda[i, k] = lambda0[k] - tau[k] * (otters[i] + delta[i, k] - h[k]) ^ 2;
    }
  }
}

model {
  
  // Occupancy random effects
  mu_psi ~ normal(0, 5.0);
  sigma_eta ~ normal(0, 2.5);
  sigma_eps ~ normal(0, 2.5);
  to_vector(eta_raw) ~ std_normal();
  to_vector(epsilon_raw) ~ std_normal();

  // Maximum count
  lambda0 ~ normal(0, 2.5);
  
  // Prey response to otter abundance
  tau ~ gamma(2, 0.1);

  // Cumulative otter abundance at maximum count
  h_raw ~ std_normal();
  
  // Otter abundance shift regression coefficients
  to_vector(gamma_raw) ~ std_normal();
  mu_gamma ~ normal(0, 2.0);
  sigma_gamma ~ normal(0, 1.0);

  // Overdispersion
  phi_inv ~ normal(0, 5);
  
  // Data model 
  for(k in 1:np){
    for(i in 1:nb){
      if(y[k, i] == 0){
        target += log_sum_exp(log1m_inv_logit(logit_psi[i, k]), 
                              log_inv_logit(logit_psi[i, k]) + neg_binomial_2_log_lpmf(0 | log_nd[i] + log_lambda[i, k], phi[k]));
      } else{
        target += log_inv_logit(logit_psi[i, k]) + neg_binomial_2_log_lpmf(y[k, i] | log_nd[i] + log_lambda[i, k], phi[k]);
      }
    }
  }
}
generated quantities{
  
  array[np, nb] int yhat;
  matrix[npred_x, np] lambda_x;
  matrix[npred_t, np] lambda_t;
  matrix[npred_t, np] delta_pred;

  // Simulate observations
  for(k in 1:np){
    for(i in 1:nb){

      int zhat = bernoulli_logit_rng(logit_psi[i, k]);

      yhat[k][i] = zhat * neg_binomial_2_log_rng(log_nd[i] + log_lambda[i, k], phi[k]);
    }
  }
  
  // Calculate predicted lambdas
  delta_pred = X_t_pred * gamma;

  for(k in 1:np){
    for(i in 1:npred_x){
      lambda_x[i, k] = lambda0[k] - tau[k] * (otter_pred[i] - h[k]) ^ 2;
    }
    for(t in 1:npred_t){
      lambda_t[t, k] = lambda0[k] - tau[k] * (delta_pred[t, k] - h[k]) ^ 2;
    }
  }
}
