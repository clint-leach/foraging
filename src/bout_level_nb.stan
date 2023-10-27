data {
  int<lower = 1> nb;
  vector[nb] nd;
  int<lower = 1> np;
  int<lower = 1> nsp;
  int<lower = 1> nr;
  int<lower = 1> ny;
  int y[np, nb];
  matrix[nb, 2] X;
  int bout_idx[nb];
  int sp_idx[np];
  int y_idx[nb];
  matrix[nb, nb] dists;
  matrix[nb, nb] L;
  int blocks[2, nr];
  int block_size[nr];
}
transformed data{
  vector[nb] log_nd;
  log_nd = log(nd);
}
parameters {
  vector[np] lambda0;
  vector[np] psi0;
  vector[np] r;
  vector<lower = 0>[2] beta;
  vector<lower = 0>[np] tau;
  vector<lower = 0>[np] nu;
  vector<lower = 0>[np] phi;
  real<lower = 0> sigma;
  matrix[nb, np] omega_raw;
}
transformed parameters {
  
  // matrix[nb, np] omega;
  matrix[nb, np] logit_psi_niche;
  matrix[nb, np] logit_psi;
  matrix[nb, np] log_lambda;
  vector[nb] c;
  
  // Otter feeding niche
  c = - X * beta;
  
  // Counts by species/size
  for(k in 1:np){
    for(i in 1:nb){
      log_lambda[i, k] = lambda0[k] - nu[k] * (c[i] - r[k]) ^ 2;
      logit_psi_niche[i, k] = psi0[k] - tau[k] * (c[i] - r[k]) ^ 2;
      logit_psi[i, k] = logit_psi_niche[i, k] + sigma * omega_raw[i, k];
    }
  }
}

model {
  
  // Static occupancy regression coefficients
  beta ~ normal(0, 1.0);
  
  // Otter niche width
  tau ~ normal(0, 10);
  nu ~ normal(0, 10);
  
  // Peak counts
  lambda0 ~ normal(0, 2.5);
  
  // Peak occupancy probability
  psi0 ~ normal(0, 5);
  
  // Prey niche axis positions
  r ~ normal(0, 1);

  // Sds of bout-level random effects
  sigma ~ normal(0, 1.0);
  
  // Unscaled bout-level random effects
  for(k in 1:np){
    omega_raw[, k] ~ std_normal();
  }

  // Overdispersion
  phi ~ inv_gamma(2, 2);
  
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

  int yhat[np, nb];

  for(k in 1:np){
    for(i in 1:nb){

      int zhat = bernoulli_logit_rng(logit_psi[i, k]);

      yhat[k][i] = zhat * neg_binomial_2_log_rng(log_nd[i] + log_lambda[i, k], phi[k]);
    }
  }
}
