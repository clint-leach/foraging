data {
  int<lower = 1> nb;
  int nd[nb];
  int<lower = 1> np;
  int<lower = 1> nsp;
  int<lower = 0> nr;
  int<lower = 0> P_env;
  int<lower = 0> P_dyn;
  int y[np, sum(nd)];
  int y_occ[np, nb];
  matrix[nb, P_env] X_env;
  matrix[nb, P_dyn] X_dyn;
  int bout_idx[nb];
  int sp_idx[np];
  matrix[nsp, np] spxsz;
  matrix[nb, nb] dists;
  matrix[nb, nb] L;
  int blocks[2, nr];
  int block_size[nr];
}
parameters {
  matrix[P_env, nsp] beta;
  matrix[P_dyn, np] alpha;
  row_vector[np] alpha0;
  matrix[P_dyn, np] gamma;
  row_vector[np] gamma0;
  matrix[P_dyn, np] delta;
  vector<lower = 0>[nsp] sigma;
  matrix[nb, nsp] omega_raw;
}
transformed parameters {
  
  matrix[nb, np] ll = rep_matrix(0.0, nb, np);
  matrix[nb, np] logit_psi;
  
  {
    matrix[nb, nsp] omega;
    matrix[nb, nsp] spsp;
    matrix[nb, np] log_lambda;
    matrix[nb, np] logit_phi;

    // Spatial-GP
    // for (k in 1:nsp){
    //   for (i in 1:nr){
    //     matrix[block_size[i], block_size[i]] K;
    //     matrix[block_size[i], block_size[i]] L;
    // 
    //     for(j in 1:block_size[i]){
    //       for(l in 1:block_size[i]){
    //         K[j, l] = sigma[k] ^ 2 * exp(-dists[blocks[1, i] - 1 + j, blocks[1, i] - 1 + l] / rho[k]);
    //       }
    // 
    //       K[j, j] += 1e-8;
    //     }
    // 
    //     L = cholesky_decompose(K);
    //     omega[blocks[1, i]:blocks[2, i], k] = L * omega_raw[blocks[1, i]:blocks[2, i], k];
    //   }
    // }
  
    for (k in 1:nsp){
      omega[, k] = sigma[k] * omega_raw[, k];
    }
  
    // Static species-level occupancy effects
    spsp = X_env * beta + omega;
  
    // Latent bout-level species occupancy
    logit_psi = spsp * spxsz + X_dyn * delta;
  
    // Dive-level occupancy by species/size
    logit_phi = rep_matrix(gamma0, nb) + X_dyn * gamma;
  
    // Counts by species/size
    log_lambda = rep_matrix(alpha0, nb) + X_dyn * alpha;
  
    // log-pmf of data model
    for (k in 1:np){
      for (i in 1:nb){
        for (j in 1:nd[i]){
          if(y[k][bout_idx[i] + j - 1] == 0){
            ll[i, k] += log_sum_exp(log1m_inv_logit(logit_phi[i, k]), 
                                    poisson_log_lpmf(0 | log_lambda[i, k]) + log_inv_logit(logit_phi[i, k]));
          } else{
            ll[i, k] += poisson_log_lpmf(y[k][bout_idx[i] + j - 1] | log_lambda[i, k]) + log_inv_logit(logit_phi[i, k]);
          }
        }
      }
    }
  }
}

model {
  
  // Static occupancy regression coefficients
  to_vector(beta) ~ normal(0, 1.25);
  
  // Dynamic occupancy regression coefficients
  to_vector(delta) ~ normal(0, 1.25);
  
  // Dive occupancy coefficients
  gamma0 ~ normal(0, 2.5);
  to_vector(gamma) ~ normal(0, 1.25);
  
  // Count regression coefficients
  alpha0 ~ normal(0, 2.5);
  to_vector(alpha) ~ normal(0, 1.25);
  
  // Sds of bout-level random effects
  sigma ~ normal(0, 1.0);
  
  // Unscaled bout-level random effects
  for(k in 1:nsp){
    omega_raw[, k] ~ std_normal();
  }
  
  // Data model 
  for(k in 1:np){
    for(i in 1:nb){
      if(y_occ[k, i] == 0){
        target += log_sum_exp(log1m_inv_logit(logit_psi[i, k]), 
                              log_inv_logit(logit_psi[i, k]) + ll[i, k]);
      } else{
        target += log_inv_logit(logit_psi[i, k]) + ll[i, k];
      }
    }
  }
}
generated quantities{
  
  int yhat[np, sum(nd)];

  for(k in 1:np){
    for(i in 1:nb){

      int zhat = bernoulli_logit_rng(logit_psi[i, k]);

      for(j in 1:nd[i]){
        yhat[k][bout_idx[i] + j - 1] = zhat * bernoulli_logit_rng(logit_phi[i, k]) * poisson_log_rng(log_lambda[i, k]);
      }
    }
  }
}
