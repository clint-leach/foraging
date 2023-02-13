data {
  int<lower = 1> nb;
  int nd[nb];
  int<lower = 1> np;
  int<lower = 1> nsp;
  int<lower = 0> nr;
  int<lower = 0> P_psi;
  int<lower = 0> P_lambda;
  int y[np, sum(nd)];
  matrix[nb, P_psi] X_psi;
  matrix[nb, P_lambda] X_lambda;
  int bout_idx[nb];
  int sp_idx[np];
  matrix[nb, nb] dists;
  matrix[nb, nb] L;
  int blocks[2, nr];
  int block_size[nr];
  vector[nsp] rho;
}
parameters {
  matrix[P_psi, nsp] beta;
  matrix[P_lambda, np] alpha;
  vector<lower = 0>[nsp] sigma;
  matrix[nb, nsp] omega_raw;
}
transformed parameters {
  
  matrix[nb, nsp] omega;
  matrix[nb, np] log_lambda;
  matrix[nb, nsp] logit_psi;

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
  
  // Latent species occupancy
  logit_psi = X_psi * beta + omega;
  
  // Counts by species/size
  log_lambda = X_lambda * alpha;
}

model {
  
  // Occupancy regression coefficients
  to_vector(beta) ~ normal(0, 2.5);
  
  // Count regression coefficients
  to_vector(alpha) ~ normal(0, 2.5);
  
  // Sds of bout-level random effects
  sigma ~ normal(0, 1.0);
  
  // Unscaled bout-level random effects
  for(k in 1:nsp){
    omega_raw[, k] ~ std_normal();
  }
  
  // Data model 
  for(k in 1:np){
    for(i in 1:nb){
      for(j in 1:nd[i]){
        if(y[k][bout_idx[i] + j - 1] == 0){
          target += log_sum_exp(log1m_inv_logit(logit_psi[i, sp_idx[k]]), 
                                log_inv_logit(logit_psi[i, sp_idx[k]]) + poisson_log_lpmf(0 | log_lambda[i, k]));
          
        } else {
          target += log_inv_logit(logit_psi[i, sp_idx[k]]) + poisson_log_lpmf(y[k][bout_idx[i] + j - 1] | log_lambda[i, k]);
        }
      }
    }
  }
}
generated quantities{
  
  int yhat[np, sum(nd)];
  int zhat[nsp, nb];

  for(k in 1:nsp){
    for(i in 1:nb){
      zhat[k][i] = bernoulli_logit_rng(logit_psi[i, sp_idx[k]]);
    }
  }
  
  for(k in 1:np){
    for(i in 1:nb){
      for(j in 1:nd[i]){
        yhat[k][bout_idx[i] + j - 1] = zhat[sp_idx[k]][i] * poisson_log_rng(log_lambda[i, k]);
      }
    }
  }
}
