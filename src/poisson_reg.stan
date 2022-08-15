data {
  int<lower=0> nb;
  int nd[nb];
  int<lower = 0> np;
  int<lower = 0> ny;
  int<lower = 0> nr;
  int<lower = 0> P;
  int y[sum(nd), np];
  matrix[nb, P] X;
  matrix[nb, ny] Z;
  int bout_idx[nb];
  matrix[nb, nb] dists;
  int blocks[2, nr];
  int block_size[nr];
  vector[np] rho;
}
parameters {
  matrix[P, np] beta;
  vector<lower = 0>[np - 1] a;
  vector<lower = 0>[np] h;
  vector<lower = 0>[np] sigma;
  matrix[nb, np] omega_raw;
  matrix[ny, np] xi;
}
transformed parameters {
  
  matrix[nb, np] omega;
  vector[np] a_all;
  matrix[nb, np] eta;
  matrix[nb, np] lambda;

  // Spatial-GP
  for (k in 1:np){
    for (i in 1:nr){
      matrix[block_size[i], block_size[i]] K;
      matrix[block_size[i], block_size[i]] L;
      
      for(j in 1:block_size[i]){
        for(l in 1:block_size[i]){
          K[j, l] = sigma[k] ^ 2 * exp(-dists[blocks[1, i] - 1 + j, blocks[1, i] - 1 + l] / rho[k]);
        }
        
        K[j, j] += 1e-8;
      }
      
      L = cholesky_decompose(K);
      omega[blocks[1, i]:blocks[2, i], k] = L * omega_raw[blocks[1, i]:blocks[2, i], k];
    } 
  }
  
  // Otter attack rate (with one value pinned to 1)
  a_all[1] = 1.0;
  a_all[2:np] = a;
  
  // Latent prey density
  eta = exp(X * beta + Z * xi + omega);
  
  // Functional response
  for(k in 1:np){
    for(i in 1:nb){
      lambda[i, k] = a_all[k] * eta[i, k] / (1 + sum(a_all .* h .* eta[i, ]'));
    }
  }
}

model {
  
  // Priors on prey species params
  
  // Attack rate
  a ~ gamma(1.0, 1.0);
  
  // Handling time
  h ~ gamma(1.0, 1.0);
  
  // Spatial marginal variance
  sigma ~ normal(0, 2);
  
  // Regression coefficients
  to_vector(beta) ~ normal(0, 2.0);
  
  // Temporal random effects
  to_vector(xi) ~ normal(0, 2.0);
  
  // Data model 
  for(k in 1:np){
    for(i in 1:nb){
      for(j in 1:nd[i]){
        y[bout_idx[i] + j, k] ~ poisson(lambda[i, k]);
      }
    }
  }
}
generated quantities{
  
  int yhat[sum(nd), np];
  
  for(k in 1:np){
    for(i in 1:nb){
      for(j in 1:nd[i]){
        yhat[bout_idx[i] + j, k] = poisson_rng(lambda[i, k]);
      }
    }
  }
}

