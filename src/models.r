zip <- nimbleCode({
  
  # Priors on bout-level environmental occupancy coefficients
  for(k in 1:nsp){
    for(l in 1:P_psi){
      beta[l, k] ~ dnorm(0, sd = 2.5)
    }
  }
  
  # Priors on bout-level sea otter occupancy coefficients
  for(k in 1:np){
    for(l in 1:(P_lambda - 1)){
      delta[l, k] ~ dnorm(0, sd = 2.5)
    }
  }
  
  # Priors on dive-level occupancy coefficients
  for(k in 1:np){
    for(l in 1:P_lambda){
      gamma[l, k] ~ dnorm(0, sd = 2.5)
    }
  }
  
  # Priors on count coefficients
  for(k in 1:np){
    for(l in 1:P_lambda){
      alpha[l, k] ~ dnorm(0, sd = 5.0)
    }
  }
  
  # Priors on random effect standard deviation
  for(k in 1:nsp){
    sigma[k] ~ T(dnorm(0, sd = 2.5), 0, Inf)
  }
  
  # Bout-level Occupancy random effect
  for(k in 1:nsp){
    for(i in 1:nb){
      omega[i, k] ~ dnorm(0, sd = sigma[k])
    }
  }  
  
  # Species-level environmental effects on bout-level occupancy
  logit_psi_env[1:nb, 1:nsp] <- X_psi[1:nb, 1:P_psi] %*% beta[1:P_psi, 1:nsp] + omega[1:nb, 1:nsp]
  
  # Adding size-level otter effects on bout-level occupancy
  logit(psi[1:nb, 1:np]) <- logit_psi_env[1:nb, 1:nsp] %*% spxsz[1:nsp, 1:np] + X_lambda[1:nb, 2:P_lambda] %*% delta[1:(P_lambda - 1), 1:np]
  
  # Bout-level counts
  log(lambda[1:nb, 1:np]) <- X_lambda[1:nb, 1:P_lambda] %*% alpha[1:P_lambda, 1:np]
  
  # Dive-level occupancy
  logit(phi[1:nb, 1:np]) <- X_lambda[1:nb, 1:P_lambda] %*% gamma[1:P_lambda, 1:np]

  # Likelihood
  for(k in 1:np){
    for(i in 1:nb){
      
      z[i, k] ~ dbern(psi[i, k])
      zhat[i, k] ~ dbern(psi[i, k])
      
      y[i, k] ~ dpois(nd[i] * lambda[i, k] * z[i, k])
      yhat[i, k] ~ dpois(nd[i] * lambda[i, k] * zhat[i, k])
      
      for(j in 1:nd[i]){
        w[bout_idx[i] + j - 1, k] ~ dbern(z[i, k] * phi[i, k])
        y[bout_idx[i] + j - 1, k] ~ dpois(z[i, k] * lambda[i, k])

        what[bout_idx[i] + j - 1, k] ~ dbern(zhat[i, k] * phi[i, k])
        yhat[bout_idx[i] + j - 1, k] ~ dpois(z[i, k] * lambda[i, k])
      }
    }
  }
})

expcov <- nimbleFunction(     
  run = function(dists = double(2), rho = double(0), sigma = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    sigma2 <- sigma*sigma
    for(i in 1:n)
      for(j in 1:n)
        result[i, j] <- sigma2*exp(-dists[i,j]/rho) + 1e-8 * (i == j)
    return(result)
  })

block_chol <- nimbleFunction(     
  run = function(dists = double(2), blocks = integer(2), rho = double(0), sigma = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    nr <- dim(blocks)[2]
    
    Sigma <- matrix(nrow = n, ncol = n, init = FALSE)
    L <- matrix(value = 0.0, nrow = n, ncol = n)
    
    sigma2 <- sigma*sigma
    for(i in 1:nr){
      for(j in blocks[1, i]:blocks[2, i]){
        for(k in blocks[1, i]:blocks[2, i]){
         Sigma[j, k] <- sigma2 * exp(-dists[j, k] / rho) + 1e-8 * (k == j)
        }
      }
      
      L[blocks[1, i]:blocks[2, i], blocks[1, i]:blocks[2, i]] = chol(Sigma[blocks[1, i]:blocks[2, i], blocks[1, i]:blocks[2, i]])
    }
    
    return(L)
  })

