library(rstan)
library(ggplot2)
library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(gtools)
library(stringr)
library(splines2)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load in foraging data
all <- readRDS("data/processed.rds")

bouts <- all$bout_dat

# Covariates
otters <- bouts$cumulative / sd(bouts$cumulative)
otters_pred <- seq(0, max(otters), length.out = 100)
t_occ <- bouts$occ_time - 1993

# Basis function matrices for time of occupancy
X_t <- bSpline(t_occ, knots = seq(min(t_occ) + 5, max(t_occ) - 5, length.out = 4), degree = 3)
X_t_pred <- bSpline(seq(min(t_occ), max(t_occ), by = 1), knots = seq(min(t_occ) + 5, max(t_occ) - 5, length.out = 4), degree = 3, intercept = F)

# Cholesky matrix for spatial covariance
rho <- 100 ^ 2
K <- exp(-all$d ^ 2/ rho) + diag(1e-8, nrow = dim(bouts)[1])
L = chol(K)

# Defining constants
const <- list(nb = nrow(bouts), 
              nd = bouts$nd, 
              np = ncol(all$y) - 1,
              nr = ncol(all$K_r),
              nt = ncol(all$K_t),
              npred_x = length(otters_pred),
              npred_t = nrow(X_t_pred),
              P_t = ncol(X_t),
              y = t(all$y[, -1]),
              otters = otters,
              otter_pred = otters_pred,
              X_t = X_t,
              X_t_pred = X_t_pred,
              dists = all$d,
              L = L,
              K_r = all$K_r,
              K_t = all$K_t)

inits <- list(
  list(lambda0 = rep(0.0, const$np),
       h_raw = rep(0.0, const$np),
       tau = rep(0.1, const$np),
       gamma_raw = matrix(0.0, const$P_t, const$np),
       mu_gamma = rep(0.0, const$P_t),
       sigma_gamma = 0.1,
       phi_inv = rep(1.0, const$np),
       mu_psi = rep(0.0, const$np),
       sigma_eta = rep(1.0, const$np),
       sigma_eps = rep(1.0, const$np),
       eta_raw = matrix(0, nrow = const$nt, ncol = const$np),
       epsilon_raw = matrix(0, nrow = const$nr, ncol = const$np)),
  list(lambda0 = rnorm(const$np),
       h_raw = rnorm(const$np),
       tau = rep(1.0, const$np),
       gamma_raw = matrix(0.0, const$P_t, const$np),
       mu_gamma = rnorm(const$P_t),
       sigma_gamma = 0.5,
       phi_inv = rep(0.1, const$np),
       mu_psi = rnorm(const$np),
       sigma_eta = rep(0.1, const$np),
       sigma_eps = rep(0.1, const$np),
       eta_raw = matrix(0, nrow = const$nt, ncol = const$np),
       epsilon_raw = matrix(0, nrow = const$nr, ncol = const$np)))

mcmc <- stan(file = "src/bout_level_nb.stan",
             data = const,
             init = inits,
             iter = 4000,
             chains = 2, 
             pars = c("logit_psi", "delta"), 
             include = FALSE)


saveRDS(mcmc, "output/bout_level_nb_chain.rds")
