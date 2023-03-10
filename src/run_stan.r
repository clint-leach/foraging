library(rstan)
library(ggplot2)
library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(gtools)
library(stringr)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load in foraging data
all <- readRDS("data/processed.rds")

# Bout dimensions
bouts <- all$bout_length

# Observations
y <- all$y[, -c(1:2)]

y_occ <- all$y %>% 
  pivot_longer(!(bout_id:dive_num),
               names_to = "prey", 
               values_to = "y") %>% 
  ddply(.(prey, bout_id), summarise,
        occ = sum(y) > 0) %>% 
  pivot_wider(names_from = bout_id, values_from = occ)

# Covariates
X <- all$bout_dat %>% 
  dplyr::select(lat, depth, slope, current, density) %>% 
  colwise(scale)() %>%
  mutate(density_sq = density ^ 2) %>%
  as.matrix()

X_env <- X[, 1:4]
X_dyn <- X[, 5:6]

# Spatial covariance indexing
blocks = all$domains %>% 
  dplyr::select(first, last) %>% 
  as.matrix() %>% 
  t()

# Species id of each prey type
sp_idx <- all$sp_idx

spxsz <- matrix(0, nrow = max(sp_idx), ncol = length(sp_idx))
for(i in 1:nrow(spxsz)){
  spxsz[i, ] <- sp_idx == i
}

# Cholesky matrix for spatial covariance
rho <- 10 ^ 2
K <- exp(-all$d ^ 2/ rho) + diag(1e-6, nrow = dim(all$bout_length)[1])
L <- K * 0

for(i in 1:ncol(blocks)){
  L[blocks[1, i]:blocks[2, i], blocks[1, i]:blocks[2, i]] <- chol(K[blocks[1, i]:blocks[2, i], blocks[1, i]:blocks[2, i]])
}

# Defining constants
const <- list(nb = nrow(bouts), 
              nd = bouts$nd, 
              np = ncol(y),
              nsp = max(sp_idx),
              nr = ncol(blocks),
              P_env = ncol(X_env) + 1, 
              P_dyn = ncol(X_dyn),
              y = y %>% t(),
              y_occ = y_occ[, -1],
              X_env = cbind(1, X_env), 
              X_dyn = X_dyn,
              sp_idx = sp_idx,
              spxsz = spxsz,
              bout_idx = bouts$idx,
              dists = all$d,
              L = L,
              blocks = blocks,
              block_size = all$domains$nbouts,
              rho = rep(10, max(sp_idx))
              )

inits <- list(beta = matrix(0, nrow = const$P_env, ncol = const$nsp),
              delta = matrix(0, nrow =const$P_dyn, ncol = const$np),
              alpha = matrix(0, nrow = const$P_dyn, ncol = const$np),
              alpha0 = rep(0, const$np),
              gamma = matrix(0, nrow = const$P_dyn, ncol = const$np),
              gamma0 = rep(0, const$np),
              sigma = rep(0.1, const$nsp),
              omega = matrix(0, nrow = const$nb, ncol = const$nsp))


mcmc <- stan(file = "src/poisson_reg.stan",
             data = const,
             init = list(inits),
             iter = 1000,
             chains = 1, 
             control = list(adapt_delta = 0.65))

saveRDS(mcmc, "output/multilevel.rds")

# Inference ====================================================================

# betas
beta <- rstan::extract(mcmc, "beta", permuted = TRUE)[[1]] %>% 
  melt(varnames = c("iter", "coefficient", "prey"), value.name = "beta") %>% 
  mutate(coefficient = colnames(const$X_psi)[coefficient])

beta %>% 
  ggplot(aes(iter, beta)) + 
  geom_line() + 
  facet_grid(coefficient ~ prey, scales = "free")

beta %>% 
  # subset(coefficient == "") %>%
  ggplot(aes(beta)) + 
  geom_histogram(bins = 50) + 
  facet_grid(prey ~ coefficient, scales = "free_x") + 
  geom_vline(xintercept = 0) +
  theme_classic()

# alphas
alpha <- rstan::extract(mcmc, "alpha", permuted = TRUE)[[1]] %>% 
  melt(varnames = c("iter", "coefficient", "prey"), value.name = "alpha") %>% 
  mutate(prey = names(y)[prey], coefficient = colnames(const$X_dyn)[coefficient])

alpha %>% 
  ggplot(aes(iter, alpha)) + 
  geom_line() + 
  facet_grid(coefficient ~ prey, scales = "free")

alpha %>% 
  # subset(coefficient == "") %>%
  ggplot(aes(alpha)) + 
  geom_histogram(bins = 50) + 
  facet_grid(prey ~ coefficient, scales = "free_x") + 
  geom_vline(xintercept = 0) +
  theme_classic()

# Sigmas
sigma <- rstan::extract(mcmc, "sigma")[[1]] %>% 
  melt(varnames = c("iter", "prey"), value.name = "sigma")

sigma %>% 
  ggplot(aes(iter, sigma)) + 
  geom_line() + 
  facet_wrap(~prey, scales = "free")

sigma %>% 
  ggplot(aes(sigma)) + 
  geom_histogram() + 
  facet_grid(prey ~., scales = "free") + 
  theme_classic() + 
  xlab("standard deviation of bout random effect")

# Goodness of fit
y_long <- all$y %>% 
  pivot_longer(clam_1:urchin_3,
               names_to = "prey", 
               values_to = "y") %>% 
  join(all$bout_length) %>% 
  join(all$bout_dat)

y_long %>% 
  ddply(.(prey), summarise,
        total = sum(y))

y_sum <- y_long %>% 
  ddply(.(bout_id, prey), summarise,
        bout_mean = mean(y),
        prop_zero = sum(y == 0) / nd[1]) %>% 
  join(all$bout_dat)

y_sum %>% 
  # subset(prey == "clam_1") %>% 
  ggplot(aes(year, bout_mean)) + 
  geom_point() + 
  facet_grid(prey ~ site, scales = "free_y") +  
  theme_classic() +
  theme(axis.text.x = element_blank()) 
    

# yhat
yhat_sum <- rstan::extract(mcmc, "yhat")[[1]] %>% 
  melt(varnames = c("iter", "prey", "dive"), value.name = "yhat") %>%
  mutate(prey = colnames(y)[prey], 
         bout_id = all$y$bout_id[dive], 
         dive_num = all$y$dive_num[dive]) %>% 
  ddply(.(bout_id, dive_num, prey), summarise,
        med = median(yhat),
        mean = mean(yhat),
        low = quantile(yhat, 0.05), 
        high = quantile(yhat, 0.95),
        pzero = mean(yhat == 0)) %>% 
  join(y_long)

# Dynamics at a particular site
yhat_sum %>% 
  subset(site == "Boulder") %>%
  ggplot(aes(bout_id, mean)) +
  geom_point(size = 1.0) +
  geom_linerange(aes(ymin = low, ymax = high), alpha = 0.25, size = 0.5) +
  geom_point(aes(bout_id, y), color = "blue", size = 2.0, shape = "-") +
  facet_grid(prey ~ year, scales = "free") +
  theme_classic()+
  xlab("bout") + 
  ylab("y") + 
  theme(axis.text.x = element_blank()) + 
  ggtitle("Boulder")

# Dynamics against otter density
yhat_sum %>% 
  ggplot(aes(density, med)) +
  # geom_point(size = 1.0) +
  geom_linerange(aes(ymin = low, ymax = high), alpha = 0.7, size = 0.5) +
  geom_point(aes(density, bout_mean), color = "blue", size = 1.0) +
  facet_grid(prey ~ ., scales = "fixed") + 
  theme_classic()+
  xlab("otter density")

# counts per dive
yhat_sum %>% 
  ggplot(aes(mean, y, color = dive_num)) + 
  geom_point() + 
  facet_wrap(~prey, scales = "free") + 
  geom_abline(intercept = 0, slope = 1) + 
  xlab("posterior mean lambda") + 
  ylab("mean count per dive")

# mean counts per bout
yhat_sum %>% 
  ddply(.(bout_id, dive_num), summarise,
        y_mean = mean(y),
        yhat_mean = mean(mean)) %>% 
  ggplot(aes(yhat_mean, y_mean)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1)


# p(presence)
yhat_sum %>% 
  ggplot(aes(pzero, prop_zero)) + 
  geom_point() + 
  facet_wrap(~prey, scales = "free") + 
  geom_abline(intercept = 0, slope = 1)

# psi
psi <- rstan::extract(mcmc, "logit_psi")[[1]] %>% 
  melt(varnames = c("iter", "bout", "prey_sp"), value.name = "logit_psi") %>% 
  mutate(bout_id = bouts$bout_id[bout],
         psi = inv.logit(logit_psi)) %>% 
  ddply(.(bout_id, prey_sp), summarise,
        med = median(psi),
        high = quantile(psi, 0.9),
        low = quantile(psi, 0.1)) %>% 
  join(all$bout_dat)

# Dynamics at a single site
psi %>% 
  ggplot(aes(bout_id, med)) + 
  geom_point() + 
  geom_linerange(aes(ymin = low, ymax = high)) +
  facet_grid(prey_sp ~ year, scales = "free") +
  theme_classic()+
  xlab("bout") + 
  ylab("psi")

# Dynamics against otter density
psi %>% 
  ggplot(aes(density, med, color = year)) + 
  geom_point() + 
  geom_linerange(aes(ymin = low, ymax = high)) +
  facet_grid(prey_sp ~ ., scales = "fixed") +
  theme_classic()+
  xlab("otter density") + 
  ylab("psi")
  

# lambda
lambda <- rstan::extract(mcmc, "log_lambda")[[1]] %>% 
  melt(varnames = c("iter", "bout", "prey"), value.name = "log_lambda") %>%
  mutate(bout_id = bouts$bout_id[bout], prey = colnames(y)[prey]) %>% 
  mutate(lambda = exp(log_lambda))

lambda_sum <- lambda %>% 
  ddply(.(bout_id, prey), summarise,
        med = median(lambda), 
        low = quantile(lambda, 0.05), 
        high = quantile(lambda, 0.95)) %>% 
  join(all$bout_dat) %>% 
  join(y_sum)

# Lambda per bout at a given site
lambda_sum %>% 
  ggplot(aes(bout_id, med)) +
  geom_point(size = 1.0) +
  geom_linerange(aes(ymin = low, ymax = high), alpha = 0.7, size = 0.5) +
  geom_point(aes(bout_id, bout_mean), color = "blue", size = 3.0, shape = "-") +
  facet_grid(prey ~ year, scales = "free") + 
  theme_classic()+
  xlab("bout") + 
  ylab(expression(lambda))

# Prey counts as a function of otter density
lambda_sum %>% 
  ggplot(aes(density, med)) +
  geom_point(size = 1.0) +
  geom_linerange(aes(ymin = low, ymax = high), alpha = 0.5, size = 0.5) +
  geom_point(aes(density, bout_mean), color = "blue", size = 1.0) +
  facet_grid(prey ~., scales = "free") + 
  theme_classic()+
  xlab("otter density") + 
  ylab(expression(lambda))

# lambda against covariates
lambda_sum %>% 
  ggplot(aes(current, med)) +
  geom_point(size = 1.0) +
  geom_linerange(aes(ymin = low, ymax = high), alpha = 0.5, size = 0.5) +
  facet_grid(prey ~., scales = "free") + 
  theme_classic()+
  xlab("otter density") + 
  ylab(expression(lambda))

# Omegas
omega <- rstan::extract(mcmc, "omega_raw")[[1]] %>% 
  melt(varnames = c("iter", "bout", "prey"), value.name = "omega")

plot(omega$omega[omega$prey == 1], omega$omega[omega$prey == 2])

# Model checking ===============================================================

# Data summaries by dive
y_check <- y_long %>% 
  ddply(.(bout_id, dive_num), summarise,
        success = sum(y) > 0,
        nprey = sum(y > 0))

# Posterior predictive realizations
yhat <- rstan::extract(mcmc, "yhat")[[1]] %>% 
  melt(varnames = c("iter", "prey", "dive"), value.name = "yhat") %>%
  mutate(prey = colnames(y)[prey], 
         bout_id = all$y$bout_id[dive],
         dive_num = all$y$dive_num[dive]) 

# Summaries by dive
post_check <- yhat %>% 
  subset(iter %in% sample(1:2000, 200)) %>% 
  ddply(.(iter, dive), summarise,
        success = sum(yhat) > 0,
        nprey = sum(yhat > 0))

post_sum <- post_check %>% 
  ddply(.(iter), summarise,
        p_success = mean(success),
        nprey = mean(nprey[success]))

# Proportion of successful dives
post_sum %>% 
  ggplot(aes(p_success)) + 
  geom_histogram(bins = 20) + 
  geom_vline(xintercept = mean(y_check$success))

# Mean number of prey types per successful dive
post_sum %>% 
  ggplot(aes(nprey)) + 
  geom_histogram(bins = 20) + 
  geom_vline(xintercept = mean(y_check$nprey[y_check$success]))

# Posterior distribution of proportion of zeros by species
y_pzero <- y_long %>% 
  ddply(.(prey), summarise,
        pzero = mean(y == 0))

yhat %>% 
  ddply(.(iter, prey), summarise,
        pzero = mean(yhat == 0)) %>% 
  ggplot(aes(pzero)) + 
  geom_histogram() +
  geom_vline(aes(xintercept = pzero), data = y_pzero) +
  facet_grid(prey ~.)


# Posterior distribution of per-bout variance? Or overall variance by prey species?


# Mechanistic simulation =======================================================

a <- c(5.0, 1.0, 0.1)
h <- c(1.0, 1.0, 0.1)
alpha <- c(0.2, 0.1, 0.01)
gamma <- c(2, 5, 8)
N <- seq(0, 10, length.out = 99)
R <- matrix(10.0, nrow = 3, ncol = 100)
lambda <- matrix(0, nrow = 3, ncol = 99)
for(i in 1:99){
  # R[, i] <- 10 * sigmoid(gamma - N[i])
  lambda[, i] <- a * R[, i] / (1 + sum(a * h * R[, i]))
  R[, i + 1] <- R[, i] * exp(-alpha * lambda[, i] * N[i])
}

plot(R[1, ], type = "l")
lines(R[2, ], col = "blue")
lines(R[3, ], col = "red")

plot(lambda[1, ], type = "l", ylim = c(0, 1))
lines(lambda[2, ], col = "blue")
lines(lambda[3, ], col = "red")

plot(N, log(R[1, 1:99]), type = "l")
lines(N, log(R[2, 1:99]), col = "blue")
lines(N, log(R[3, 1:99]), col = "red")

plot(N, log(lambda[1, ]), type = "l")
lines(N, log(lambda[2, ]), col = "blue")
lines(N, log(lambda[3, ]), col = "red")


omega <- 1
alpha = 0.8
y <- rep(NA, 100)
y[1] <- 0
for(i in 2:100){
  y[i] = rpois(1, omega + alpha * y[i -1])
}
plot(y, type = "l")


register_google(key = "AIzaSyD1cMY6KL2aj7EiXaqfe85Iir_dMTFEPwQ")

map <- get_googlemap(center = c(lon = 144.92517277926538, lat = 13.565700972935051), zoom = 14, maptype = "satellite")

ggmap(map)
