library(rstan)
library(ggplot2)
library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(gtools)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load in foraging data
all <- readRDS("data/processed.rds")

# Prey counts data
data <- list(y = all$y[, -c(1:2)])

# Bout dimensions
bouts <- all$bout_length

# Covariates
X <- all$bout_dat %>% 
  dplyr::select(lat, density, depth, slope, current) %>% 
  colwise(scale)() %>% 
  as.matrix()

# Spatial covariance indexing
blocks = all$domains %>% 
  dplyr::select(first, last) %>% 
  as.matrix() %>% 
  t()

# Defining constants
const <- list(nb = nrow(bouts), 
              nd = bouts$nd, 
              np = ncol(data$y), 
              ny = ncol(all$K),
              nr = ncol(blocks),
              P = ncol(X) + 1, 
              y = all$y[, -c(1:2)],
              X = cbind(1, X), 
              Z = all$K * 1,
              bout_idx = bouts$idx,
              dists = all$d,
              blocks = blocks,
              block_size = all$domains$nbouts,
              rho = rep(100, ncol(data$y)))


inits <- list(beta = matrix(0, nrow = const$P, ncol = const$np), 
              a = rep(1.0, const$np - 1),
              h = rep(1.0, const$np),
              sigma = rep(1.0, const$np),
              omega_raw = matrix(0, nrow = const$nb, ncol = const$np),
              xi = matrix(0, nrow = const$ny, ncol = const$np))


mcmc <- stan(file = "src/poisson_reg.stan",
             data = const,
             init = list(inits),
             iter = 2000,
             chains = 1)

pairs(mcmc, pars = c("sigma", "a", "h"), inclue = T)

# Inference ====================================================================

# betas
beta <- rstan::extract(mcmc, "beta")[[1]] %>% 
  melt(varnames = c("iter", "coefficient", "prey"), value.name = "beta") %>% 
  mutate(prey = names(data$y)[prey], coefficient = colnames(const$X)[coefficient])

beta %>% 
  ggplot(aes(iter, beta)) + 
  geom_line() + 
  facet_grid(coefficient ~ prey, scales = "free")

beta %>% 
  ggplot(aes(beta)) + 
  geom_histogram(bins = 50) + 
  facet_grid(coefficient ~ prey) + 
  geom_vline(xintercept = 0) +
  theme_classic()

# Attack rates
a <- rstan::extract(mcmc, "a_all")[[1]] %>% 
  melt(varnames = c("iter", "prey"), value.name = "a") %>% 
  mutate(prey = colnames(data$y)[prey])

a %>% 
  ggplot(aes(iter, a)) + 
  geom_line() + 
  facet_wrap(~prey, scales = "free")

a %>% 
  ggplot(aes(a)) + 
  geom_histogram() + 
  geom_vline(xintercept = 1) +
  facet_wrap(~prey, scales = "free") + 
  theme_classic()

beta <- join(beta, a)

beta %>% 
  subset(coefficient == "") %>% 
  ggplot(aes(beta, a)) + 
  geom_point() + 
  facet_wrap(~prey, scales = "free")


# Handling times
h <- rstan::extract(mcmc, "h")[[1]] %>% 
  melt(varnames = c("iter", "prey"), value.name = "h") %>% 
  mutate(prey = colnames(data$y)[prey])

h %>% 
  ggplot(aes(iter, h)) + 
  geom_line() + 
  facet_wrap(~prey, scales = "free")

h %>% 
  ggplot(aes(h)) + 
  geom_histogram() + 
  facet_wrap(~prey, scales = "free") + 
  theme_classic()

# Sigmas
sigma <- rstan::extract(mcmc, "sigma")[[1]] %>% 
  melt(varnames = c("iter", "prey"), value.name = "sigma") %>% 
  mutate(prey = colnames(data$y)[prey])

sigma %>% 
  ggplot(aes(iter, sigma)) + 
  geom_line() + 
  facet_wrap(~prey, scales = "free")

sigma %>% 
  ggplot(aes(sigma)) + 
  geom_histogram() + 
  facet_wrap(~prey, scales = "free") + 
  theme_classic()

# Xi
xi <- rstan::extract(mcmc, "xi")[[1]] %>% 
  melt(varnames = c("iter", "year", "prey"), value.name = "xi") %>% 
  mutate(prey = colnames(data$y)[prey])

xi %>% 
  ggplot(aes(iter, xi)) + 
  geom_line() + 
  facet_grid(year ~ prey)

xi %>% 
  ddply(.(year, prey), summarise,
        med = median(xi),
        low = quantile(xi, 0.1),
        high = quantile(xi, 0.9)) %>% 
  ggplot(aes(year, med)) + 
  geom_point() + 
  geom_linerange(aes(ymin = low, ymax = high)) +
  facet_grid(prey ~.) + 
  theme_classic() + 
  geom_hline(yintercept = 0)

# lambda
lambda <- rstan::extract(mcmc, "lambda")[[1]] %>% 
  melt(varnames = c("iter", "bout", "prey"), value.name = "lambda") %>%
  mutate(bout_id = bouts$bout_id[bout], prey = colnames(data$y)[prey]) %>% 
  join(all$bout_dat)

lambda_sum <- lambda %>% 
  ddply(.(bout_id, prey), summarise,
        med = median(lambda), 
        low = quantile(lambda, 0.05), 
        high = quantile(lambda, 0.95),
        density = density[1], 
        year = year[1]) 

lambda_sum %>% 
  # subset(year == 2000) %>%
  ggplot(aes(bout_id, med)) +
  # geom_point() + 
  geom_linerange(aes(ymin = low, ymax = high), alpha = 0.7, size = 1.0) +
  geom_point(aes(bout_id, y), data = subset(y_long, year < 2020), color = "blue", size = 3.0, shape = "-") +
  facet_grid(prey ~ year, scales = "free") + 
  theme_classic()+
  xlab("bout") + 
  ylab(expression(lambda)) + 
  theme(axis.text.x = element_blank())
