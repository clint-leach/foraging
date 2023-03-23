library(rstan)
library(ggplot2)
library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(gtools)
library(stringr)
library(patchwork)

# Load in foraging data
all <- readRDS("data/processed.rds")

# Long format data
y_long <- all$y %>% 
  # select(!c(clam_4, crab_4, modiolus_4, mussel_3, urchin_3)) %>% 
  select(c(bout_id, dive_num, clam_1, clam_2, clam_3, urchin_1, urchin_2)) %>% 
  pivot_longer(clam_1:urchin_2,
               names_to = "prey", 
               values_to = "y") %>% 
  mutate(prey_sp = str_split(prey, pattern = "_", simplify = TRUE)[, 1], 
         prey_sz = str_split(prey, pattern = "_", simplify = TRUE)[, 2]) %>% 
  join(all$bout_dat)

# Species lookup
prey_sp <- ddply(y_long, .(prey), summarise,
                 prey_sp = prey_sp[1],
                 prey_sz = prey_sz[1]) %>% 
  mutate(prey_idx = c(1:length(prey)))
  
sp_idx <- prey_sp$prey_sp %>% as.factor() %>% as.numeric()

spxsz <- matrix(0, nrow = max(sp_idx), ncol = length(sp_idx))
for(i in 1:nrow(spxsz)){
  spxsz[i, ] <- sp_idx == i
}

# Read in MCMC output
mcmc <- readRDS("output/poisson_horse.rds")

# Plotting data ================================================================
 
y_long %>% 
  ggplot(aes(density, y)) + 
  geom_point() + 
  facet_grid(prey ~ ., scales = "free_y")

y_long %>% 
  # subset(density > 150) %>% 
  # subset(site %in% c("Inner Beardslee", "Outer Beardslee")) %>% 
  ggplot(aes(density, y)) + 
  geom_point() + 
  geom_hline(yintercept = 0) +
  facet_grid(prey ~ site, scales = "free_y")

occ <- y_long %>% 
  ddply(.(bout_id, prey), summarise,
        z = any(y > 0),
        p = mean(y > 0),
        success_mean = mean(y[y > 0]),
        n = length(y),
        bout_mean = mean(y)) %>% 
  join(all$bout_dat)

occ %>% 
  ggplot(aes(density, z)) + 
  geom_point() + 
  facet_grid(prey ~ .)

occ %>% 
  subset(prey == "urchin_2") %>% 
  # subset(z > 0) %>% 
  ggplot(aes(density, p, alpha = n)) + 
  geom_point() + 
  facet_grid(prey ~ .)

occ %>% 
  subset(z > 0) %>% 
  ggplot(aes(density, success_mean)) + 
  geom_point() + 
  facet_grid(prey ~ ., scales = "free_y")

occ %>% 
  # subset(z > 0) %>% 
  ggplot(aes(density, bout_mean)) + 
  geom_point() + 
  facet_grid(prey ~ ., scales = "free_y")
  
all$bout_dat %>% 
  ggplot(aes(year, density, color = site)) + 
  geom_point()

# Plotting main effects of otters ==============================================

# Covariate matrix
X <- all$bout_dat %>% 
  dplyr::select(lat, depth, slope, current, density) %>% 
  colwise(function(x){(x - mean(x)) / (2 * sd(x))})() %>%
  mutate(density_sq = density ^ 2) %>%
  as.matrix()

X_env <- X[, 1:4]
X_dyn <- X[, 5:6]

# Coefficients
beta <- rstan::extract(mcmc, "beta", permute = T)[[1]]
delta <- rstan::extract(mcmc, "delta", permute = T)[[1]]
gamma <- rstan::extract(mcmc, "gamma", permute = T)[[1]]
alpha <- rstan::extract(mcmc, "alpha", permute = T)[[1]]

gamma0 <- rstan::extract(mcmc, "gamma0", permute = T)[[1]]
alpha0 <- rstan::extract(mcmc, "alpha0", permute = T)[[1]]
delta0 <- rstan::extract(mcmc, "delta0", permute = T)[[1]]

lambda_alpha <- rstan::extract()
tau <- rstan::extract(mcmc, "tau", permute = T)[[1]]
plot(tau[, 5], type = "l")

# lambda (counts given capture)
log_lambda <- aaply(alpha, c(1), function(a) X_dyn %*% a) %>% 
  aaply(c(2), function(x) x + alpha0) %>%
  melt(varnames = c("bout_idx", "iter", "prey_idx"), value.name = "log_lambda") %>%
  mutate(lambda = exp(log_lambda)) %>% 
  ddply(.(bout_idx, prey_idx), summarise,
        lambda_median = median(lambda),
        lambda_low = quantile(lambda, 0.1),
        lambda_high = quantile(lambda, 0.9)) %>% 
  join(all$bout_dat) %>% 
  join(prey_sp)

log_lambda %>% 
  ggplot(aes(density, lambda_median, color = prey_sz)) + 
  geom_ribbon(aes(density, ymin = lambda_low, ymax = lambda_high, fill = prey_sz), alpha = 0.5, inherit.aes = FALSE) +
  geom_line() + 
  scale_color_discrete(guide = "none") +
  scale_fill_discrete(guide = "none") +
  facet_grid(. ~ prey_sp) + 
  theme_classic() + 
  xlab("otter density") +
  ylab(expression(log(lambda))) -> figa

# phi (dive probability given bout presence)
logit_phi <- aaply(alpha, c(1), function(a) X_dyn %*% a) %>% 
  aaply(c(2), function(x) x + alpha0) %>%
  melt(varnames = c("bout_idx", "iter", "prey_idx"), value.name = "logit_phi") %>%
  # melt(varnames = c("iter", "bout_idx", "prey_idx"), value.name = "logit_phi") %>% 
  mutate(phi = inv.logit(logit_phi)) %>% 
  ddply(.(bout_idx, prey_idx), summarise,
        phi_median = median(phi),
        phi_low = quantile(phi, 0.1),
        phi_high = quantile(phi, 0.9)) %>% 
  join(all$bout_dat) %>% 
  join(prey_sp)

logit_phi %>% 
  # subset(prey_sp == "urchin") %>%
  ggplot(aes(density, phi_median, color = prey_sz)) + 
  geom_ribbon(aes(density, ymin = phi_low, ymax = phi_high, fill = prey_sz), alpha = 0.5, inherit.aes = FALSE) +
  geom_line() + 
  # geom_point(aes(density, p, alpha = n), data = subset(occ, prey == "urchin_2"), inherit.aes = FALSE) +
  scale_color_discrete(guide = "none") +
  scale_fill_discrete(guide = "none") +
  facet_grid(. ~ prey_sp) +
  theme_classic() + 
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  xlab("otter density") +
  ylab(expression(phi)) -> figb

# psi (bout-level occupancy probability)
logit_psi <- aaply(delta, c(1), function(a) X_dyn %*% a) %>% 
  aaply(c(2), function(x) x + delta0) %>%
  melt(varnames = c("bout_idx", "iter", "prey_idx"), value.name = "logit_psi") %>%
  # melt(varnames = c("iter", "bout_idx", "prey_idx"), value.name = "logit_psi") %>%
  mutate(psi = inv.logit(logit_psi)) %>% 
  ddply(.(bout_idx, prey_idx), summarise,
        psi_median = median(psi),
        psi_low = quantile(psi, 0.1),
        psi_high = quantile(psi, 0.9)) %>% 
  join(all$bout_dat) %>% 
  join(prey_sp)

logit_psi %>% 
  # subset(prey_sp == "urchin") %>% 
  ggplot(aes(density, psi_median, color = prey_sz)) + 
  geom_ribbon(aes(density, ymin = psi_low, ymax = psi_high, fill = prey_sz), alpha = 0.5, inherit.aes = FALSE) +
  geom_line() +
  facet_grid(. ~ prey_sp, scales = "free_y") + 
  theme_classic() + 
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  xlab("otter density") +
  ylab(expression(psi)) -> figc

figa + figc + plot_layout(ncol = 1)

# Hierarchical plots ===========================================================

lambda <- aaply(alpha, c(1), function(a) X_dyn %*% a) %>% 
  aaply(c(2), function(x) exp(x + alpha0)) 

psi <- aaply(delta, c(1), function(a) X_dyn %*% a) %>% 
  aaply(c(2), function(x) inv.logit(x + delta0))

z <- rbinom(length(psi), size = 1, psi) %>% 
  array(dim = dim(psi))

lambdaz <- lambda * z

foo <- lambdaz %>% 
  melt(varnames = c("bout_idx", "iter", "prey_idx"), value.name = "lambda") %>%
  subset(lambda > 0) %>%
  ddply(.(bout_idx, prey_idx), summarise,
        lambda_var = var(lambda),
        lambda_median = median(lambda),
        lambda_mean = mean(lambda),
        lambda_low = quantile(lambda, 0.1),
        lambda_high = quantile(lambda, 0.9)) %>% 
  join(all$bout_dat) %>% 
  join(prey_sp)



joint %>% 
  # subset(prey_sp == "urchin") %>% 
  ggplot(aes(density, phi_median, color = psi_median)) + 
  geom_line(size = 1) + 
  geom_line(aes(density, phi_low, color = psi_median)) +
  geom_line(aes(density, phi_high, color = psi_median)) +
  scale_color_distiller(palette = "Greys", direction = 1) +
  facet_grid(prey_sz ~ prey_sp) +
  theme_classic()

foo %>% 
  # subset(psi_high > 0.01) %>%
  ggplot(aes(density, lambda_median, color = prey_sz)) + 
  geom_ribbon(aes(density, ymin = lambda_low, ymax = lambda_high, fill = prey_sz), alpha = 0.5, inherit.aes = FALSE) +
  geom_line() + 
  scale_color_brewer(palette = "Dark2", guide = "none") +
  scale_fill_brewer(palette = "Dark2", guide = "none") +
  facet_grid(. ~ prey_sp) +
  theme_classic() + 
  theme(strip.background = element_blank()) +
  xlab("otter density") +
  ylab(expression(lambda))

foo %>% 
  subset(prey == "urchin_2") %>% 
  ggplot(aes(density, lambda_var, color = prey_sz)) + 
  geom_line() +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  facet_grid(. ~ prey_sp) +
  theme_classic() + 
  theme(strip.background = element_blank()) +
  xlab("otter density") +
  ylab(expression(lambda))


joint %>% 
  subset(psi_high > 0.1) %>%
  ggplot(aes(density, phi_median, color = prey_sz)) + 
  geom_ribbon(aes(density, ymin = phi_low, ymax = phi_high, fill = prey_sz), alpha = 0.5, inherit.aes = FALSE) +
  geom_line() + 
  scale_color_brewer(palette = "Dark2", guide = "none") +
  scale_fill_brewer(palette = "Dark2", guide = "none") +
  facet_grid(. ~ prey_sp) +
  theme_classic() + 
  theme(strip.background = element_blank()) +
  xlab("otter density") +
  ylab(expression(phi)) -> figa
 

joint %>% 
  # subset(prey_sp == "urchin") %>% 
  ggplot(aes(density, psi_median, color = prey_sz)) + 
  geom_ribbon(aes(density, ymin = psi_low, ymax = psi_high, fill = prey_sz), alpha = 0.5, inherit.aes = FALSE) +
  geom_line() +
  scale_color_brewer("prey size", palette = "Dark2") +
  scale_fill_brewer("prey size", palette = "Dark2") +
  facet_grid(. ~ prey_sp, scales = "free_y") + 
  theme_classic() + 
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  xlab("otter density") +
  ylab(expression(psi)) -> figb

figa + figb + plot_layout(ncol = 1)
