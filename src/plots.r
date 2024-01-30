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
library(sf)
library(lubridate)

# Load in foraging data
all <- readRDS("data/processed.rds")

# Bout info
boutdat <- all$bout_dat

# Long format data (summarised to bout level)
y_long <- all$y %>% 
  pivot_longer(clam_1:urchin_2,
               names_to = "prey", 
               values_to = "y") %>% 
  join(all$bout_dat) %>% 
  mutate(prey_sp = str_split(prey, pattern = "_", simplify = TRUE)[, 1], 
         prey_sz = str_split(prey, pattern = "_", simplify = TRUE)[, 2],
         ypd = y / nd)

# Species lookup
prey_sp <- ddply(y_long, .(prey), summarise,
                 prey_sp = prey_sp[1],
                 prey_sz = as.numeric(prey_sz[1])) %>% 
  mutate(prey_idx = c(1:length(prey)),
         prey_sp = factor(prey_sp, labels = c("clam", "crab", "italic(Modiolus)", "italic(Mytilus)", "scallop", "snail", "star", "urchin")))
  
# Read in MCMC output
mcmc <- readRDS("output/bout_level_nb_chain.rds")

# Figure 1: Bout locations =====================================================

# Reading in Glacier Bay shapefile for mapping
gb <- st_read(dsn = "data/PH6502/historicl1.shp")

ggplot(gb) + 
  geom_sf(size = 0.25) + 
  geom_point(aes(long, lat, color = bout_year), data = arrange(boutdat, desc(bout_year)), size = 0.05) + 
  scale_color_viridis_c(guide = "none") +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.text = element_text(size = 8),
        panel.grid = element_blank()) +
  scale_y_continuous(expand = c(0.0, 0)) +
  xlim(c(-136.79, -135.8)) +
  xlab("longitude") + 
  ylab("latitude") + 
  ggtitle("(a)") + 
  annotation_scale(location = "bl") -> fig1a

ggplot(gb) + 
  geom_sf(size = 0.25) + 
  geom_point(aes(long, lat, color = occ_time), data = arrange(boutdat, desc(bout_year)), size = 0.05) + 
  scale_color_viridis_c("year") +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.text.y = element_blank(),
        axis.text = element_text(size = 8),
        panel.grid = element_blank()) +
  scale_y_continuous(expand = c(0.0, 0)) +
  xlim(c(-136.79, -135.8)) +
  xlab("longitude") +
  ylab("")  + 
  ggtitle("(b)") +
  annotation_north_arrow(location = "bl", height = unit(0.5, "in"), width = unit(0.4, "in")) -> fig1b

pdf(file = "output/figures/Figure1.pdf",
    width = 8.5, height = 5)

fig1a + fig1b + plot_layout(ncol = 2)

dev.off()

# Figure 2: distribution of cumulative occupancy over time =====================

all$bout_dat %>%
  ggplot(aes(bout_year, cumulative)) + 
  geom_jitter(width = 0.25, size = 0.2) + 
  theme_classic() + 
  scale_y_continuous(expand = c(0, 2)) +
  theme(axis.text = element_text(size = 8)) +
  ylab("local cumulative otter abundance") + 
  xlab("year") -> fig2

pdf(file = "output/figures/Figure2.pdf",
    width = 4, height = 4)

fig2

dev.off()

# Figure 3: Diet across otter gradient at a fixed t_occ ========================

otters <- all$bout_dat$cumulative / sd(all$bout_dat$cumulative)

otters_pred <- data.frame(otter_idx = 1:100,
                          otters = sd(boutdat$cumulative) * seq(0, max(otters), length.out = 100))

lambda_x <- rstan::extract(mcmc, "lambda_x", permute = TRUE)[[1]] %>% 
  melt(varnames = c("iter", "otter_idx", "prey_idx"), value.name = "log_lambda") %>% 
  mutate(lambda = exp(log_lambda)) %>%
  ddply(.(otter_idx, prey_idx), summarise,
        mean_lambda = mean(lambda),
        low_lambda = quantile(lambda, 0.1),
        high_lambda = quantile(lambda, 0.9)) %>% 
  join(prey_sp) %>% 
  join(otters_pred)

lambda_x %>% 
  ggplot(aes(otters, mean_lambda, color = factor(prey_sz))) +
  geom_line() + 
  geom_ribbon(aes(otters, ymin = low_lambda, ymax = high_lambda, fill = factor(prey_sz)), alpha = 0.5, inherit.aes = FALSE) +
  scale_color_brewer("prey size", palette = "Dark2") +
  scale_fill_brewer("prey size", palette = "Dark2") +
  facet_grid(prey_sp ~ ., scales = "free_y", labeller = label_parsed) + 
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(n.breaks = 3) +
  xlab("local cumulative sea otter abundance") +
  ylab(expression(lambda)) -> fig3


pdf(file = "output/figures/Figure3.pdf",
    width = 4, height = 7)

fig3

dev.off()

# Figure 4: Diet at invasion front =============================================

# Estimated lambda for x = 0 across t_occ
lambda_t <- rstan::extract(mcmc, "lambda_t", permute = TRUE)[[1]] %>% 
  melt(varnames = c("iter", "occ_year", "prey_idx"), value.name = "log_lambda") %>% 
  mutate(lambda = exp(log_lambda), occ_year = occ_year + 1992) %>%
  ddply(.(occ_year, prey_idx), summarise,
        mean_lambda = mean(lambda),
        low_lambda = quantile(lambda, 0.1),
        high_lambda = quantile(lambda, 0.9)) %>%
  join(prey_sp)

lambda_t %>%   
  ggplot(aes(occ_year, mean_lambda, color = factor(prey_sz))) +
  geom_path() +
  geom_ribbon(aes(occ_year, ymin = low_lambda, ymax = high_lambda, fill = factor(prey_sz)), alpha = 0.5, inherit.aes = FALSE) +
  scale_color_brewer("prey size", palette = "Dark2") +
  scale_fill_brewer("prey size", palette = "Dark2") +
  facet_grid(prey_sp ~ ., scales = "free_y", labeller = label_parsed) + 
  theme_classic() +
  theme(strip.background = element_blank()) +
  ylab(expression(lambda)) +
  scale_x_continuous("year", expand = c(0, 0)) +
  scale_y_continuous(n.breaks = 4) -> fig4


pdf(file = "output/figures/Figure4.pdf",
    width = 4, height = 7)

fig4

dev.off()

# Figure 5: temporal occupancy trends ==========================================

# Temporal random effects for occupancy
sigma_eta <- rstan::extract(mcmc, "sigma_eta", permute = TRUE)[[1]]
mu_psi <-  rstan::extract(mcmc, "mu_psi", permute = TRUE)[[1]]

years <- sort(unique(all$bout_dat$bout_year))

psi_t <- rstan::extract(mcmc, "eta_raw", permute = TRUE)[[1]] %>% 
  aaply(2, function(x) mu_psi + x * sigma_eta) %>% 
  melt(varnames = c("year_idx", "iter", "prey_idx"), value.name = "logit_psi") %>% 
  mutate(occ_year = years[year_idx],
         psi = inv.logit(logit_psi)) %>% 
  ddply(.(occ_year, prey_idx), summarise,
        mean_psi = mean(psi),
        low_psi = quantile(psi, 0.1),
        high_psi = quantile(psi, 0.9)) %>% 
  join(prey_sp)

psi_t %>% 
  ggplot(aes(occ_year, mean_psi, color = factor(prey_sz))) +
  geom_point(size = 0.75) + 
  geom_path(alpha = 0.5) +
  geom_linerange(aes(occ_year, ymin = low_psi, ymax = high_psi, color = factor(prey_sz)), alpha = 0.5, inherit.aes = FALSE) +
  scale_color_brewer("prey size", palette = "Dark2") +
  scale_fill_brewer("prey size", palette = "Dark2") +
  facet_grid(prey_sp ~ ., scales = "fixed", labeller = label_parsed) + 
  theme_classic() +
  theme(strip.background = element_blank()) +
  ylab(expression(psi)) +
  scale_x_continuous("year", expand = c(0, 0.2), breaks = seq(1995, 2015, by = 5)) +
  scale_y_continuous(limits = c(0.0, 1.0), n.breaks = 3) -> fig5


pdf(file = "output/figures/Figure5.pdf",
    width = 4, height = 7)

fig5

dev.off()

# Figure 6: trends in average size  ============================================

bout_ppd <- rstan::extract(mcmc, "yhat")[[1]] %>% 
  magrittr::extract(seq(1, 2000, by = 10), , ) %>%
  melt(varnames = c("iter", "prey_idx", "bout_idx"), value.name = "yhat") %>%
  join(prey_sp) %>% 
  ddply(.(bout_idx, iter), summarise,
        success = any(yhat > 0),
        richness = ifelse(success, length(unique(prey_sp[yhat > 0])), NA),
        avg_size = ifelse(success, sum(prey_sz * yhat) / sum(yhat), NA)) %>% 
  ddply(.(bout_idx), summarise,
        mean_richness = mean(richness, na.rm = T),
        mean_size = mean(avg_size, na.rm = T),
        p_success = mean(success)) %>% 
  join(boutdat)


# Mean size per bout
bout_ppd %>% 
  subset(occ_time == 1993) %>% 
  ggplot(aes(cumulative, mean_size, color = bout_year)) + 
  geom_point(size = 0.25) + 
  scale_color_viridis_c(guide = "none") +
  theme_classic() + 
  ylim(c(1, 3)) +
  xlim(c(0, max(boutdat$cumulative))) +
  ylab("posterior predicted average prey size") + 
  xlab("local cumulative \notter abundance") +
  ggtitle("(a)") -> fig6a

bout_ppd %>% 
  subset(bout_year > 2017) %>% 
  ggplot(aes(cumulative, mean_size, color = bout_year)) + 
  geom_point(size = 0.25) + 
  scale_color_viridis_c(limits = c(1993, 2019), guide = "none") +
  theme_classic() + 
  scale_y_continuous("", limits = c(1, 3), labels = NULL) +
  xlab("local cumulative \notter abundance") +
  ggtitle("(b)") -> fig6b

bout_ppd %>% 
  ggplot(aes(bout_year, mean_size, color = bout_year)) + 
  geom_point(size = 0.25) + 
  scale_color_viridis_c(guide = "none") +
  theme_classic() + 
  scale_y_continuous("", limits = c(1, 3), labels = NULL) +
  xlab("bout year") +
  ggtitle("(c)") -> fig6c

pdf(file = "output/figures/Figure6.pdf",
    width = 8, height = 4)

fig6a + fig6b + fig6c + plot_layout(nrow = 1)

dev.off()


# Snapshot of yhat in a fixed year =============================================

yhat_mean <- rstan::extract(mcmc, "yhat")[[1]] %>% 
  aaply(c(2, 3), mean) %>% 
  melt(varnames = c("prey_idx", "bout_idx"), value.name = "yhat") %>%
  join(boutdat) %>% 
  mutate(yhat = yhat / nd) %>%
  join(prey_sp)

yhat_mean %>% 
  subset(bout_year > 2017) %>%
  ggplot(aes(cumulative, yhat, color = occ_time)) + 
  geom_point() + 
  scale_color_viridis_c("time of \noccupancy") +
  facet_wrap(~ prey, scales = "fixed") + 
  theme_classic() + 
  ylab("mean posterior predicted count/dive") + 
  xlab("cumulative local otter abundance")


# Snapshot of lambda across gradient in a fixed year ===========================

lambda <- rstan::extract(mcmc, "log_lambda", permute = TRUE)[[1]] %>% 
  melt(varnames = c("iter", "bout_idx", "prey_idx"), value.name = "log_lambda") %>% 
  mutate(lambda = exp(log_lambda)) %>%
  ddply(.(bout_idx, prey_idx), summarise,
        mean_lambda = mean(lambda),
        low_lambda = quantile(lambda, 0.1),
        high_lambda = quantile(lambda, 0.9)) %>% 
  join(prey_sp) %>% 
  join(all$bout_dat)

lambda %>% 
  subset(bout_year > 2017) %>% 
  ggplot(aes(cumulative, mean_lambda, color = occ_time)) + 
  geom_point() + 
  scale_color_viridis_c() +
  facet_wrap(~ prey, scales = "free_y")


# Plotting spatial occupancy ===================================================

sigma_eps <- rstan::extract(mcmc, "sigma_eps", permute = TRUE)[[1]]
mu_psi <-  rstan::extract(mcmc, "mu_psi", permute = TRUE)[[1]]

regions <- sort(unique(all$bout_dat$site))

psi_r <- rstan::extract(mcmc, "epsilon_raw", permute = TRUE)[[1]] %>% 
  aaply(2, function(x) mu_psi + x * sigma_eps) %>% 
  melt(varnames = c("region_idx", "iter", "prey_idx"), value.name = "logit_psi") %>% 
  mutate(site = regions[region_idx],
         psi = inv.logit(logit_psi)) %>% 
  ddply(.(site, prey_idx), summarise,
        mean_psi = mean(psi),
        low_psi = quantile(psi, 0.1),
        high_psi = quantile(psi, 0.9)) %>% 
  join(prey_sp)

psi_r %>% 
  ggplot(aes(site, mean_psi, color = factor(prey_sz))) +
  geom_point(size = 0.75) + 
  geom_linerange(aes(site, ymin = low_psi, ymax = high_psi, color = factor(prey_sz)), alpha = 0.5, inherit.aes = FALSE) +
  scale_color_brewer("prey size", palette = "Dark2") +
  scale_fill_brewer("prey size", palette = "Dark2") +
  facet_grid(prey ~ ., scales = "fixed") + 
  theme_classic() +
  theme(strip.background = element_blank()) +
  ylab(expression(psi)) +
  xlab("region") + 
  theme(axis.text.x=element_text(angle=45, hjust=1),
        strip.text = element_text(size = 7.0))

# Variation in h ===============================================================

h <- rstan::extract(mcmc, "h", permute = TRUE)[[1]] %>% 
  melt(varnames = c("iter", "prey_idx"), value.name = "h0")

delta <- rstan::extract(mcmc, "delta_pred", permute = T)[[1]] %>% 
  melt(varnames = c("iter", "occ_year", "prey_idx"), value.name = "delta") %>% 
  join(h) %>% 
  mutate(h = h0 - delta, occ_year = occ_year + 1992) %>% 
  ddply(.(occ_year, prey_idx), summarise,
        mean_h = mean(h),
        low_h = quantile(h, 0.1),
        high_h = quantile(h, 0.9)) %>% 
  join(prey_sp)

delta %>% 
  ggplot(aes(occ_year, mean_h, color = prey_sz, group = prey_sz)) + 
  geom_line() +
  geom_ribbon(aes(ymin = low_h, ymax = high_h, fill = prey_sz), alpha = 0.5) + 
  scale_color_brewer("prey size", palette = "Dark2") +
  scale_fill_brewer("prey size", palette = "Dark2") +
  facet_grid(.~ prey_sp) + 
  theme_classic() + 
  ylab("h") + 
  xlab("year of first occupancy") 


# Plotting data ================================================================

y_long %>% 
  ddply(.(prey), summarise,
        total = sum(y),
        pbouts = mean(y > 0))

# Counts of all prey at a specific site
y_long %>% 
  subset(site %in% c("Boulder")) %>%
  ggplot(aes(cumulative, ypd)) + 
  geom_point() + 
  facet_grid(prey_sp ~ prey_sz, scales = "free_y")

# Counts of a specific prey (by site)
y_long %>% 
  subset(prey == "scallop_1") %>% 
  ggplot(aes(cumulative, ypd, color = bout_year)) + 
  geom_point() + 
  facet_wrap(. ~ site, scales = "free_y")

# Counts of a specific prey (by bout year)
y_long %>% 
  subset(prey == "scallop_1") %>%
  ggplot(aes(cumulative, ypd, color = occ_time)) + 
  geom_point() + 
  facet_grid(bout_year ~ ., scales = "free_y")

# Counts of a specific prey (by occupancy year)
y_long %>% 
  subset(prey == "scallop_1") %>%
  ggplot(aes(cumulative, ypd, color = bout_year)) + 
  geom_point() + 
  facet_grid(occ_time ~ ., scales = "free_y")

# Plotting y_ppd ===============================================================

yhat_sum <- rstan::extract(mcmc, "yhat")[[1]] %>% 
  magrittr::extract(seq(1, 4000, by = 5), , ) %>%
  melt(varnames = c("iter", "prey_idx", "bout_idx"), value.name = "yhat") %>%
  join(boutdat) %>% 
  mutate(yhat = yhat / nd) %>%
  ddply(.(bout_idx, prey_idx), summarise,
        med = median(yhat),
        mean = mean(yhat),
        low = quantile(yhat, 0.1), 
        high = quantile(yhat, 0.9),
        pzero = mean(yhat == 0),
        max = max(yhat)) %>% 
  join(select(prey_sp, prey, prey_idx)) %>% 
  join(y_long, by = c("prey", "bout_idx"))

# Estimated vs observed catch per dive
yhat_sum %>% 
  ggplot(aes(ypd, mean, color = bout_year)) + 
  geom_point() + 
  scale_color_viridis_c() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ prey, scales = "free")


# Estimated catch per dive vs otters at a given site
yhat_sum %>% 
  subset(site == "Boulder") %>% 
  ggplot(aes(cumulative, mean)) +
  geom_point(color = "gray50", alpha = 0.5) + 
  geom_linerange(aes(ymin = low, ymax = high), color = "gray50", alpha = 0.5) +
  geom_point(aes(cumulative, ypd), color = "black", size = 0.5) +
  facet_grid(prey_sp ~ prey_sz, scales = "free_y") + 
  theme_classic()+
  xlab("otter density")

# Estimated catch per dive vs time at a given site
yhat_sum %>% 
  subset(site == "Boulder") %>% 
  ggplot(aes(bout_year, mean)) +
  geom_point(color = "gray50", alpha = 0.5) + 
  geom_linerange(aes(ymin = low, ymax = high), color = "gray50", alpha = 0.5) +
  geom_point(aes(bout_year, ypd), color = "black", size = 0.5) +
  facet_grid(prey_sp ~ prey_sz, scales = "free_y") + 
  theme_classic()+
  xlab("year")

# Estimated catch per dive vs otters for a given prey (by size and year)
yhat_sum %>% 
  subset(prey_sp == "urchin") %>% 
  ggplot(aes(cumulative, mean)) +
  geom_point(color = "gray50", alpha = 0.5) + 
  geom_linerange(aes(ymin = low, ymax = high), color = "gray50", alpha = 0.5) +
  geom_point(aes(cumulative, ypd), color = "black", size = 0.5) +
  facet_grid(bout_year ~ prey_sz, scales = "free_y") + 
  theme_classic()+
  xlab("otter density")

# Specific prey type through time and by site
yhat_sum %>% 
  subset(prey == "modiolus_2") %>% 
  ggplot(aes(bout_year, mean)) +
  geom_point(color = "gray50", alpha = 0.5) + 
  geom_linerange(aes(ymin = low, ymax = high), color = "gray50", alpha = 0.5) +
  geom_point(aes(bout_year, ypd), color = "black", size = 0.5) +
  facet_wrap(~site, scales = "free_y") + 
  theme_classic()+
  xlab("year")

# Counts for a specific time of occupancy
yhat_sum %>% 
  subset(occ_time == 1993) %>% 
  ggplot(aes(cumulative, mean)) +
  geom_point(color = "gray50", alpha = 0.5) + 
  geom_linerange(aes(ymin = low, ymax = high), color = "gray50", alpha = 0.5) +
  geom_point(aes(cumulative, ypd), color = "black", size = 0.5) +
  facet_grid(prey ~ ., scales = "free_y")

# Model checking ===============================================================

# Data summaries
y_sum <- y_long %>% 
  ddply(.(prey), summarise,
        pzero = mean(y == 0),
        yhat = mean(y[y > 0]),
        var = var(y[y > 0])) %>% 
  join(prey_sp)

# PPD summaries
yhat_check <- rstan::extract(mcmc, "yhat")[[1]] %>% 
  melt(varnames = c("iter", "prey_idx", "bout_idx"), value.name = "y") %>% 
  ddply(.(iter, prey_idx), summarise,
        pzero = mean(y == 0),
        yhat = mean(y[y > 0]),
        var = var(y[y > 0])) %>% 
  join(prey_sp)
  
yhat_check %>%
  join(y_sum, by = "prey") %>% 
  ddply(.(prey), summarise,
        pvar = mean(var > y_sum$var[y_sum$prey == prey]))

# Proportion of zeros
yhat_check %>% 
  ggplot(aes(pzero)) + 
  geom_histogram() + 
  geom_vline(aes(xintercept = pzero), data = y_sum) +
  facet_wrap(~ prey, scales = "free_x") + 
  xlab("proportion of zeros")

# Mean of nonzero observations
yhat_check %>% 
  ggplot(aes(yhat)) + 
  geom_histogram() + 
  geom_vline(aes(xintercept = yhat), data = y_sum) +
  facet_wrap(~ prey, scales = "free_x") + 
  xlab("mean of positive counts")

# Variance of nonzero observations
yhat_check %>% 
  ggplot(aes(var)) + 
  geom_histogram() + 
  geom_vline(aes(xintercept = var), data = y_sum) +
  facet_wrap(~ prey, scales = "free_x")

# Number of prey types per bout
bout_check <- rstan::extract(mcmc, "yhat")[[1]] %>% 
  is_greater_than(0) %>% 
  apply(c(1), colSums) %>% 
  apply(c(2), mean)

bout_sum <- y_long %>% 
  ddply(.(bout_id), summarise,
        nprey = sum(y > 0))

data.frame(nprey = bout_check, iter = 1:length(bout_check)) %>% 
  ggplot(aes(nprey)) + 
  geom_histogram() + 
  geom_vline(xintercept = mean(bout_sum$nprey))
