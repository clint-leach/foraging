library(magrittr)
library(plyr)
library(lubridate)
library(ggplot2)
library(rgdal)
library(sf)

# Reading in prey lookup
prey <- read.csv("data/LU_preytype.csv")

# Reading in foraging data and formatting dates
data <- read.csv("data/GLBA_SEOT_Forage_Data_1993-2019.csv") %>% 
  mutate(bout_date = mdy(bout_date), 
         year = year(bout_date),
         month = month(bout_date)) %>% 
  join(prey)


# Reading in Glacier Bay shapefile for mapping
gb <- st_read(dsn = "data/PH6502/historicl1.shp")

# Summaries ====================================================================

# How many bouts
data$bout_id %>% unique() %>% length()

# Dives by year and site
data %>% 
  ddply(.(year, site_name), summarise,
        total = length(success_cd)) %>% 
  ggplot(aes(site_name, year, label = total)) + 
  geom_label() +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  xlab("site") + 
  ggtitle("number of observed dives")

# Missing otter coords
sum(is.na(data$otter_lat))

# Missing observer coords
sum(is.na(data$obs_lat))

# Missing prey type
sum(data$preytype_cd == "" & data$success_cd == "y")

# Missing prey size
sum(data$preysize_cd == "" & data$success_cd == "y")

# Unidentified prey type
sum(data$preytype_cd == "uni" & data$success_cd == "y")

# How many successful dives
sum(data$success_cd == "y")

# Dive success codes
data %>% 
  ggplot(aes(success_cd)) + 
  geom_bar()

# Where are unsuccessful dives distributed within a bout?
bout_sum <- ddply(data, .(bout_id), summarise,
                  ndives = max(dive_num))
data %>% 
  join(bout_sum) %>% 
  subset(success_cd == "n") %>% 
  subset(ndives > 1) %>% 
  mutate(bout_pos = dive_num / ndives) %>% 
  ggplot(aes(bout_pos)) + 
  geom_histogram()

# How often do bouts end on an unsuccessful dive?
data %>% 
  join(bout_sum) %>% 
  subset(dive_num == ndives) %>% 
  ggplot(aes(success_cd)) + 
  geom_bar()

# Multiple prey ================================================================
 
bout_ids = unique(data$bout_id)

# How many dives have a prey type with multiple individuals
data %>% 
  ggplot(aes(prey_qty)) + 
  geom_bar()

# Prey quantity by prey type
data %>% 
  ggplot(aes(prey_qty)) + 
  geom_histogram() + 
  facet_wrap(~SOFA_category)

# Total number of prey items per dive (summed across type)
data %>% 
  subset(success_cd %in% c("y", "n")) %>% 
  mutate(prey_qty = ifelse(is.na(prey_qty), 0, prey_qty)) %>% 
  ddply(.(bout_id, dive_num), summarise,
        nprey = sum(prey_qty)) %>% 
  ggplot(aes(nprey)) + 
  geom_bar()

data %>% 
  subset(success_cd %in% c("y", "n")) %>% 
  subset(bout_id %in% sample(bout_ids, 20)) %>% 
  mutate(prey_qty = ifelse(is.na(prey_qty), 0, prey_qty)) %>% 
  ddply(.(bout_id, dive_num), summarise,
        nprey = sum(prey_qty)) %>% 
  ggplot(aes(nprey)) + 
  geom_bar() + 
  facet_wrap(~bout_id)

# How many dives have multiple types of prey (excluding unidentified)
data %>% 
  subset(SOFA_category != "unidentified") %>% 
  ddply(.(bout_id, dive_num), summarise,
        ntypes = length(unique(SOFA_category))) %>% 
  ggplot(aes(ntypes)) + 
  geom_bar()

nprey <- data %>% 
  ddply(.(bout_id, dive_num), summarise,
        ntypes = length(unique(SOFA_category)))

data <- join(data, nprey)

subset(data, ntypes > 1) %>% 
  ggplot(aes(SOFA_category)) + 
  geom_bar()

# Predominant prey categories ==================================================

data %>% 
  subset(success_cd == "y") %>% 
  ddply(.(preytype_cd), summarise,
        total = length(bout_id)) %>% 
  mutate(preytype_cd = factor(preytype_cd, levels = preytype_cd[order(total, decreasing = TRUE)])) %>% 
  ggplot(aes(preytype_cd, total)) + 
  geom_col() + 
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  scale_y_continuous(expand = c(0, 0))

# Prey categories by groups
data %>% 
  ddply(.(SOFA_category), summarise,
        total = length(bout_id)) %>% 
  mutate(SOFA_category = factor(SOFA_category, levels = SOFA_category[order(total, decreasing = TRUE)])) %>% 
  ggplot(aes(SOFA_category, total)) + 
  geom_col() + 
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  scale_y_continuous(expand = c(0, 0))
  

# Basic mapping ================================================================

sites <- read.csv("data/sites.csv")

ggplot(gb) + 
  geom_sf() + 
  geom_point(aes(otter_long, otter_lat, color = site_name), data = data, size = 0.1) + 
  geom_point(aes(longitude, latitude), data = sites, pch = 4, color = "black", size = 2.0) +
  theme_classic()

subset(data, otter_long > -135.6)

# Plotting observation sites for a single prey species
ggplot(gb) + 
  geom_sf() + 
  geom_point(aes(otter_long, otter_lat), data = subset(data, preytype_cd == "std"), size = 1.0, color = "blue") + 
  # geom_point(aes(longitude, latitude), data = sites, pch = 4, color = "black", size = 2.0) +
  theme_classic() + 
  facet_wrap(~year)

# Success rates ================================================================

data %>% 
  ddply(.(year, site_name), summarise,
        success = sum(success_cd == "y") / length(success_cd)) %>% 
  ggplot(aes(year, success)) + 
  geom_point() +
  geom_line() + 
  facet_grid(site_name ~.)

# Foraging effort and dive times ===============================================

# Overall number of dives per bout
bouts <- data %>% 
  ddply(.(bout_id), summarise,
        ndives = length(dive_num), 
        year = year[1],
        site_name = site_name[1])

bouts %>% 
  ggplot(aes(ndives)) + 
  geom_histogram()


# average dive time by site and time
data %>% 
  subset(success_cd = "y") %>% 
  ddply(.(year, site_name), summarise,
        avg_dive = mean(dive_time, na.rm = T)) %>% 
  ggplot(aes(year, avg_dive, group = site_name)) + 
  geom_point() + 
  facet_wrap(~site_name)

# average dive time by site and time and bout
data %>% 
  subset(success_cd = "y") %>% 
  ddply(.(bout_id), summarise,
        avg_dive = mean(dive_time, na.rm = T),
        year = year[1],
        site_name = site_name[1]) %>% 
  ggplot(aes(year, avg_dive, group = site_name)) + 
  geom_point() + 
  facet_wrap(~site_name)

# average bout length
data %>% 
  ddply(.(bout_id), summarise,
        bout_length = sum(dive_time, na.rm = T), 
        year = year[1],
        site_name = site_name[1]) %>% 
  ddply(.(year, site_name), summarise,
        avg_length = mean(bout_length)) %>% 
  ggplot(aes(year, avg_length)) + 
  geom_point() + 
  facet_wrap(~site_name)

# overall average dive time by prey
data %>% 
  ddply(.(SOFA_category), summarise, 
        avg_dive = mean(dive_time, na.rm = T),
        sd = sd(dive_time, na.rm = T),
        n = length(bout_id)) %>% 
  ggplot(aes(SOFA_category, avg_dive)) + 
  geom_point() + 
  geom_linerange(aes(ymin = avg_dive - sd, ymax = avg_dive + sd)) +
  geom_label(aes(SOFA_category, avg_dive + sd + 10, label = n)) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# average dive time by prey through time
data %>% 
  ddply(.(year, SOFA_category), summarise, 
        avg_dive = mean(dive_time, na.rm = T)) %>% 
  ggplot(aes(year, avg_dive)) + 
  geom_point() + 
  geom_line() +
  facet_wrap(~SOFA_category, scales = "free_y")

# average dive time by prey across sites and time
data %>% 
  ddply(.(year, SOFA_category, site_name), summarise, 
        avg_dive = mean(dive_time, na.rm = T)) %>% 
  ggplot(aes(year, avg_dive)) + 
  geom_point() + 
  geom_line() +
  facet_grid(site_name~SOFA_category, scales = "free_y")

# Prey switching within bouts ==================================================

multiple_dives <- subset(bouts, ndives > 1)

preyseq <- data %>% 
  subset(bout_id %in% multiple_dives$bout_id) %>%
  ddply(.(bout_id), summarise,
        prey = SOFA_category[1:(length(SOFA_category) - 1)],
        followed_by = SOFA_category[2:length(SOFA_category)]) %>% 
  ddply(.(prey, followed_by), summarise,
        count = length(prey)) 

# Counts of prey type following capture of each prey type
preyseq %>% 
  ggplot(aes(0, 0, label = count)) + 
  geom_text() + 
  facet_grid(prey ~ followed_by) + 
  theme_classic()

preyseq %>% 
  ggplot(aes(followed_by, count)) + 
  geom_col() + 
  facet_grid(prey ~., scales = "free_y") + 
  theme_classic()

clamfirst <- data %>% 
  subset(dive_num == 1 & SOFA_category == "clam")

clambouts <- data %>% 
  subset(bout_id %in% clamfirst$bout_id) %>% 
  mutate(isclam = SOFA_category == "clam")

clambouts %>% 
  subset(!is.na(SOFA_category)) %>% 
  subset(SOFA_category != "unidentified") %>% 
  ddply(.(bout_id), summarise,
        pclam = sum(isclam) / length(isclam)) %>% 
  ggplot(aes(pclam)) +
  geom_histogram()

clambouts %>% 
  subset(!is.na(SOFA_category)) %>% 
  subset(SOFA_category != "unidentified") %>% 
  subset(bout_id == "20100709BPW351") %>% 
  ggplot(aes(dive_num, SOFA_category)) + 
  geom_point()

  
  
# Within bout diet diversity ===================================================

# Number of prey types per bout
data %>% 
  subset(!is.na(SOFA_category)) %>% 
  subset(SOFA_category != "unidentified") %>% 
  ddply(.(bout_id), summarise,
        nprey = length(unique(SOFA_category))) %>% 
  ggplot(aes(nprey)) + 
  geom_bar()

# Number of prey types per bout through time
data %>% 
  subset(!is.na(SOFA_category)) %>% 
  subset(SOFA_category != "unidentified") %>% 
  ddply(.(year, site_name, bout_id), summarise,
        nprey = length(unique(SOFA_category))) %>% 
  ddply(.(year, site_name), summarise,
        nprey = mean(nprey)) %>% 
  ggplot(aes(year, nprey)) + 
  geom_point() + 
  facet_wrap(~site_name)

# Local variation in diet within a sampling event ==============================

data %>% 
  ddply(.(year, site_name), summarise,
        ndays = length(unique(bout_date)))

# Looking at temporal trends in prey composition ===============================

ndives <- ddply(data, .(year, site_name), summarise,
                ndives = length(bout_id))

nprey <- ddply(data, .(year, site_name, preytype_cd), summarise,
               nobs = length(success_cd)) %>% 
  join(ndives) %>% 
  mutate(pdives = nobs/ ndives)


# Looking at a subset of species across time and sites
nprey %>% 
  subset(preytype_cd %in% c("sag", "std", "mtr", "mom", "cla", "mus", "myt")) %>% 
  ggplot(aes(year, pdives)) + 
  geom_point() +
  geom_line(color = "black", alpha = 0.5) +
  facet_grid(preytype_cd ~ site_name) + 
  scale_color_viridis_c() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  theme(strip.text.x = element_text(angle = 90, hjust = 0),
        strip.background = element_blank()) + 
  ylab("proportion of dives")

nprey %>% 
  subset(preytype_cd %in% c("sag", "std", "mtr", "mom", "cla", "mus", "myt")) %>% 
  ggplot(aes(year, pdives, color = preytype_cd)) + 
  geom_point() +
  geom_line() +
  facet_wrap(~site_name) + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab("proportion of dives")

# Looking at all species across time at a single site
nprey %>% 
  subset(site_name == "Berg") %>% 
  # subset(preytype_cd %in% c("", "uni", "cla", "mom", "myt", "sag", "std")) %>% 
  ggplot(aes(preytype_cd, ymax = pdives)) + 
  geom_linerange(ymin = 0, size = 2) +
  facet_grid(year ~.)

# Looking at all species and sites in a single year
nprey %>% 
  subset(year == 2018) %>% 
  ggplot(aes(preytype_cd, ymax = pdives)) + 
  geom_linerange(ymin = 0, size = 2) + 
  facet_grid(site_name ~.)

# Looking at time series of each species
nprey %>% 
  ggplot(aes(year, pdives, group = site_name)) + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~preytype_cd)

# Looking at community data ====================================================

# Summing over quadrats and size
comm <- read.csv("data/zero_augmented.csv") %>% 
  ddply(.(site, year, spp), summarise,
        total = sum(count),
        nquad = length(unique(quad)),
        density = total / nquad)

# Looking at species dynamics across a subset of sites
comm %>% 
  # subset(spp %in% c("STD", "SAG")) %>% 
  subset(site %in% c("Puffin", "Triangle_PCH", "Berg", "67", "Boulder_PCH", "Strawberry", "Secret", "Secret_PCH")) %>% 
  ggplot(aes(year, density, color = spp)) + 
  geom_line() + 
  facet_grid(site~spp, scales = "free")

# Looking at a subset of species at all sites
comm %>% 
  subset(spp %in% c("LES", "SAG", "STD")) %>% 
  ggplot(aes(year, density, color = spp)) + 
  geom_line() + 
  facet_wrap(~site, scales = "fixed")

# Summing over quadrats to get size distributions
sizedist <- read.csv("data/zero_augmented.csv") %>% 
  ddply(.(site, year, size, spp), summarise,
        total = sum(count),
        nquad = length(unique(quad)),
        density = total / nquad)

# Looking at size dist dynamics for a subset of species at subset of sites
sizedist %>% 
  subset(spp %in% c("SAG")) %>% 
  subset(site %in% c("Boulder_PCH", "Strawberry", "Secret", "Secret_PCH", "Puffin", "Triangle_PCH", "Berg", "67")) %>%
  ggplot(aes(size, density, color = spp, group = spp)) + 
  geom_line() + 
  geom_vline(xintercept = 25) +
  facet_grid(year ~ site, scales = "fixed")

# Aggregating urchins into 1cm bins and looking at dynamics of each bin
sizedist %>% 
  subset(spp == "STD") %>% 
  mutate(cm = floor(size / 10)) %>% 
  ddply(.(site, year, cm), summarise,
        density = sum(density)) %>% 
  ggplot(aes(year, density, group = site)) + 
  geom_line(alpha = 0.5) + 
  # geom_vline(xintercept = 2.5) +
  facet_grid(cm ~., scales = "free_y")

 # Aggregating Saxidomus into 1cm bins and looking at dynamics at each site
sizedist %>%
  subset(spp == "LES") %>% 
  mutate(cm = floor(size / 10)) %>% 
  ddply(.(site, year, spp, cm), summarise,
        density = sum(density)) %>% 
  ggplot(aes(cm, density, color = year, group = year)) + 
  geom_line() + 
  facet_wrap(~site, scales = "fixed")

# Saxidomus dynamics at unresponsive sites
sizedist %>% 
  subset(spp == "SAG") %>% 
  subset(site %in% c("Puffin", "Triangle_PCH", "Berg", "67")) %>% 
  ggplot(aes(size, density, color = spp, group = spp)) + 
  geom_line() + 
  facet_grid(year ~ site, scales = "fixed")
