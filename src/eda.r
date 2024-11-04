library(magrittr)
library(plyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(sf)

# Reading in prey lookup
prey <- read.csv("data/seaOtterForage_GlacierBay_Kloecker/GLBA_SEOT_forage_taxonomy_1993-2019.csv")

# Grouping prey into taxonomic classes
classify_prey <- function(phylum, class, order, family, genus){
  
  category <- "other"
  
  if(class == "Bivalvia" & order != "Mytiloida" & family != "Pectinidae"){
    category <- "clam"
  }
  
  if(genus == "Modilous"){
    category <- "modiolus"
  }
  
  if(order == "Mytiloida" & genus != "Modilous"){
    category <- "mytilus"
  }
  
  if(family == "Pectinidae"){
    category <- "scallop"
  }
  
  if(order == "Decapoda" & genus != "Pandalus"){
    category <- "crab"
  }
  
  if(class == "Gastropoda"){
    category <- "snail"
  }
  
  if(class %in% c("Asteroidea", "Ophiuroidea")){
    category <- "star"
  }
  
  if(class == "Echinoidea"){
    category <- "urchin"
  }
  
  if(class == "Polyplacophora"){
    category <- "chiton"
  }
  
  if(class == "Holothuroidea"){
    category <- "cucumber"
  }
  
  if(phylum %in% c("Annelida", "Sipuncula", "Annelida, Sipuncula")){
    category <- "worm"
  }
  
  if(phylum == ""){
    category <- "unidentified"
  }
  
  return(category)
}

prey <- prey %>% 
  rowwise() %>% 
  mutate(category = classify_prey(phylum, class, order, family, genus))

# Reading in foraging data and formatting dates
data <- read.csv("data/seaOtterForage_GlacierBay_Kloecker/GLBA_SEOT_forage_bouts_1993-2019.csv") %>% 
  join(prey) %>% 
  mutate(bout_date = ymd(bout_date), 
         year = year(bout_date),
         month = month(bout_date),
         preysize_cd = dplyr::recode(preysize_cd, 
                                     "9Z" = "9z",
                                     "2B" = "2b", 
                                     "2C" = "2c", 
                                     "1" = "1z",
                                     "4a" = "4z"))

# Dropping regions outside the domain of the otter survey
data <- subset(data, !(site_name %in% c("Dundas", "Althorp", "Inian", "Lemesurier", "")))

# Giving each dive a unique id
data <- mutate(data, dive_id = stringr::str_c(bout_id, "_", dive_num))

# Reading in Glacier Bay shapefile for mapping
gb <- st_read(dsn = "data/PH6502/historicl1.shp")

# Breakdown by sex and pup status ==============================================

obs_otter <- data %>% 
  ddply(.(bout_id), summarise,
        sex_cd = sex_cd[1],
        pupsize_cd = pupsize_cd[1])

obs_otter %>% 
  subset(sex_cd != "u ") %>% 
  ggplot(aes(sex_cd)) + 
  geom_bar() + 
  ylab("number of bouts") + 
  theme_bw() + 
  scale_y_continuous(expand = c(0, 10))

obs_otter %>% 
  subset(sex_cd == "f") %>%
  ggplot(aes(pupsize_cd)) + 
  geom_bar() + 
  ylab("number of bouts") + 
  theme_bw() + 
  scale_y_continuous(expand = c(0, 10))

# Summaries ====================================================================

# How many bouts
data$bout_id %>% unique() %>% length()

# How many dives
dive_sum <- data %>% 
  ddply(.(bout_id), summarise,
        ndives = length(unique(dive_num)),
        year = year[1],
        site_name = site_name[1])

dive_sum %>% 
  dplyr::select(ndives) %>% 
  sum()

site_dives <- ddply(dive_sum, .(site_name, year), summarise,
                    ndives = sum(ndives),
                    nbouts = length(unique(bout_id)))

# Bouts by year
data %>% 
  ddply(.(year), summarise,
        total = length(unique(bout_id))) %>% 
  ggplot(aes(year, total)) + 
  geom_point() +
  geom_line() +
  theme_classic() +
  xlab("year") + 
  ylab("number of observed bouts")

# Bouts by year and site
data %>% 
  ddply(.(year, site_name), summarise,
      total = length(unique(bout_id))) %>% 
  ggplot(aes(year, total)) + 
  geom_point() + 
  facet_wrap(~site_name) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ylab("number of observed bouts")


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

# Missing dive coords
sum(is.na(data$dive_otter_lat))

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

# Overall distribution of prey quantity
data %>% 
  ggplot(aes(prey_qty)) + 
  geom_bar()

# Prey quantity by prey type
data %>% 
  ggplot(aes(prey_qty)) + 
  geom_histogram() + 
  facet_wrap(~category)

# Prey quantity by prey size
data %>% 
  ggplot(aes(prey_qty)) + 
  geom_histogram() + 
  facet_wrap(~preysize_cd)

# Total number of prey items per dive (summed across type)
data %>% 
  subset(success_cd %in% c("y", "n")) %>% 
  mutate(prey_qty = ifelse(is.na(prey_qty), 0, prey_qty)) %>% 
  ddply(.(bout_id, dive_num), summarise,
        nprey = sum(prey_qty)) %>% 
  ggplot(aes(nprey)) + 
  geom_bar() + 
  xlab("prey quantity") + 
  theme_classic()

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
  subset(category != "unidentified") %>% 
  ddply(.(bout_id, dive_num), summarise,
        ntypes = length(unique(category))) %>% 
  ggplot(aes(ntypes)) + 
  geom_bar()

nprey <- data %>% 
  ddply(.(bout_id, dive_num), summarise,
        ntypes = length(unique(category)))

data <- join(data, nprey)

subset(data, ntypes > 1) %>% 
  ggplot(aes(category)) + 
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
  ddply(.(category), summarise,
        total = length(bout_id)) %>% 
  mutate(category = factor(category, levels = category[order(total, decreasing = TRUE)])) %>% 
  ggplot(aes(category, total)) + 
  geom_col() + 
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("SOFA category") + 
  ylab("number of dives")
  
# Proportion of all prey accounted for by main prey types
main <- data %>% 
  subset(category %in% c("clam", "crab", "mussel", "modiolus", "urchin"))

sum(main$prey_qty, na.rm = TRUE) / sum(data$prey_qty, na.rm = TRUE)

# Number of dives each prey observed on
total_dives <- data$dive_id %>% unique() %>% length()

data %>% 
  ddply(.(category), summarise,
        pbouts = length(unique(bout_id)) / length(unique(data$bout_id))) %>% 
  arrange(desc(pbouts)) %>% 
  as_tibble()

prey_dives <- data %>% 
  ddply(.(preytype_cd), summarise,
        pbouts = length(unique(bout_id)) / length(unique(data$bout_id))) %>% 
  join(prey) %>% 
  arrange(category) %>% 
  as_tibble()

prey_dives %>% subset(category == "crab") %>% 
  arrange(desc(pbouts)) %>% 
  print(n = Inf)

# Basic mapping ================================================================

ggplot(gb) + 
  geom_sf() + 
  geom_point(aes(otter_long, otter_lat, color = site_name), data = data, size = 1.0) + 
  theme_classic() + 
  scale_color_discrete()

# Plotting observation sites for a single prey species
ggplot(gb) + 
  geom_sf() + 
  geom_point(aes(otter_long, otter_lat), data = subset(data, preytype_cd == "std"), size = 1.0, color = "blue") + 
  theme_classic() + 
  facet_wrap(~year)

# Plotting observations for a single site
ggplot(gb) + 
  geom_sf() + 
  geom_point(aes(otter_long, otter_lat, color = factor(year)), data = subset(data, site_name == "Geikie"), size = 2) + 
  theme_classic()

ggplot(gb) + 
  geom_sf() + 
  geom_point(aes(otter_long, otter_lat, color = site_name), data = subset(data, site_name %in% c("Scidmore/HughMiller", "Russel/Reid")), size = 1.0) + 
  theme_classic()


# Plotting observations through time (facets)
ggplot(gb) + 
  geom_sf() + 
  geom_point(aes(otter_long, otter_lat, color = site_name), data = subset(data, otter_long < -135.6), size = 2.0) +
  facet_wrap(~year) + 
  theme_classic()

# Plotting observations through time (color)
ggplot(gb) + 
  geom_sf() + 
  geom_point(aes(otter_long, otter_lat, color = year), data = subset(data, otter_long < -135.6 & site_name != "Dundas"), size = 0.2) +
  scale_color_distiller(palette = "YlGnBu", direction = 1) +
  theme_classic()

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
  subset(success_cd == "y") %>% 
  ddply(.(year, site_name), summarise,
        avg_dive = mean(dive_time, na.rm = T)) %>% 
  ggplot(aes(year, avg_dive, group = site_name)) + 
  geom_point() + 
  facet_wrap(~site_name)

# average dive time by site and time and bout
data %>% 
  subset(success_cd == "y") %>% 
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
  ddply(.(category), summarise, 
        avg_dive = mean(dive_time, na.rm = T),
        sd = sd(dive_time, na.rm = T),
        n = length(bout_id)) %>% 
  ggplot(aes(category, avg_dive)) + 
  geom_point() + 
  geom_linerange(aes(ymin = avg_dive - sd, ymax = avg_dive + sd)) +
  geom_label(aes(category, avg_dive + sd + 10, label = n)) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# average dive time by prey through time
data %>% 
  ddply(.(year, category), summarise, 
        avg_dive = mean(dive_time, na.rm = T)) %>% 
  ggplot(aes(year, avg_dive)) + 
  geom_point() + 
  geom_line() +
  facet_wrap(~category, scales = "free_y")

# average dive time by prey across sites and time
data %>% 
  ddply(.(year, category, site_name), summarise, 
        avg_dive = mean(dive_time, na.rm = T)) %>% 
  ggplot(aes(year, avg_dive)) + 
  geom_point() + 
  geom_line() +
  facet_grid(site_name~category, scales = "free_y")

# Prey switching within bouts ==================================================

multiple_dives <- subset(bouts, ndives > 1)

preyseq <- data %>% 
  subset(bout_id %in% multiple_dives$bout_id) %>%
  ddply(.(bout_id), summarise,
        prey = category[1:(length(category) - 1)],
        followed_by = category[2:length(category)]) %>% 
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
  subset(dive_num == 1 & category == "clam")

clambouts <- data %>% 
  subset(bout_id %in% clamfirst$bout_id) %>% 
  mutate(isclam = category == "clam")

clambouts %>% 
  subset(!is.na(category)) %>% 
  subset(category != "unidentified") %>% 
  ddply(.(bout_id), summarise,
        pclam = sum(isclam) / length(isclam)) %>% 
  ggplot(aes(pclam)) +
  geom_histogram()

# Within bout diet diversity ===================================================

# Number of prey types per bout
preyperbout <- data %>% 
  subset(!is.na(category)) %>%
  subset(category != "unidentified") %>%
  ddply(.(bout_id), summarise,
        nprey = length(unique(category))) 

preyperbout %>% 
  ggplot(aes(nprey)) + 
  geom_bar()

mean(preyperbout$nprey)

# Number of prey types per bout through time
data %>% 
  subset(!is.na(category)) %>% 
  subset(category != "unidentified") %>% 
  ddply(.(year, site_name, bout_id), summarise,
        nprey = length(unique(category))) %>% 
  ddply(.(year, site_name), summarise,
        nprey = mean(nprey)) %>% 
  ggplot(aes(year, nprey)) + 
  geom_point() + 
  facet_wrap(~site_name)


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

# Size dynamics within species =================================================

data <- data %>% 
  mutate(preysize_cat = stringr::str_sub(preysize_cd, 1, 1),
         prey_cat = stringr::str_c(category, "_", preysize_cat))

totals <- data %>% 
  ddply(.(year, site_name), summarise,
        total_prey = sum(prey_qty, na.rm = TRUE),
        total_sag = sum(prey_qty[preytype_cd == "sag"], na.rm = TRUE),
        prop_sag = total_sag / total_prey) %>% 
  join(site_dives)

# Overall proportion of saxidomus through time
totals %>% 
  ggplot(aes(year, prop_sag, alpha = log10(total_prey))) + 
  geom_line() + 
  geom_point() +
  facet_wrap(~site_name)


# Summarizing and filling in zeros for each dive
sag <- subset(data, preytype_cd == "sag") %>% 
  ddply(.(site_name, year, bout_id, dive_id, preysize_cat), summarise,
        prey_qty = sum(prey_qty, na.rm = TRUE)) %>% 
  pivot_wider(id_cols = c(site_name, year, bout_id, dive_id), 
              names_from = preysize_cat,
              values_from = prey_qty,
              values_fill = 0) %>% 
  pivot_longer(!(site_name | year | bout_id | dive_id), names_to = "preysize_cat", values_to = "count")

sag_sum <- sag %>% 
  ddply(.(site_name, year, preysize_cat), summarise,
        total = sum(count, na.rm = TRUE),
        mean_captured = mean(count[count > 0], na.rm = TRUE),
        unique_dives = length(unique(dive_id[count > 0])),
        unique_bouts = length(unique(bout_id[count > 0]))) %>% 
  join(totals) %>% 
  mutate(prop_of_sag = total / total_sag,
         prop_of_all = total / total_prey,
         prop_of_dives = unique_dives / ndives,
         prop_of_bouts = unique_bouts / nbouts,
         mean_per_dive = total / ndives)

# Proportion of saxidomus in each size class through time
sag_sum %>% 
  subset(preysize_cat != "9") %>% 
  ggplot(aes(preysize_cat, prop_of_sag, color = year, group = year)) + 
  geom_point() +
  geom_line() +
  scale_color_viridis_c() +
  facet_wrap(~ site_name) + 
  theme_classic()

sag_sum %>% 
  subset(preysize_cat != "9") %>% 
  ggplot(aes(year, prop_of_sag, color = site_name, group = site_name)) + 
  geom_point() +
  geom_line() +
  facet_grid(preysize_cat ~ .) + 
  theme_classic()

sag_sum %>% 
  subset(preysize_cat != "9") %>% 
  ggplot(aes(year, prop_of_sag, color = preysize_cat, group = preysize_cat)) + 
  geom_point() +
  geom_line() +
  facet_wrap(~ site_name) + 
  theme_classic()

# Mean count of saxidomus per dive in each size class through time
sag_sum %>% 
  subset(preysize_cat != "9") %>% 
  ggplot(aes(preysize_cat, mean_per_dive, color = year, group = year)) + 
  geom_point() +
  geom_line() +
  scale_color_viridis_c() +
  facet_wrap(~ site_name) + 
  theme_classic()

sag_sum %>% 
  subset(preysize_cat != "9") %>% 
  ggplot(aes(year, mean_per_dive, color = site_name, group = site_name)) + 
  geom_point() +
  geom_line() +
  facet_grid(preysize_cat ~ .) + 
  theme_classic()

sag_sum %>% 
  subset(preysize_cat != "9") %>% 
  ggplot(aes(year, mean_per_dive, color = preysize_cat, group = preysize_cat)) + 
  geom_point() +
  geom_line() +
  facet_wrap(~ site_name) + 
  theme_classic()

# Proportion of dives with saxidomus of given size
sag_sum %>% 
  subset(preysize_cat != "9") %>% 
  ggplot(aes(preysize_cat, prop_of_dives, color = year, group = year)) + 
  geom_point() +
  geom_line() +
  scale_color_viridis_c() +
  facet_wrap(~ site_name) + 
  theme_classic()

sag_sum %>% 
  subset(preysize_cat != "9") %>% 
  ggplot(aes(year,prop_of_dives, color = site_name, group = site_name)) + 
  geom_point() +
  geom_line() +
  facet_grid(preysize_cat ~ .) + 
  theme_classic()

sag_sum %>% 
  subset(preysize_cat != "9") %>% 
  ggplot(aes(year, prop_of_dives, color = preysize_cat, group = preysize_cat)) + 
  geom_point() +
  geom_line() +
  facet_wrap(~ site_name) + 
  theme_classic()

# Proportion of bouts with saxidomus of given size
sag_sum %>% 
  subset(preysize_cat != "9") %>% 
  ggplot(aes(preysize_cat, prop_of_bouts, color = year, group = year)) + 
  geom_point() +
  geom_line() +
  scale_color_viridis_c() +
  facet_wrap(~ site_name) + 
  theme_classic()

sag_sum %>% 
  subset(preysize_cat != "9") %>% 
  ggplot(aes(year, prop_of_bouts, color = site_name, group = site_name)) + 
  geom_point() +
  geom_line() +
  facet_grid(preysize_cat ~ .) + 
  theme_classic()

sag_sum %>% 
  subset(preysize_cat != "9") %>% 
  ggplot(aes(year, prop_of_bouts, color = preysize_cat, group = preysize_cat)) + 
  geom_point() +
  geom_line() +
  facet_wrap(~ site_name) + 
  theme_classic()

# Count per dive given selected
sag_sum %>% 
  subset(preysize_cat != "9") %>% 
  ggplot(aes(year, mean_captured, color = preysize_cat, group = preysize_cat)) + 
  geom_point() +
  geom_line() +
  facet_wrap(~ site_name) + 
  theme_classic()

# Average proportion of dives within bout (on bouts where prey appears)

sag_bouts <- data %>% 
  subset(preytype_cd == "sag") %>% 
  ddply(.(site_name, year, bout_id, dive_id, preysize_cat), summarise,
        prey_qty = sum(prey_qty, na.rm = TRUE)) %>% 
  ddply(.(site_name, year, bout_id, preysize_cat), summarise,
        sag_dives = length(unique(dive_id))) %>% 
  join(bout_sum) %>% 
  mutate(pdives = sag_dives / ndives) %>% 
  ddply(.(site_name, year, preysize_cat), summarise,
        mean_pdive = mean(pdives))

sag_bouts %>% 
  subset(preysize_cat != "9") %>% 
  ggplot(aes(year, mean_pdive, color = preysize_cat, group = preysize_cat)) + 
  geom_point() +
  geom_line() +
  facet_wrap(~ site_name) + 
  theme_classic()

