library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(raster)
library(rgdal)
library(sf)
library(stringr)

# Reading in prey lookup
prey <- read.csv("data/LU_preytype.csv")

# Reading in foraging data and formatting dates
data <- read.csv("data/GLBA_SEOT_Forage_Data_1993-2019.csv") %>% 
  join(prey) %>% 
  mutate(bout_date = mdy(bout_date), 
         year = year(bout_date),
         month = month(bout_date),
         preysize_cd = dplyr::recode(preysize_cd, 
                                     "9Z" = "9z",
                                     "2B" = "2b", 
                                     "2C" = "2c", 
                                     "1" = "1z",
                                     "4a" = "4z"))

# Dropping outlier longitude (check if these can be corrected)
data <- subset(data, otter_long < -135.6)

# Dropping regions outside the domain of the otter survey
data <- subset(data, !(site_name %in% c("Dundas")))

# Dropping non-foraging or interrupted dives
data <- subset(data, success_cd %in% c("y", "n"))

# Dropping successful dives with missing prey quantity (note these are largely unidentified prey with unknown size)
data <- subset(data, !(success_cd == "y" & is.na(prey_qty)))

# Filling in 4 for successful dives with missing prey size category
data$preysize_cd[data$success_cd == "y" & data$preysize_cd == ""] <- 4

# Dropping bouts without coordinates
data <- subset(data, !is.na(otter_lat))

# Ordering the data and bout_ids
data <- arrange(data, year, site_name, bout_id) %>% 
  mutate(bout_id = factor(bout_id, levels = unique(bout_id)))

# Generating names for sp x size categories
data <- data %>% 
  mutate(preysize_cat = stringr::str_sub(preysize_cd, 1, 1),
         prey_cat = stringr::str_c(SOFA_category, "_", preysize_cat))

# Aggregating by bout and SOFA taxonomic categories (summing over species and fine-grained size)
agg <- data %>% 
  ddply(.(bout_id, prey_cat), summarise,
        prey_qty = sum(prey_qty, na.rm = TRUE))

# Reformatting to a wide format (while filling in zeros)
complete <- agg %>% 
  pivot_wider(names_from = prey_cat, values_from = prey_qty, values_fill = 0, names_sort = TRUE) 

# Evaluating prevalence of prey categories to set cut-offs for inclusion
pbouts <- agg %>% 
  subset(!is.na(prey_cat)) %>% 
  mutate(prey_sp = str_split(prey_cat, pattern = "_", simplify = TRUE)[, 1],
         prey_sz = str_split(prey_cat, pattern = "_", simplify = TRUE)[, 2]) %>% 
  subset(prey_sp != "unidentified") %>% 
  ddply(.(prey_cat), summarise,
        total = sum(prey_qty),
        nbouts = length(bout_id),
        pbouts = nbouts / length(unique(agg$bout_id))) %>% 
  subset(nbouts > 50)

# Keeping only major prey  (collected on at least 50 bouts)
prey_sub <- complete %>% 
  dplyr::select(bout_id, pbouts$prey_cat) 

# Spatial processing ===========================================================

# Extracting lat-long for each bout
bout_coords <- data %>% 
  ddply(.(bout_id), summarise,
        long = otter_long[1],
        lat = otter_lat[1], 
        site = site_name[1],
        bout_year = year[1])

## Aligning with otter abundance raster

# Making spatial object from bout points
crdref <- CRS(SRS_string = "EPSG:4269")
bout_pts <- SpatialPoints(bout_coords[, 2:3], proj4string = crdref)

# Reading in mean otter abundance raster
otter_grid <- brick("data/otter_mean.grd")

# Aligning the CRS
crs(otter_grid) <- "EPSG:26708"
bout_pts_aligned <- spTransform(bout_pts, wkt(otter_grid))

# Computing kernel-weighted otter abundance at each bout
otter_pts <- rasterToPoints(otter_grid, spatial = TRUE)
otter_values <- extract(otter_grid, otter_pts)

distances <- pointDistance(otter_pts, bout_pts_aligned)
K <- exp(- distances ^ 2 / 500 ^ 2)

otters <- t(otter_values) %*% K %>% t()

# Adding column of zeros for 1992 for lags and differencing
otters <- cbind(0, otters)
colnames(otters) <- c(1992:2018)

# Reformatting to dataframe
otter_df <- otters %>% 
  as.data.frame() %>% 
  mutate(bout_id = bout_coords$bout_id) %>% 
  pivot_longer(!bout_id, 
               names_to = "year", 
               values_to = "density",
               names_transform = list(year = as.integer)) %>% 
  left_join(bout_coords)

# Otter metrics
otter_x <- otter_df %>% 
  ddply(.(bout_id), summarise,
        lagged = density[year == (bout_year[1] - 1)],
        cumulative = sum(density[year < bout_year[1]]),
        occ_time = min(year[density > 0.25]),
        end = max(density)) %>% 
  left_join(bout_coords) %>% 
  mutate(occ_time = ifelse(occ_time < Inf, occ_time, bout_year))

# ## Bathymetry 
# 
# # (Note this was extracted from ArcGIS gdb using ArcRasterRescue:
# # Barnes, Richard. 2020. Arc Raster Rescue. Software. doi: 10.5281/zenodo.4128479.)
# 
# # Depth
# depth <- raster("data/bathymetry.tif")
# bout_depth <- raster::extract(depth, spTransform(bout_pts, crs(depth))) 
# bout_depth[bout_depth < -1000] <- NA
# 
# # Slope
# slope <- raster("data/slope.tif")
# bout_slope <- raster::extract(slope, spTransform(bout_pts, crs(slope)))
# bout_slope[bout_slope < -1000] <- NA
# 
# ## Current
# current <- raster("data/current/w001001.adf") 
# bout_current <- raster::extract(current, spTransform(bout_pts, crs(current)))
# 
# # Combining all covariates and dropping bouts with missing data
# bout_dat <- otter_x %>% 
#   mutate(depth = bout_depth,
#          slope = bout_slope,
#          current = bout_current) %>% 
#   subset(!(is.na(depth) | is.na(slope) | is.na(current))) %>% 
#   mutate(bout_idx = 1:length(bout_id))

##  Spatial random effects

# Mapping bouts to regions
regions <- unique(bout_coords$site) %>% sort()
K_r <- daply(bout_coords, .(bout_id), function(x) regions %in% x$site)

# Mapping bouts to years
years <- unique(bout_coords$bout_year) %>% sort()
K_t <- daply(bout_coords, .(bout_id), function(x) years %in% x$bout_year)

# Computing pairwise bout distance matrix
d <- pointDistance(cbind(otter_x$long, otter_x$lat), cbind(otter_x$long, otter_x$lat), lonlat = TRUE, allpairs = TRUE)

# Saving =======================================================================

# Computing number of dives in each bout and combining with other bout data
bouts <- ddply(data, .(bout_id), summarise,
               nd = dive_num %>% unique %>% length) %>% 
  join(otter_x) %>% 
  mutate(bout_idx = 1:nrow(otter_x))

# Packaging and saving
all <- list(y = prey_sub,
            bout_dat = bouts, 
            d = d,
            K_r = K_r,
            K_t = K_t)

saveRDS(all, "data/processed.rds")
