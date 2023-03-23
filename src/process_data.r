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
data <- arrange(data, site_name, year, bout_id) %>% 
  mutate(bout_id = factor(bout_id, levels = unique(bout_id)))

# Generating names for sp x size categories
data <- data %>% 
  mutate(preysize_cat = stringr::str_sub(preysize_cd, 1, 1),
         prey_cat = stringr::str_c(SOFA_category, "_", preysize_cat))

# Lumping into SOFA taxonomic categories (summing over species and fine-grained size)
agg <- data %>% 
  ddply(.(bout_id, dive_num, prey_cat), summarise,
        prey_qty = sum(prey_qty, na.rm = TRUE))

# Reformatting to a wide format (while filling in zeros)
complete <- agg %>% 
  pivot_wider(names_from = prey_cat, values_from = prey_qty, values_fill = 0, names_sort = TRUE) 

# Keeping only major prey 
prey_sub <- complete %>% 
  dplyr::select(bout_id, dive_num,
                (starts_with("clam") | starts_with("crab") | starts_with("modiolus") | starts_with("mussel") | starts_with("urchin")) &
                !ends_with("9")) 

# Species index of each of the prey categories
sp_idx <- prey_sub %>% 
  dplyr::select(!(bout_id:dive_num)) %>% 
  names() %>% 
  str_split(pattern = "_", simplify = TRUE) %>% 
  magrittr::extract(, 1) %>% 
  as.factor() %>% 
  as.numeric()

# Spatial processing ===========================================================

# Extracting lat-long for each bout
bout_coords <- data %>% 
  ddply(.(bout_id), summarise,
        long = otter_long[1],
        lat = otter_lat[1], 
        site = site_name[1],
        year = year[1])

## Aligning with otter abundance raster

# Making spatial object from bout points
crdref <- CRS(SRS_string = "EPSG:4269")
bout_pts <- SpatialPoints(bout_coords[, 2:3], proj4string = crdref)

# Reading in mean otter abundance raster
otter_grid <- brick("data/otter_mean.grd")

# Aligning the CRS
crs(otter_grid) <- "EPSG:26708"
bout_pts_aligned <- spTransform(bout_pts, wkt(otter_grid))

# # Extracting otter density in 800m buffer
# otters <- raster::extract(otter_grid, bout_pts_aligned, buffer = 800) %>% 
#   laply(colMeans, na.rm = T)

# Computing kernel-weighted otter abundance at each bout
otter_pts <- rasterToPoints(otter_grid, spatial = TRUE)
otter_values <- extract(otter_grid, otter_pts)

distances <- pointDistance(otter_pts, bout_pts_aligned)
K <- exp(- distances ^ 2 / 1000 ^ 2)

otters <- t(otter_values) %*% K %>% t()

# Adding column of zeros for 1992 for lags and differencing
otters <- cbind(0, otters)
colnames(otters) <- c(1992:2018)

# Reformatting to dataframe (and shift year indexing by 1)
otter_df <- otters %>% 
  as.data.frame() %>% 
  mutate(bout_id = bout_coords$bout_id) %>% 
  pivot_longer(!bout_id, 
               names_to = "year", 
               values_to = "density",
               names_transform = list(year = as.integer)) %>% 
  mutate(year = year + 1) %>% 
  right_join(bout_coords)

## Bathymetry 

# (Note this was extracted from ArcGIS gdb using ArcRasterRescue:
# Barnes, Richard. 2020. Arc Raster Rescue. Software. doi: 10.5281/zenodo.4128479.)

# Depth
depth <- raster("data/bathymetry.tif")
bout_depth <- raster::extract(depth, spTransform(bout_pts, crs(depth))) 
bout_depth[bout_depth < -1000] <- NA

# Slope
slope <- raster("data/slope.tif")
bout_slope <- raster::extract(slope, spTransform(bout_pts, crs(slope)))
bout_slope[bout_slope < -1000] <- NA

## Current
current <- raster("data/current/w001001.adf") 
bout_current <- raster::extract(current, spTransform(bout_pts, crs(current)))

# Combining all covariates and dropping bouts with missing data
bout_dat <- otter_df %>% 
  mutate(depth = bout_depth,
         slope = bout_slope,
         current = bout_current) %>% 
  subset(!(is.na(depth) | is.na(slope) | is.na(current))) %>% 
  mutate(bout_idx = 1:length(bout_id))

##  Spatial covariance

# # Mapping bouts to years
# years <- unique(bout_coords$year) %>% sort()
# K_t <- daply(bout_coords, .(bout_id), function(x) years %in% x$year)
# 
# # Mapping bouts to regions
# regions <- unique(bout_coords$site) %>% sort()
# K_r <- daply(bout_coords, .(bout_id), function(x) regions %in% x$site)
# 
# K <- daply(bout_coords, .(bout_id), function(x) domains$site == x$site & domains$year == x$year)

## Calculating distance matrix for spatial covariance

# Computing pairwise bout distance matrix
d <- pointDistance(cbind(bout_dat$long, bout_dat$lat), cbind(bout_dat$long, bout_dat$lat), lonlat = TRUE, allpairs = TRUE)

# Zero-ing out correlation among sites and years
time_mask <- outer(bout_dat$year, bout_dat$year, function(x, y) x == y)
region_mask <- outer(bout_dat$site, bout_dat$site, function(x, y) x == y)

d <- d * time_mask * region_mask

# Creating distance matrix indices for site x year
domains <- bout_dat %>% 
  ddply(.(site, year), summarise,
        nbouts = length(bout_id),
        first = bout_idx[1],
        last = bout_idx[nbouts])

# Saving =======================================================================

# # Aside: looking at GB contour data from Ben
# bath <- st_read("data/bathy/glba_bathy")
# 
# bout_pts <- st_as_sf(bout_coords, coords = c("long", "lat"), crs = 4269)
# bout_pts <- st_transform(bout_pts, crs = st_crs(bath))
# 
# bout_pts <- bout_pts %>% 
#   mutate(missing = bout_id %in% missing$bout_id)
# 
# foo <- st_nearest_feature(bout_pts, bath)
# 
# depth_contour <- bath[foo, ]$DEPTH
# 
# bout_dat <- bout_dat %>% 
#   mutate(contour = depth_contour,
#          diff = abs(contour - depth))
# 
# bout_dat %>% 
#   ggplot(aes(depth, contour)) + 
#   geom_point() + 
#   geom_abline(intercept = 0, slope = 1)
# 
# subset(bout_dat, diff == max(bout_dat$diff, na.rm = TRUE))
# 
# bath %>%
#   ggplot() + 
#   geom_sf(aes(color = DEPTH), size = 0.01)  + 
#   geom_sf(data = bout_pts, inherit.aes = FALSE, size = 0.1, shape = 1)

# Subsetting to bouts with complete data
y <- subset(prey_sub, bout_id %in% bout_dat$bout_id)

# Identifying the row index at which eat new bout starts
bout_length <- ddply(y, .(bout_id), summarise,
                     nd = dive_num %>% unique %>% length) %>% 
  mutate(idx = cumsum(c(1, nd[1:(length(nd) - 1)])))

# Packaging and saving
all <- list(y = y,
            sp_idx = sp_idx,
            bout_length = bout_length,
            bout_dat = bout_dat, 
            domains = domains,
            d = d)

saveRDS(all, "data/processed.rds")
