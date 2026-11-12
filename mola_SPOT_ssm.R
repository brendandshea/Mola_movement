# SETUP ----

library(tidyverse)
library(aniMotum)
library(sf)
library(scales)
library(patchwork)
library(zoo)
library(ggpubr)
library(marmap)
library(rnaturalearth)

## read in data and wrangle ----
df <- read.csv("mola_all.csv") %>%
  mutate( 
    date = mdy_hms(date, tz = "UTC"),
    id = as.character(id),
    tag_time = case_when( #set tag times to drop first 24 hours
      id == "40170" ~ ymd_hm("2025-10-09 14:31", tz = "America/New_York"),
      id == "40197" ~ ymd_hm("2025-09-19 15:19", tz = "America/New_York"),
      id == "40199" ~ ymd_hm("2025-09-18 12:10", tz = "America/New_York"),
      TRUE ~ as.POSIXct(NA, tz = "America/New_York")
    ),
    
    start_time = with_tz(tag_time + hours(24), "UTC") #convert to UTC to match tag time and then filter
  ) %>%
  filter(date >= start_time) %>%
  select(-tag_time, -start_time)


## bathy and land maps ----
# depth
bathy <- marmap::getNOAA.bathy(-85, -5, 25, 50, resolution = 1, keep = T)
bathy_df <- marmap::fortify.bathy(bathy)

# land
land <- ne_countries(scale = "large", returnclass = "sf")


# Split tracks into segments and format for aniMotum ----
source("~/Documents/R Functions/segment_tracks.R") #load function
#arguments: gap = max # of days to interpolate across before splitting; ping_threshold = number of pings within each segment needed to try to fit
res <- segment_tracks(df,gap=10, ping_threshold = 10) #max gap 10 days, #minimum 10 pings per segment

fulldat <- res$fulldat #all data
filtered <- res$filtered #segmented, with only passing segments
summ_segments <- res$summ_segments #segment by segment summaries
id_summary <- res$id_summary #individual summaries

#animotum format
filtered <- format_data(filtered, id='id',
                        date='date', #date
                        lc = 'lc', #argos location class
                        coord=c("Longitude","Latitude"),
                        epar=c('smaj','smin', 'eor'))

# MM32 ####
MM32 <- filtered %>%
  filter(
    id == "40170-1")

## fit MP model ----
fit.MM32 <- fit_ssm(MM32, model='mp', time.step=24, vmax=4,
                  map = list(psi = factor(NA)))
summary(fit.MM32)

#plot model
plot(fit.MM32, what='predicted') 
plot(fit.MM32, what='fitted') 
map(fit.MM32, what='predicted') 
map(fit.MM32, what='fitted') 

## One-Step-Ahead Residuals ----
res.MM32 <- osar(fit.MM32)
plot(res.MM32, "ts") # Plot time-series of prediction residuals 
plot(res.MM32, "qq") # qq plot
plot(res.MM32, "acf") # Plot autocorrelation functions of prediction residuals 

## Extract track and MP estimates ----
dailyMM32 <- grab(fit.MM32, what = "predicted") %>% 
  mutate(
    g_scaled = scales::rescale(g, to = c(0, 1), na.rm = TRUE) #scale for plotting
  )

## ROM summaries ----
dailyMM32 <- dailyMM32 %>%
  arrange(id, date) %>%
  group_by(id) %>%
  mutate(
    step_dist_km = sqrt((x - lag(x))^2 + (y - lag(y))^2),
    days_elapsed = as.numeric(difftime(date, lag(date), units = "days")),
    disp_rate_km_day = step_dist_km / days_elapsed,
    disp_rate_3d = rollmean(
      disp_rate_km_day,
      k = 3,
      fill = NA,
      align = "right"
    ),
    disp_rate_7d = rollmean(
      disp_rate_km_day,
      k = 7,
      fill = NA,
      align = "right"
    )
  ) %>%
  ungroup()

#get depth
dailyMM32$depth <- get.depth(bathy, x=dailyMM32$lon, y=dailyMM32$lat, locator = F)$depth
dailyMM32$depth = dailyMM32$depth*-1

# get distance to land
pts <- st_as_sf(dailyMM32, #get points
                coords = c("lon", "lat"),
                crs = 4326)

# project to an equal-area / metric CRS
pts_proj  <- st_transform(pts, 6933)
land_proj <- st_transform(land, 6933)

# distance to nearest land (meters)
dist_m <- st_distance(pts_proj, land_proj)

dailyMM32$dist_to_coast_km <- apply(dist_m, 1, min) / 1000 #merge back to data

#save dataframe for modeling
write_csv(dailyMM32, "MP_track_MM32.csv")


# MM31 ####
MM31 <- filtered %>%
  filter(
    id == "40197-1")

## fit MP model ----
fit.MM31 <- fit_ssm(MM31, model='mp', time.step=24, vmax=4,
                    map = list(psi = factor(NA)))
summary(fit.MM31)

#plot model
plot(fit.MM31, what='predicted') 
plot(fit.MM31, what='fitted') 
map(fit.MM31, what='predicted') 
map(fit.MM31, what='fitted') 

## One-Step-Ahead Residuals ----
res.MM31 <- osar(fit.MM31)
plot(res.MM31, "ts") # Plot time-series of prediction residuals 
plot(res.MM31, "qq") # qq plot
plot(res.MM31, "acf") # Plot autocorrelation functions of prediction residuals 

## Extract track and MP estimates ----
dailyMM31 <- grab(fit.MM31, what = "predicted") %>% 
  mutate(
    g_scaled = scales::rescale(g, to = c(0, 1), na.rm = TRUE) #scale for plotting
  )

## ROM summaries ----
dailyMM31 <- dailyMM31 %>%
  arrange(id, date) %>%
  group_by(id) %>%
  mutate(
    step_dist_km = sqrt((x - lag(x))^2 + (y - lag(y))^2),
    days_elapsed = as.numeric(difftime(date, lag(date), units = "days")),
    disp_rate_km_day = step_dist_km / days_elapsed,
    disp_rate_3d = rollmean(
      disp_rate_km_day,
      k = 3,
      fill = NA,
      align = "right"
    ),
    disp_rate_7d = rollmean(
      disp_rate_km_day,
      k = 7,
      fill = NA,
      align = "right"
    )
  ) %>%
  ungroup()

#get depth
dailyMM31$depth <- get.depth(bathy, x=dailyMM31$lon, y=dailyMM31$lat, locator = F)$depth
dailyMM31$depth = dailyMM31$depth*-1

# get distance to land
pts <- st_as_sf(dailyMM31, #get points
                coords = c("lon", "lat"),
                crs = 4326)

# project to an equal-area / metric CRS
pts_proj  <- st_transform(pts, 6933)
land_proj <- st_transform(land, 6933)

# distance to nearest land (meters)
dist_m <- st_distance(pts_proj, land_proj)

dailyMM31$dist_to_coast_km <- apply(dist_m, 1, min) / 1000 #merge back to data

#save dataframe for modeling
write_csv(dailyMM31, "MP_track_MM31.csv")


#MM29 ####
#isolate one fish and remove 24 hrs after tagging
#need to fit first segment
MM29 <- filtered %>%
  filter(
    id == "40199-2")

fit.MM29 <- fit_ssm(MM29, model='mp', time.step=24, vmax=4,
                    map = list(psi = factor(NA)))
summary(fit.MM29)
plot(fit.MM29)
map(fit.MM29)

## One-Step-Ahead Residuals ----
res.MM29 <- osar(fit.MM29)
plot(res.MM29, "ts") # Plot time-series of prediction residuals 
plot(res.MM29, "qq") # qq plot
plot(res.MM29, "acf") # Plot autocorrelation functions of prediction residuals 

#track part 1
MM29_early <- filtered %>%
  filter(
    id == "40199-1")
MM29_early <- MM29_early[-1,]

fit.MM29_early <- fit_ssm(MM29_early, model='mp', time.step=24, vmax=4,
                    map = list(psi = factor(NA)))

summary(fit.MM29_early)
plot(fit.MM29_early)
map(fit.MM29_early)

dailyMM29 <- grab(fit.MM29, what='predicted') 
dailyMM29_early <- grab(fit.MM29_early, what='predicted') 

dailyMM29 <- dailyMM29 %>%
  bind_rows(dailyMM29_early) %>%
  mutate(id = '40199',
         g_scaled = scales::rescale(g, to = c(0, 1), na.rm = TRUE)) 

dailyMM29 <- dailyMM29 %>%
  arrange(id, date) %>%
  group_by(id) %>%
  mutate(
    step_dist_km = sqrt((x - lag(x))^2 + (y - lag(y))^2),
    days_elapsed = as.numeric(difftime(date, lag(date), units = "days")),
    disp_rate_km_day = step_dist_km / days_elapsed,
    disp_rate_3d = rollmean(
      disp_rate_km_day,
      k = 3,
      fill = NA,
      align = "right"
    ),
    disp_rate_7d = rollmean(
      disp_rate_km_day,
      k = 7,
      fill = NA,
      align = "right"
    )
  ) %>%
  ungroup()

#get depth
dailyMM29$depth <- get.depth(bathy, x=dailyMM29$lon, y=dailyMM29$lat, locator = F)$depth
dailyMM29$depth = dailyMM29$depth*-1

# distance to land
pts <- st_as_sf(dailyMM29,
                coords = c("lon", "lat"),
                crs = 4326)

# project to an equal-area / metric CRS
pts_proj  <- st_transform(pts, 6933)
land_proj <- st_transform(land, 6933)

# distance to nearest land (meters)
dist_m <- st_distance(pts_proj, land_proj)

dailyMM29$dist_to_coast_km <- apply(dist_m, 1, min) / 1000

#save dataframe for modeling
write_csv(dailyMM29, "MP_track_MM29.csv")
