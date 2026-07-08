library(spdep)
library(sf)
library(brms)
library(tidyverse)

step1 <- readRDS("brmsFittedModels/manuscript/step1.rds")

# Extract fitted probabilities and compute raw residuals manually
fitted_vals <- fitted(step1)[, "Estimate"]
observed    <- step1$data$is_monomorphic

resid_step1 <- observed - fitted_vals

# Average residuals per site (since multiple lots per site)
## generate site IDs
## read in data
mdf_scaled <- read.csv("data/mdf_scaled_July2026Update.csv") %>% 
    mutate(n_sc = scale(n),
           mpd_sc = scale(MPD)) %>% 
    filter(correctedspeciesname != "Auriculella cerea")

uniqueSites <- distinct(mdf_scaled, decimallatitudeUSE, decimallongitudeUSE)

uniqueSites <- uniqueSites %>% 
    mutate(siteID = 1:nrow(uniqueSites))

resid_by_site <- step1$data %>%
    left_join(uniqueSites) %>% 
    mutate(resid = resid_step1) %>%
    group_by(siteID, decimallongitudeUSE, decimallatitudeUSE) %>%
    summarise(mean_resid = mean(resid, na.rm = TRUE),
              .groups = "drop")

# Build spatial weights from coordinates
library(spdep)

coords <- cbind(resid_by_site$decimallongitudeUSE,
                resid_by_site$decimallatitudeUSE)

nb <- knn2nb(knearneigh(coords, k = 5))
lw <- nb2listw(nb, style = "W")

# Run Moran's I
moran_result <- moran.test(resid_by_site$mean_resid, lw)
moran_result
# Build spatial weights from your site coordinates
coords <- step1$data %>% 
    left_join(uniqueSites) %>% 
    distinct(siteID, decimallongitudeUSE, decimallatitudeUSE) %>%
    arrange(siteID)

# Average residuals per site
resid_by_site <- step1$data %>% 
    left_join(uniqueSites) %>% 
    mutate(resid = resid_step1) %>%
    group_by(siteID, decimallongitudeUSE, decimallatitudeUSE) %>%
    summarise(mean_resid = mean(resid), .groups = "drop")

# Moran's I test
nb <- knn2nb(knearneigh(
    cbind(resid_by_site$decimallongitudeUSE, resid_by_site$decimallatitudeUSE), 
    k = 5
))
lw <- nb2listw(nb)
moran.test(resid_by_site$mean_resid, lw)

## now for step 2
step2 <- readRDS("brmsFittedModels/manuscript/step2.rds")

# Extract fitted probabilities and compute raw residuals manually
fitted_vals <- fitted(step2)[, "Estimate"]
observed    <- step2$data$dex

resid_step2 <- observed - fitted_vals

# Average residuals per site (since multiple lots per site)
## generate site IDs
uniqueSites <- distinct(mdf_scaled, decimallatitudeUSE, decimallongitudeUSE)

uniqueSites <- uniqueSites %>% 
    mutate(siteID = 1:nrow(uniqueSites))

resid_by_site <- step2$data %>%
    left_join(uniqueSites) %>% 
    mutate(resid = resid_step2) %>%
    group_by(siteID, decimallongitudeUSE, decimallatitudeUSE) %>%
    summarise(mean_resid = mean(resid, na.rm = TRUE),
              .groups = "drop")

# Build spatial weights from coordinates
library(spdep)

coords <- cbind(resid_by_site$decimallongitudeUSE,
                resid_by_site$decimallatitudeUSE)

nb <- knn2nb(knearneigh(coords, k = 5))
lw <- nb2listw(nb, style = "W")

# Run Moran's I
moran_result <- moran.test(resid_by_site$mean_resid, lw)
moran_result

