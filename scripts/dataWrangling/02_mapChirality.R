library(terra)
library(tidyverse)
library(sf)
library(glmmTMB)


r <- rast("predictors/oahu.bil")
plot(r)
r
crs(r) <- "+init=epsg:26704"

preds <- rast("predictors/baseline_bioclims.tif")
preds
elev <- rast("predictors/elevation.tif")
elev # elevation layer is at coarser resolution than pred vars, try to get finer HI elevation raster layer
elev <- crop(elev, preds)

bio1 <- preds[[1]]
bio12 <- preds[[12]]
elev <- r %>% 
    project(bio1)

bbox <- st_bbox(c(xmin = -158.4, xmax = -157.5, ymin = 21.2, ymax = 21.8))

tempRainElev <- c(bio1,bio12,elev) %>% 
    crop(bbox)

plot(tempRainElev)

## richness per grid cell 
# add richness to model # 
rm_files <- list.files(path = "output_polygons/", pattern = "*.shp", full.names = T)
rasterize_toCSV <- function(fz){
    temp_vec <-  vect(fz) %>% 
        project(tempRainElev)
    temp_rast <- rasterize(temp_vec, tempRainElev) 
    temp_df <- terra::as.data.frame(temp_rast, xy = TRUE)
    temp_df <- temp_df %>%            
        mutate(binomial = word(fz, start = 2, end = 3, sep = fixed("_")) %>% str_remove("polygons/"))
    return(temp_df)
}

rm <- lapply(X = rm_files, FUN = rasterize_toCSV)
rm2 <- bind_rows(rm)

rich_df <- rm2 %>% 
    group_by(x,y) %>% 
    summarise(richness = length(unique(binomial))) %>% 
    ungroup()


ggplot(rich_df, mapping = aes(x = x, y = y, fill = richness)) +
    geom_tile()

tre_df <- as.data.frame(tempRainElev, xy = TRUE)

preds_df <- left_join(tre_df, rich_df)

preds <- preds_df %>% 
    rename(bio1 = baseline_bioclims_1,
           bio12 = baseline_bioclims_12,
           elev = oahu) %>% 
    select(x,y,richness,bio1,bio12,elev) %>% 
    rast()
plot(preds)

#read in chiral data
cd <- read.csv("data/Chirality_data-for-mapping-fromChandra.csv")
head(cd)

#read in HI shapefile
hi <- vect("shapefiles/HI.shp") %>% 
    st_as_sf()

ggplot() +
    geom_sf(hi, mapping = aes()) +
    geom_point(cd, mapping = aes(x = decimallongitudeUSE, y = decimallatitudeUSE)) +
    coord_sf(xlim = c(-158.4,-157.5),
             ylim = c(21.2,21.8)) +
    theme_classic()

# now let's plot chirality by species, make into one pdf
# then we'll show average chirality across all species

# let's restrict the points and shapefile to only oahu
cd_oahu <- cd %>% 
    filter(decimallongitudeUSE > -158.4 & decimallongitudeUSE < -157.5,
           decimallatitudeUSE > 21.2 & decimallatitudeUSE < 21.8)

bbox <- st_bbox(c(xmin = -158.4, xmax = -157.5, ymin = 21.2, ymax = 21.8))
oahu <- hi %>% 
    st_crop(bbox)


ggplot() +
    geom_sf(oahu, mapping = aes(), fill = NA) +
    geom_point(cd_oahu, 
               mapping = aes(x = decimallongitudeUSE, y = decimallatitudeUSE,
                             color = X.dex)) +
    scale_color_viridis_c() + 
    theme_classic()

ggsave(filename = "figures/allPoints_dext.png", dpi = 450, width = 8.5, height = 11)

ggplot() +
    geom_sf(oahu, mapping = aes(), fill = NA) +
    geom_point(cd_oahu, 
               mapping = aes(x = decimallongitudeUSE, y = decimallatitudeUSE,
                             color = X.dex)) +
    scale_color_viridis_c() + 
    theme_classic() +
    facet_wrap(~correctedspeciesname, ncol = 1)

ggsave(filename = "figures/species_dextMap.pdf", 
       width = 8.5, height = 90, limitsize =  FALSE)

# transform data to not be percentage
cd_oahu_ind <- cd_oahu %>% 
    mutate(dex = if_else(condition = X.sin == 1 & is.na(dex), true = 0, false = dex)) %>% 
    filter(!is.na(n),
           !is.na(sin),
           !is.na(dex)) %>% 
    select(correctedspeciesname, decimallatitudeUSE, decimallongitudeUSE, coordinateuncertaintyinMetersUSE, n, sin) %>% 
    mutate(dex = n-sin)
dex_df <- uncount(cd_oahu_ind, dex) %>% 
    mutate(dex = 1) %>% 
    select(correctedspeciesname, decimallatitudeUSE, decimallongitudeUSE,dex)
sin_df <- uncount(cd_oahu_ind, sin) %>% 
    mutate(dex = 0) %>% 
    select(correctedspeciesname, decimallatitudeUSE, decimallongitudeUSE,dex)

mdf <- bind_rows(dex_df, sin_df)

mdf2 <- st_as_sf(mdf, 
                 coords = c("decimallongitudeUSE", "decimallatitudeUSE"),
                 crs = st_crs(tempRainElev)) %>% 
    vect()

e <- terra::extract(preds, mdf2)

mdf <- cbind(mdf, e)
mdf <- mdf %>% 
    select(-ID) %>% 
    mutate(Species = word(correctedspeciesname,2,2))

## add body size to model
traits <- read.csv("data/auriculella_summary_all_characters_use_raw.csv") %>% 
    select(Species, meanWeight_multipleIndv, Ratio)

mdf <- mdf %>% 
    left_join(traits)

mdf_scaled <- mdf %>% 
    #  mutate(richness = if_else(is.na(richness), true = 0, false = richness)) %>% 
    mutate(bio1_sc = scale(bio1),
           bio12_sc = scale(bio12),
           elev_sc = scale(elev),
           meanMass_sc = scale(meanWeight_multipleIndv),
           ratio_sc = scale(Ratio),
           richness_sc = scale(richness)) 

mdf_scaled <- na.omit(mdf_scaled)


m1 <- glmmTMB(formula = dex ~ bio1_sc + bio12_sc + elev_sc + richness_sc + meanMass_sc + ratio_sc + 
                  (1 | Species) +
                  (0 + bio1_sc|Species) + 
                  (0 + bio12_sc|Species) + 
                  (0 + elev_sc|Species),
              data = mdf_scaled,
              family = binomial)

summary(m1)   

## run with lasso! - regularisation 

# USe cross validation to select value of regulator (which governs parametrization)

library(DHARMa)
so <- simulateResiduals(m1)
DHARMa::plotSimulatedResiduals(so)

# predict model
pred_df <- terra::as.data.frame(preds, xy = TRUE) %>% 
    filter(!is.na(bio1)) %>% 
    mutate(correctedspeciesname = NA,
           ratio_sc = 0,
           meanMass_sc = 0) %>% 
    # mutate(richness = if_else(is.na(richness), true = 0, false = richness)) %>% 
    rename(decimallatitudeUSE = y, 
           decimallongitudeUSE = x)

pred_df_sc <- pred_df %>% 
    mutate(richness_sc = scale(richness, center = mean(mdf$richness, na.rm = T), scale = sd(mdf$richness, na.rm = T)),
           bio1_sc = scale(bio1, center = mean(mdf$bio1, na.rm = T), scale = sd(mdf$bio1, na.rm = T)),
           bio12_sc = scale(bio1, center = mean(mdf$bio12, na.rm = T), scale = sd(mdf$bio12, na.rm = T)),
           elev_sc = scale(elev, center = mean(mdf$elev, na.rm = T), scale = sd(mdf$elev, na.rm = T)),
           Species = NA
    )


pOut <- predict(object = m1, newdata = pred_df_sc, type = "response")
pOut_df <- pred_df %>% 
    mutate(estimate = pOut)

ggplot() +
    geom_tile(data = pOut_df, mapping = aes(x = decimallongitudeUSE, y = decimallatitudeUSE, fill = estimate)) +
    scale_fill_viridis_c() +
    theme_classic()


## predict species-specific stuff 
binomial <- "Auriculella auricula"
s <- word(binomial, 2,2)
rsc <- filter(mdf_scaled, correctedspeciesname == binomial)$ratio_sc[1]
rsc <- filter(mdf_scaled, correctedspeciesname == binomial)$ratio_sc[1]
mmsc <- filter(mdf_scaled, correctedspeciesname == binomial)$meanMass_sc[1]
spp_poly <- vect(paste0("output_polygons/", str_replace(binomial, pattern = " ", "_"),
                        "_distribution.shp")) %>% 
    project(tempRainElev)

pred_spp <- mask(preds, spp_poly)

pred_df <- terra::as.data.frame(pred_spp, xy = TRUE) %>% 
    filter(!is.na(bio1)) %>% 
    mutate(Species = s,
           ratio_sc = rsc,
           meanMass_sc = mmsc) %>% 
    # mutate(richness = if_else(is.na(richness), true = 0, false = richness)) %>% 
    rename(decimallatitudeUSE = y, 
           decimallongitudeUSE = x)

pred_df_sc <- pred_df %>% 
    mutate(richness_sc = scale(richness, center = mean(mdf$richness, na.rm = T), scale = sd(mdf$richness, na.rm = T)),
           bio1_sc = scale(bio1, center = mean(mdf$bio1, na.rm = T), scale = sd(mdf$bio1, na.rm = T)),
           bio12_sc = scale(bio1, center = mean(mdf$bio12, na.rm = T), scale = sd(mdf$bio12, na.rm = T)),
           elev_sc = scale(elev, center = mean(mdf$elev, na.rm = T), scale = sd(mdf$elev, na.rm = T))
    )


pOut <- predict(object = m1, newdata = pred_df_sc, type = "response")
pOut_df <- pred_df %>% 
    mutate(estimate = pOut) 

ggplot() +
    geom_tile(data = pOut_df, mapping = aes(x = decimallongitudeUSE, y = decimallatitudeUSE, fill = estimate)) +
    scale_fill_viridis_c() +
    theme_classic()




# function to predict for each species
predFun <- function(binomial){
    
    s <- word(binomial, 2,2) 
    rsc <- filter(mdf_scaled, correctedspeciesname == binomial)$ratio_sc[1]
    mmsc <- filter(mdf_scaled, correctedspeciesname == binomial)$meanMass_sc[1]
    spp_poly <- vect(paste0("output_polygons/", str_replace(binomial, pattern = " ", "_"),
                            "_distribution.shp")) %>% 
        project(tempRainElev)
    pred_spp <- mask(preds, spp_poly)
    
    pred_df <- terra::as.data.frame(pred_spp, xy = TRUE) %>% 
        filter(!is.na(bio1)) %>% 
        mutate(Species = s,
               ratio_sc = rsc,
               meanMass_sc = mmsc) %>% 
        # mutate(richness = if_else(is.na(richness), true = 0, false = richness)) %>% 
        rename(decimallatitudeUSE = y, 
               decimallongitudeUSE = x)
    
    pred_df_sc <- pred_df %>% 
        mutate(richness_sc = scale(richness, center = mean(mdf$richness, na.rm = T), scale = sd(mdf$richness, na.rm = T)),
               bio1_sc = scale(bio1, center = mean(mdf$bio1, na.rm = T), scale = sd(mdf$bio1, na.rm = T)),
               bio12_sc = scale(bio1, center = mean(mdf$bio12, na.rm = T), scale = sd(mdf$bio12, na.rm = T)),
               elev_sc = scale(elev, center = mean(mdf$elev, na.rm = T), scale = sd(mdf$elev, na.rm = T))
        )
    
    
    pOut <- predict(object = m1, newdata = pred_df_sc, type = "response")
    pOut_df <- pred_df %>% 
        mutate(estimate = pOut) %>% 
        mutate(binomial = binomial)
    
    return(pOut_df)
}



pOut_df <- predFun(binomial = binomial)

ggplot() +
    geom_tile(data = pOut_df, mapping = aes(x = decimallongitudeUSE, y = decimallatitudeUSE, fill = estimate)) +
    scale_fill_viridis_c() +
    theme_classic()

pOut_df <- predFun(binomial = "Auriculella tenuis")

ggplot() +
    geom_tile(data = pOut_df, mapping = aes(x = decimallongitudeUSE, y = decimallatitudeUSE, fill = estimate)) +
    scale_fill_viridis_c() +
    theme_classic()

spp_list <- data.frame(spp = unique(cd_oahu$correctedspeciesname))

spp_list <- filter(spp_list, !spp %in% c("Auriculella cerea",
                                         "Auriculella crassula",
                                         "Auriculella armata"))$spp

lOut <- lapply(X = spp_list, predFun)

sppPOut_df <- bind_rows(lOut)

ggplot() +
    geom_sf(oahu, mapping = aes(), fill = NA) +
    geom_tile(data = sppPOut_df, 
              mapping = aes(x = decimallongitudeUSE, y = decimallatitudeUSE, fill = estimate)) +
    scale_fill_viridis_c(na.value = "transparent") +
    theme_classic() +
    facet_wrap(~binomial, ncol = 1)
