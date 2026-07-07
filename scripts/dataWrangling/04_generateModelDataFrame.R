library(terra)
library(tidyverse)
library(sf)
library(brms)

r <- rast("predictors/oahu.bil")
crs(r) <- "+init=epsg:26704"

preds <- rast("predictors/baseline_bioclims.tif")
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

tre_df <- as.data.frame(tempRainElev, xy = TRUE)

preds_df <- left_join(tre_df, rich_df)
preds_df <- preds_df %>% 
    rename(snailRichness = richness)

preds <- preds_df %>% 
    rename(bio1 = baseline_bioclims_1,
           bio12 = baseline_bioclims_12,
           elev = oahu) %>% 
    select(x,y,snailRichness,bio1,bio12,elev) %>% 
    rast()

#read in chiral data
cd <- read.csv("data/Chirality_data-for-mapping-fromChandra.csv")
head(cd)
sort(unique(cd$correctedspeciesname))

#read in HI shapefile
hi <- vect("shapefiles/HI.shp") %>% 
    st_as_sf()

# let's restrict the points and shapefile to only oahu
cd_oahu <- cd %>% 
    filter(decimallongitudeUSE > -158.4 & decimallongitudeUSE < -157.5,
           decimallatitudeUSE > 21.2 & decimallatitudeUSE < 21.8)

bbox <- st_bbox(c(xmin = -158.4, xmax = -157.5, ymin = 21.2, ymax = 21.8))
oahu <- hi %>% 
    st_crop(bbox)

## model data frame
mdf <- cd_oahu
sort(unique(mdf$correctedspeciesname))

mdf2 <- st_as_sf(mdf, 
                 coords = c("decimallongitudeUSE", "decimallatitudeUSE"),
                 crs = st_crs(tempRainElev)) %>% 
    vect()

e <- terra::extract(preds, mdf2, touches = TRUE)

mdf <- cbind(mdf, e)
mdf <- mdf %>% 
    select(-ID) %>% 
    mutate(Species = word(correctedspeciesname,2,2))

## also grab richness of birds
bird_rich <- rast("data/birdSpeciesRichness.tif")
eBird <- terra::extract(bird_rich, mdf2, touches = TRUE)
mdf <- mdf %>% 
    mutate(birdRichness = eBird$richness)

# okay i think it is fine to replace NAs with 1s
mdf <- mdf %>% 
    mutate(snailRichness = if_else(
        is.na(snailRichness), true = 1, false = snailRichness
    ))

## add body size to model
traits <- read.csv('data/shell_data.csv') %>%
    mutate(shellDim = width * width * height)

ggplot(traits) +
    geom_point(mapping = aes(x = shellDim, y = mass))

ggplot(traits) +
    geom_point(mapping = aes(x = shellDim, y = log_mass))


ggplot(traits) +
    geom_point(mapping = aes(x = log(shellDim), y = log(mass)))

t_lm <- lm(log_mass ~ shellDim, data = traits)
t_resids <- residuals(t_lm)

traits <- mutate(traits, 
                 relativeMass = t_resids)

## now get residuals of ratio to shell dim
ggplot(traits) +
    geom_point(mapping = aes(x = ratio, y = shellDim))

ratio_lm <- lm(ratio ~ shellDim, data = traits)
ratio_resids <- residuals(ratio_lm)
traits <- mutate(traits, 
                 relativeRatio = ratio_resids)

#westerlundia is armata
#spaldingi is gagneorum 

sort(unique(traits$Species))

traits <- traits %>% 
    mutate(species = case_when(
        species == "A. spaldingi" ~ "A. gagneorum",
        species == "A. westerlundia" ~ "A. armata", .default = species
    ))

mdf <- mdf %>% 
    mutate(species = paste("A.", Species))

mdf <- mdf %>% 
    left_join(traits, by = "species")

## what species don't have traits? 
noTraits <- filter(mdf, is.na(mass))
unique(noTraits$correctedspeciesname)


mdf_scaled <- mdf %>% 
    mutate(bio1_sc = scale(bio1),
           bio12_sc = scale(bio12),
           elev_sc = scale(elev),
           meanMass_sc = scale(mass),
           ratio_sc = scale(ratio),
           snailRichness_sc = scale(snailRichness),
           birdRichness_sc = scale(birdRichness),
           meanHeight_sc = scale(height),
           meanWidth_sc = scale(width),
           shellDim_sc = scale(shellDim),
           relativeMass_sc = scale(relativeMass),
           relativeRatio_sc = scale(relativeRatio)) 

sort(unique(mdf_scaled$correctedspeciesname))

mdf_scaled <- mdf_scaled %>% 
    filter(!is.na(bio1_sc),
           !is.na(bio12_sc),
           !is.na(elev_sc),
           !is.na(meanMass_sc),
           !is.na(ratio_sc),
           !is.na(snailRichness_sc),
           !is.na(birdRichness_sc),
           !is.na(meanHeight_sc),
           !is.na(meanWidth_sc)) # not so many lost

sort(unique(mdf_scaled$correctedspeciesname))

## remove the following species 
#Remove amata observations
mdf_scaled <- filter(mdf_scaled, correctedspeciesname != "Auriculella armata")

#Remove crassula observations
mdf_scaled <- filter(mdf_scaled, correctedspeciesname != "Auriculella crassula")

write.csv(mdf_scaled, file = "data/mdf_scaled_March2026Update.csv")

