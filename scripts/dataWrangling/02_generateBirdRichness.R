library(terra)
library(sf)
library(tidyverse)

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

## read in bird df
birds_sf <- vect("shapefiles/data_0.shp") %>% 
    st_as_sf()

# crop to oahu
#read in HI shapefile
hi <- vect("shapefiles/HI.shp") %>% 
    st_as_sf()

bbox <- st_bbox(c(xmin = -158.4, xmax = -157.5, ymin = 21.2, ymax = 21.8))
oahu <- hi %>% 
    st_crop(bbox)

birds_sf <- st_crop(birds_sf, oahu) 

# remove introduced birds
unique(birds_sf$LEGEND)

birds_sf <- birds_sf %>% 
    filter(!LEGEND %in% c("Extant & Introduced (resident)",
                          "Extant & Introduced (seasonality uncertain)"))

## get higher taxonomy of remaining birds
library(taxotools)
mylist <- data.frame("canonical" = birds_sf$SCI_NAME)
my_taxo_list <- taxotools::list_higher_taxo(mylist,"canonical")
my_taxo_list <- my_taxo_list %>% 
    mutate(Family = case_when(
        canonical == "Hydrobates castro" ~ "Hydrobatidae",
        canonical == "Hydrobates tristrami" ~ "Hydrobatidae",
        canonical == "Gygis candida" ~ "Laridae",
        canonical == "Gygis candida" ~ "Laridae",
        canonical == "Myadestes woahensis" ~ "Turdidae",
        .default = Family
    )) %>% 
    distinct(.keep_all = TRUE)
unique(my_taxo_list$Family)
keepFamilies <- c("Fringillidae", "Monarchidae", "Turdidae")

birds_sf <- left_join(birds_sf,my_taxo_list, by = c("SCI_NAME"="canonical")) %>% 
    filter(Family %in% keepFamilies)

bird_sf_tableS1 <- birds_sf %>% 
    st_drop_geometry() %>% 
    rename(yearCompiled = YRCOMPILED,
           citation = CITATION,
           status = LEGEND,
           scientificName = SCI_NAME) %>% 
    select(scientificName, yearCompiled, status, Family)

write.csv(x = bird_sf_tableS1, file = 'tables/manuscript/birdS1.csv', row.names = F)

# richnesss per grid
birds <- unique(birds_sf$SCI_NAME)
x <- birds[2]
birdsToCSV <- function(bn){
    
    temp_vec <- filter(birds_sf, SCI_NAME == bn)
    temp_vec <-  vect(temp_vec) %>% 
        project(tempRainElev)
    temp_rast <- rasterize(temp_vec, tempRainElev) 
    temp_df <- terra::as.data.frame(temp_rast, xy = TRUE)
    temp_df <- temp_df %>%            
        mutate(binomial = bn)
    return(temp_df)
}

bird_df <- lapply(X = birds, FUN = birdsToCSV)
bird_df <- bind_rows(bird_df)

bird_rich_df <- bird_df %>% 
    group_by(x,y) %>% 
    summarise(richness = length(unique(binomial))) %>% 
    ungroup()

write.csv(x = bird_rich_df, file = "data/birdSpeciesRichness.csv", row.names = F)

bird_rich_rast <- rast(bird_rich_df)

plot(bird_rich_rast)

writeRaster(bird_rich_rast, filename = "data/birdSpeciesRichness.tif")
