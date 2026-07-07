library(tidyverse)
library(terra)
library(sf)

rm_files <- list.files(path = "sdmOutput_polygons_oahuOnly/", pattern = "*.shp", full.names = T)

list_of_spatvectors <- lapply(rm_files, vect)
combined_spatvector <- do.call(rbind, list_of_spatvectors)

plot(combined_spatvector)

# restrict to only oahu
#read in HI shapefile
hi <- vect("shapefiles/HI.shp") %>% 
    st_as_sf()

# let's restrict the points and shapefile to only oahu
bbox <- st_bbox(c(xmin = -158.4, xmax = -157.5, ymin = 21.2, ymax = 21.8))
oahu <- hi %>% 
    st_crop(bbox) %>% 
    st_transform(st_crs(combined_spatvector))

combined_spatvector <- crop(combined_spatvector, oahu)
plot(combined_spatvector)

## read in data for models 
mdf_scaled <- read.csv("data/mdf_scaled_July2026Update.csv")

## wrangle to picante style data
mdf3 <- mdf_scaled %>% 
    group_by(correctedspeciesname, decimallatitudeUSE, decimallongitudeUSE) %>% 
    summarise(count = sum(n)) %>% 
    ungroup()

uniqueSites <- distinct(mdf3, decimallatitudeUSE, decimallongitudeUSE)

uniqueSites <- uniqueSites %>% 
    mutate(siteID = 1:nrow(uniqueSites))

uniqueSites_sf <- st_as_sf(uniqueSites, 
                           coords = c("decimallongitudeUSE",
                                      "decimallatitudeUSE"),
                           crs = "WGS84") %>% 
    st_transform(crs = st_crs(combined_spatvector))

uniqueSites_v <- vect(uniqueSites_sf)
plot(uniqueSites_v)

#' write funciton to do this one shapefile at a time
#' 
#' 
sitePA <- function(x){
    sdm <- vect(rm_files[x])
    sp <- word(sources(sdm),start = 8, end = 8, sep = fixed("\\")) %>% 
        str_remove("_distribution.shp")
    
    sdm_crop <- crop(sdm, oahu)
    
    e <- terra::extract(sdm_crop, uniqueSites_v) %>% 
        dplyr::rename(siteID = 1) %>% 
        mutate(PA = if_else(condition = is.na(Ar__SDM), 0, 1)) %>% 
        mutate(scientificName = sp) %>% 
        select(siteID, PA, scientificName)
    
    return(e)
}

PAList <- lapply(X = seq_along(rm_files), FUN = sitePA)

combinedExtracts <- do.call(rbind, PAList)

# uniqueSites PA
p <- pivot_wider(combinedExtracts, names_from = scientificName,
                 values_from = PA)
p2 <- p %>% tibble::column_to_rownames(var = "siteID")

mdf_scaled2 <- mdf_scaled %>% 
    mutate(scientificName = str_replace_all(correctedspeciesname, " ", "_")) 

# read in phylogeny
library(ape)
library(picante)
tree <- read.tree(file = "data/auriculella_species_collapsed.tre")
tree$tip.label

tree <- prune.sample(p2, tree)

comm <-p[,tree$tip.label]
comm

# get vectors to join to data.frames
Site <- c(rownames(p2), "totalPhy")

#PD
pd <- pd(comm, tree, include.root = TRUE)
pd
# what's the minimum PD above 0?
filter(pd, PD > 0) %>% min()
pd <- pd %>% 
    mutate(PD = if_else(condition = PD == 0,
                       true = 0.007, 
                       false = PD))
hist(pd$PD)

# now calculate mpd
# MPD modeling
Site <- c(rownames(p2))
phy.dist<-cophenetic(tree)
phydist <- as.matrix(phy.dist)
comm.m <- as.matrix(comm)

mpd <-mpd(comm.m, phy.dist, abundance.weighted = FALSE)
mpd

## join outputs with uniqueSites
uniqueSites_mdf <- uniqueSites %>% 
    mutate(PD = pd$PD,
           SR = pd$SR,
           MPD = mpd
           )

mdf_scaled3 <- left_join(mdf_scaled2, uniqueSites_mdf)
write.csv(x = mdf_scaled3, file = "data/mdf_scaled_July2026Update.csv",
          row.names = F)
