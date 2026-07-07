library(magrittr)
library(dplyr)
library(terra)
library(ggplot2)

my_dir_path <- getwd()
wd <- list()
wd$R       <- file.path( my_dir_path, "R" )
wd$fun     <- file.path( my_dir_path, "functions" )
wd$bin     <- file.path( my_dir_path, "bin" )
wd$data    <- file.path( my_dir_path, "data" )
wd$occs    <- file.path( wd$data, "Better_Chirality_data-for-mapping_badpoints_removed.csv")
wd$figs    <- file.path( my_dir_path, "figs" )
wd$out     <- file.path( my_dir_path, "out" )
invisible({ lapply(wd, function(i) if( dir.exists(i) != 1 ) dir.create(i) ) })

proj4_aea <- "+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-157 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
proj4_wgs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
