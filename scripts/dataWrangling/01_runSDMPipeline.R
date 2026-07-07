# Set up pipeline ---------------------------------------------------------

## Step 1: Load in the necessary libraries ----
source('R/00_setup.R')
library(rnaturalearth)
library(rgeos)
library(dismo)
library(ENMeval)
library(stringr)
library(dplyr)
library(terra)
library(sf)
library(rJava)

## Step 2: Load in the pipeline scripts ----
list.files(wd$fun, full.names = T) %>% 
  lapply(source)

## Step 3: Prepare the spatial data ----
### 3C. load occurrence records
occs <- data.table::fread("data/Better_Chirality_data-for-mapping_badpoints_removed.csv") %>% 
  dplyr::mutate(decimalLongitude = decimallongitudeUSE,
                decimalLatitude = decimallatitudeUSE) %>% 
  dplyr::filter(!is.na(decimalLongitude),
                !is.na(decimalLatitude)) %>% 
  dplyr::mutate(bw = str_replace(correctedspeciesname, " ", "_")) #bw stands for "Binomial With _ between genus and specific epithet"

### 3D. Read in basemap for visualizing
world <- vect(x = "shapefiles/HI.shp") %>% st_as_sf()

### 3E. Choose projection and project data if necessary
study_proj <- proj4_aea #aea projection
world <- st_transform(world, crs = study_proj) 

occs_sf <- st_as_sf(occs, 
                    coords = c("decimalLongitude", "decimalLatitude"),
                    crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") 
occs_sf <- st_transform(occs_sf, crs = study_proj)
occs <- occs %>% 
  mutate(x = st_coordinates(occs_sf)[,1], y = st_coordinates(occs_sf)[,2])


############################ Part 3: Execute the pipeline #############################

rerunPipeline <- FALSE

#### Step 1: Set up species list
spp_list <- unique(occs$bw)

#### Step 2: Pass pipeline through the species_season list

sdm_pipeline <- function(x){
  
  print("=========================================================")
  print(paste("============== Working on", spp_list[x], "=============="))  
  print("=========================================================")
  binomial <- spp_list[x]
  species_df_raw <-   dplyr::filter(as.data.frame(occs), bw == binomial) %>% 
    filter(!is.na(x),
           !is.na(y))
  stopifnot(nrow(species_df_raw) >= 1)

  # Clip all potential predictor variables based on all cleaned occurrence
  # records from this species.
  
  ##### Change alphaBuff for each Aricuella species based on biological location factors
  # 1  ambusta    = 5000
  # 2  armata     = 6000 
  # 3  auricula    = 1000 
  # 4  brunnea    = 4000
  # 5  canalifera = 10000 
  # 6  castanea   = 4000 
  # 7  cerea      = 7000
  # 8  crassula   = 5000
  # 9  diaphana   = 6000
  # 10 expansa    = 4000
  # 11 flavida    = 5000 
  # 12 gagneorum  = 4000
  # 13 lanaiensis = 4000
  # 14 malleata   = 3000
  # 15 minuta     = 7000
  # 16 montana    = 5000
  # 17 newcombi   = 1000
  # 18 olivacea   = 3000
  # 19 perpusilla = 5000
  # 20 perversa   = 5000
  # 21 pulchra    = 6000
  # 22 serrula    = 6000
  # 23 straminea  = 5000
  # 24 tantalus   = 4000
  # 25 tenella    = 6000
  # 26 tenuis     = 8000
  # 27 turritella = 10000
  # 28 uniplicata = 7000
  aa_shp <- define_accessibleArea(species_df = species_df_raw, alphabuff = 1000, projCRS = study_proj)
  
  mod_vars <- rast("modelVariablesStack.tif")
  # Clip environmental variable layers to the defined accessible area
  mymod_vars <- clip_variableLayers(rstack = mod_vars, accessibleArea = aa_shp)

  # Thin points based on the accessible area.
  ### 2A. Prepare the coordinates for rarefaction
  coordinates(species_df_raw) <- ~ x + y
  
  ### 2B. Perform the rarefaction
#  area_sqkm <- raster::area(aa_shp)*0.000001 
  species_df <- thinPoints(
    spp_df = species_df_raw, 
    area_sqkm = area_sqkm, 
    bio = mymod_vars[[2]], 
    method = "simple",
    simpleMult = 4
  )
  spp_df <- species_df_raw
  
  #### Step 3: Test and fine-tune the model
  
  ### 3A. Select top performing variables to reduce colinearity
  ## First run a test model
  print(">>>>>>>>> Running Maxent test model <<<<<<<<<")
  max_model <- maxent(x = mymod_vars, p = coordinates(spp_df), progress = "text", silent = TRUE) 
  ## Using the test model, iteratively test the colinearity of variables, removing highly colinear ones one at a time 
  print(">>>>>>>>> Selecting top SDM variables <<<<<<<<<")
  predictors <- select_sdmVariables(pred_vars = mymod_vars, maxent_mod = max_model, maxVIF = 5)
  
  ### 3B. Evaluate various tuning parameters of the model for fine tuning to create the best performing model with best tuned parameters
  print(">>>>>>>>> Evaluating tuning variables in model <<<<<<<<<")
  
  eval1 <- ENMeval::ENMevaluate(
    occ = coordinates(spp_df), env = predictors,
    method = "block", RMvalues = c(0.5, 1, 2, 3, 4),
    fc= c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
    parallel = TRUE, numCores = 10, algorithm = 'maxent.jar',
    quiet = FALSE, updateProgress = TRUE
  )
  
  #### Step 4: Create the output for the best performing model
  
  ### 4A. Prepare the output path and coordinates
  bw <- spp_df$bw[1] # make a speciesName string without a space for better saving
  bw <- str_replace(bw, " ", "_")
  
  resultDir <-  wd$out
  if (!dir.exists(resultDir)) { dir.create(resultDir) }
  
  ### 4B. Output the model
  ## Output the model evaluation
  #  save(eval1, file = file.path(resultDir, paste( bw, "_ENMeval", ".RData", sep = "")))
  ## Output the best performing model, including both the actual model and the presence-absence model
  print(">>>>>>>>> Saving best model <<<<<<<<<")
  
 # sp_df2 <- as.data.frame(coordinates(species_df)) %>% 
  #  dplyr::rename(x = 1, y = 2)
  
  sp_df2 <- as.data.frame(coordinates(species_df_raw)) %>% 
    dplyr::rename(x = 1, y = 2)
  
  save_SDM_results(ENMeval_output = eval1, AUCmin = 0.7, resultDir = wd$out,
                   spp = bw, occ_df = sp_df2)
  
  ### 4C. Visualize the model
  ## Load in the rasters
  r    <- raster(file.path(resultDir, paste0(bw,"_SDM.tif")))
  r_pa <- raster(file.path(resultDir, paste0(bw,"_SDM_PA.tif")))
  
  create_sdmFigure(
    spp = bw, r = r, r_pa = r_pa, occ_df = sp_df2, world = world, resultDir = resultDir, bw = bw
  )
  
  print(paste0("!!!!!! Model for ", gsub("_", " ", bw), " complete !!!!!"))
  
  # Clear tempdir
  options(java.parameters = "-Xmx8000m")
 # unixtools::set.tempdir("tmp/")
 # tmp_dir <- unixtools::set.tempdir("tmp/")
  #files <- list.files(tmp_dir, full.names = T,  all.files = T, recursive = T)
  #file.remove(files)
  gc()
  
}

sdm_pipeline(17)


-# looping through pipeline
for(i in seq_along(spp_list)){
  tryCatch(sdm_pipeline(x = i),
           error = function(e) print(paste(spp_list[i], "Error in Code Skipping for now!")))
}  

