#### Setup ####

# Set input and output directories
dataDir <- "0_data/"
fertilizerDataDir <- "0_data/FertilizerCropSpecific_Geotiff/"
pesticideDataDir <- "0_data/PEST-CHEMGRIDS_v1_01_APR/GEOTIFF/"
outDir <- "1_PrepareMapData/"

################################################################################
########################## OBTAINING DATA ######################################
################################################################################

## All the data below need to be saved in a folder called '0_data'

## COVERAGE OF NATURAL HABITAT
## Data on the coverage of natural (primary and secondary) habitats can be
## downloaded here:
## https://data.csiro.au/collection/csiro:15276v3
## These data then need to be pre-processed by multiplying by 1000 and converting
## to an integer value, with the resulting raster files saved as:
## 'pri_1km_int' and 'sec_1km_int'

## FERTILIZERS
## Global data on the density of application (kg/ha) of nitrogen, phosphorous and 
## potassium fertilizer on 17 major crops can be downloaded here:
## http://www.earthstat.org/nutrient-application-major-crops/
## The downloaded data need to be in a folder 'FertilizerCropSpecific_Geotiff'

## PESTICIDES
## Global data on the density of application (kg/ha) of different types of pesticides
## on different agricultural crops can be downloaded (in GeoTIFF format) here:
## https://sedac.ciesin.columbia.edu/data/set/ferman-v1-pest-chemgrids-v1-01/data-download
## The downloaded data need to be in the following directory structure:
## 'PEST-CHEMGRIDS_v1_01_APR/GEOTIFF'

# Create output directory, if it doesn't already exist
if(!dir.exists(outDir)) dir.create(outDir)

# Create log file
sink(paste0(outDir,"log.txt"))

# Start timer
t.start <- Sys.time()

print(t.start)

# Load required packages
suppressMessages(suppressWarnings(library(terra)))

# Print session information to log file
sessionInfo()

#### Natural Habitat Map ####
cat('Natural habitat\n')

# Define Behrmann equal-area projection
behrCRS <- '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs'

# Read in mapped estimates of primary vegetation
pri <- terra::rast(paste0(dataDir,"pri_1km_int"))

# Read in mapped estimates of secondary vegetation
sec <- terra::rast(paste0(dataDir,"sec_1km_int"))

# Combine primary and secondary vegetation to create estimates of total natural habitat
nat <- pri + sec

# Reproject to Behrmann equal-area projection with 1-km resolution
nat <- terra::project(x = nat,y = behrCRS,res = 1000)

# Aggregate natural habitat estimates into 2- 5 and 10-km blocks
nat.2 <- terra::aggregate(x = nat,fact = 2,fun = 'mean',na.rm = TRUE)
nat.5 <- terra::aggregate(x = nat,fact = 5,fun = 'mean',na.rm = TRUE)
nat.10 <- terra::aggregate(x = nat,fact = 10,fun = 'mean',na.rm = TRUE)

# Write out natural rasters showing percent natural habitat at all grid resolutions
terra::writeRaster(x = nat,filename = paste0(outDir,"PercentNaturalHabitat1k.tif"))
terra::writeRaster(x = nat.2,filename = paste0(outDir,"PercentNaturalHabitat2k.tif"))
terra::writeRaster(x = nat.5,filename = paste0(outDir,"PercentNaturalHabitat5k.tif"))
terra::writeRaster(x = nat.10,filename = paste0(outDir,"PercentNaturalHabitat10k.tif"))

#### Fertilizer Map ####
cat('Fertiliser\n')

# List directories containing estimates of fertilizer application for 17 crops
cropDirs <- dir(path = fertilizerDataDir,pattern = "Fertilizer_",full.names = TRUE)

# Loop over directories for each crop, summing resulting fertilizer maps
application.rate <- sum(rast(lapply(cropDirs,function(cd){
  
  cat(paste0(gsub("0_data/FertilizerCropSpecific_Geotiff/Fertilizer_","",cd),"\n"))
  
  # List files that contain N, P and K fertilizer application rates for the crop
  dataFiles <- dir(path = cd,pattern = "_Rate.tif",full.names = TRUE)
  dataFiles <- dataFiles[!(grepl("aux",dataFiles))]
  dataFiles <- dataFiles[!(grepl("ovr",dataFiles))]
  
  # Loop over the files for N, P and K fertilizer
  fert.application <- sum(rast(lapply(dataFiles,function(f){
    
    appl <- rast(f)
    
    return(appl)
    
  })),na.rm=TRUE)
  
  return(fert.application)
  
})),na.rm=TRUE)

writeRaster(x = application.rate,filename = paste0(outDir,"FertilizerMap.tif"))

#### Pesticide Map ####
cat("Pesticides\n")

# Read table of LD50 values for honeybees
haz_tab <- read.csv(paste0(dataDir, "/TABLE_LD50_PesticidePropertiesDatabase.csv"))

# List all files in the directory of pesticide application maps for the year 2015
pest.files <- list.files(pesticideDataDir, pattern = "2015")

# Separate out the files for the low and high estimates of application densities
pest.files_H <- pest.files[grep(pest.files, pattern = "2015_H")]
pest.files_L <- pest.files[grep(pest.files, pattern = "2015_L")]

# Stack the rasters for the high estimates
pest_H <- rast(paste0(pesticideDataDir, "/", pest.files_H))

# Reclassify the various no-data values as NA:
# -2 = Water; -1.5 = No data; -1 = Banned or not approved substances
pest_H <- classify(pest_H, matrix(
  data = c(-2, -1.5, -1, NA, NA, NA), byrow = F, ncol = 2), right = F)

# For each substance applied to each crop, divide application rate by honeybee LD50s
pest_tox_H <- rast(lapply(X = pest_H,FUN = function(pest.map){
  
  pest_name <- strsplit(x = names(pest.map),split = '_')[[1]][3]
  
  # extract the honeybee contact LD50 value from the table
  LD50.contact <- haz_tab[haz_tab$pesticide == pest_name, "Honeybees_contact"]
  LD50.oral <- haz_tab[haz_tab$pesticide == pest_name, "Honeybees_oral"]
  
  
  if (!is.na(LD50.contact)){
    # If LD50 contact estimate is available, use that preferentially
    cat(paste0(pest_name," - CONTACT","\n"))
    
    tox.map <- pest.map/LD50.contact
    
  } else if (!is.na(LD50.oral)){
    # Otherwise, if LD50 oral estimate is available, use that
    cat(paste0(pest_name," - ORAL","\n"))
    
    tox.map <- pest.map/LD50.oral
  } else {
    # Otherwise, set this pesticide to zero
    cat(paste0(pest_name," - MISSING","\n"))
    
    tox.map <- pest.map*0
  }
  
  return(tox.map)
}))

# Sum estimates of pesticide toxicity across all substances and all crops
pest_tox_H_total <- sum(pest_tox_H,na.rm = TRUE)

# Output the final map for high estimates of total pesticide application rate
writeRaster(pest_tox_H_total, filename = paste0(outDir, "/Pesticide_TotalToxicity_High.tif"))

# Stack the rasters for the low estimates
pest_L <- rast(paste0(pesticideDataDir, "/", pest.files_L))

# Reclassify the various no-data values as NA:
# -2 = Water; -1.5 = No data; -1 = Banned or not approved substances
pest_L <- classify(pest_L, matrix(
  data = c(-2, -1.5, -1, NA, NA, NA), byrow = F, ncol = 2), right = F)

# For each substance applied to each crop, divide application rate by honeybee LD50s
pest_tox_L <- rast(lapply(X = pest_L,FUN = function(pest.map){
  
  pest_name <- strsplit(x = names(pest.map),split = '_')[[1]][3]
  
  # extract the honeybee contact LD50 value from the table
  LD50.contact <- haz_tab[haz_tab$pesticide == pest_name, "Honeybees_contact"]
  LD50.oral <- haz_tab[haz_tab$pesticide == pest_name, "Honeybees_oral"]
  
  
  if (!is.na(LD50.contact)){
    # If LD50 contact estimate is available, use that preferentially
    cat(paste0(pest_name," - CONTACT","\n"))
    
    tox.map <- pest.map/LD50.contact
    
  } else if (!is.na(LD50.oral)){
    # Otherwise, if LD50 oral estimate is available, use that
    cat(paste0(pest_name," - ORAL","\n"))
    
    tox.map <- pest.map/LD50.oral
  } else {
    # Otherwise, set this pesticide to zero
    cat(paste0(pest_name," - MISSING","\n"))
    
    tox.map <- pest.map*0
  }
  
  return(tox.map)
}))

# Sum estimates of pesticide toxicity across all substances and all crops
pest_tox_L_total <- sum(pest_tox_L,na.rm = TRUE)

# Output the final map for low estimates of total pesticide application rate
writeRaster(pest_tox_L_total, filename = paste0(outDir, "/Pesticide_TotalToxicity_Low.tif"))

#### Finish ####

t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()