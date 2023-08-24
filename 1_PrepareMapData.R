#### Setup ####

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
suppressMessages(suppressWarnings(library(rgdal)))
suppressMessages(suppressWarnings(library(raster)))

# Print session information to log file
sessionInfo()

#### Natural Habitat Map ####

# Define Behrmann equal-area projection
behrCRS <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')

cat('% natural habitat\n')

# Read map of primary habitat
pri <- readGDAL(paste(
  dataDir,"pri_1km_int",sep=""),
  silent = TRUE)

# Create natural-habitat map, initially containing coverage of primary habitats
nat <- pri

# Clean-up
rm(pri)
gc()

# Read map of secondary habitat
sec <- readGDAL(paste(
  dataDir,"sec_1km_int",sep=""),
  silent = TRUE)

# Sum primary and secondary habitat as natural habitat
nat$band1 <- nat$band1 + sec$band1

rm(sec)
gc()

# Convert to raster
nat <- raster(nat)

# Project natural habitat map to Behrmann equal-area projection
nat <- raster::projectRaster(from = nat,res = 1000,crs = behrCRS)

# Round back to integer values
values(nat) <- round(values(nat))

# Create raster to depict 5-km blocks
toRas <- nat

# Define longitudinal steps to create 5-km blocks
step <- ncol(nat)/5

# Define number of blocks in the latitudinal direction
nBlocks <- ceiling(nrow(nat)/5)

# Create raster values depicting the identity of 5-km blocks
cat('Splitting map\n')

for (block in 1:nBlocks){
  
  blockStart <- ((block-1) * step) + 1

  inds <- rep(blockStart:(blockStart+step - 1),each=5)
  
  if (block != nBlocks){
    inds <- rep(inds,5)
  } else {
    inds <- rep(inds,4)
  }
  
  if(block==1){
    blockVals <- inds
  } else {
    blockVals <- c(blockVals,inds)
  }
}

# Assign block identity values to raster map
values(toRas) <- blockVals

# Summarize natural habitat values in 1-km map by 5-km blocks
cat('Calculating values\n')

zm <- zonal(nat,toRas)

# Create a new 5-km-resolution raster to hold summarized natural-habitat values
e <- extent(nat)
e[3] <- e[3]-1000

newRas <- raster(nrows=nBlocks,ncols=step,ext=e,crs=behrCRS)

# Assign summarized natural-habitat values to new 5-km raster
values(newRas) <- zm[,2]

# Write final raster to file
writeRaster(x = newRas,filename = paste(outDir,"PercentNatural.tif",sep=""),format="GTiff")

# Clean-up
rm(newRas,blockVals,nat,zm)
gc()

#### Fertilizer Map ####

# List directories containing estimates of fertilizer application for 17 crops
cropDirs <- dir(path = fertilizerDataDir,pattern = "Fertilizer_",full.names = TRUE)

# Loop over directories for each crop, summing resulting fertilizer maps
total.application <- sum(stack(lapply(cropDirs,function(cd){
  
  cat(paste0(gsub("0_data/FertilizerCropSpecific_Geotiff/Fertilizer_","",cd),"\n"))
  
  # List files that contain N, P and K fertilizer application rates for the crop
  dataFiles <- dir(path = cd,pattern = "_Total.tif",full.names = TRUE)
  dataFiles <- dataFiles[!(grepl("aux",dataFiles))]
  dataFiles <- dataFiles[!(grepl("ovr",dataFiles))]
  
  # Loop over the files for N, P and K fertilizer
  fert.application <- sum(stack(lapply(dataFiles,function(f){
    
    appl <- raster(f)
    
    return(appl)
    
  })),na.rm=TRUE)
  
  return(fert.application)
  
})),na.rm=TRUE)

writeRaster(x = total.application,filename = paste0(outDir,"FertilizerMap.tif"),format="GTiff")

#### Pesticide Map ####

# List all files in the directory of pesticide application maps for the year 2015
pest.files <- list.files(pesticideDataDir, pattern = "2015")

# Separate out the files for the low and high estimates of application densities
pest.files_H <- pest.files[grep(pest.files, pattern = "2015_H")]
pest.files_L <- pest.files[grep(pest.files, pattern = "2015_L")]

# Stack the rasters for the high estimates
pest_H <- stack(paste0(pesticideDataDir, "/", pest.files_H))

# Reclassify the various no-data values as NA:
# -2 = Water; -1.5 = No data; -1 = Banned or not approved substances
pest_H <- reclassify(pest_H, matrix(
  data = c(-2, -1.5, -1, NA, NA, NA), byrow = F, ncol = 2), right = F)

# Sum high estimates of pesticide application across all substances and all crops
pest_H_total <- calc(x = pest_H, fun = sum, na.rm = TRUE)

# Output the final map for high estimates of total pesticide application rate
writeRaster(pest_H_total, filename = paste0(outdir, "/Pesticide_totalAPR_High.tif"))

# Stack the rasters for the low estimates
pest_L <- stack(paste0(pesticideDataDir, "/", pest.files_L))

# Reclassify the various no-data values as NA:
# -2 = Water; -1.5 = No data; -1 = Banned or not approved substances
pest_L <- reclassify(pest_L, matrix(
  data = c(-2, -1.5, -1, NA, NA, NA), byrow = F, ncol = 2), right = F)

# Sum low estimates of pesticide application across all substances and all crops
pest_L_total <- calc(x = pest_L, fun = sum, na.rm = TRUE)

# Output the final map for low estimates of total pesticide application rate
writeRaster(pest_L_total, filename = paste0(outdir, "/Pesticide_totalAPR_Low.tif"))

#### Finish ####

t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()
