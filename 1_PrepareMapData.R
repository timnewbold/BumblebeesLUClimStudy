#### Setup ####

dataDir <- "0_data/"
fertilizerDataDir <- "0_data/FertilizerCropSpecific_Geotiff/"
outDir <- "1_PrepareMapData/"

################################################################################
########################## OBTAINING DATA ######################################
################################################################################

## FERTILIZERS
## Global data on density of application (kg/ha) of nitrogen, phosphorous and 
## potassium fertilizer on 17 major crops can be downloaded here:
## http://www.earthstat.org/nutrient-application-major-crops/
## The downloaded data are contained in a folder 'FertilizerCropSpecific_Geotiff'

## COVERAGE OF NATURAL HABITAT
## Data on the coverage of natural (primary and secondary) habitats can be
## downloaded here:
## https://data.csiro.au/collection/csiro:15276v3
## These data then need to be pre-processed by multiplying by 1000 and converting
## to an integer value, with the resulting raster files saved as:
## 'pri_1km_int' and 'sec_1km_int'

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

#### PROBABLY NOT NEEDED - CHECK!

# cat('Habitat diversity\n')
# 
# clc.files <- dir(paste(dataDir,"consensuslandcover",sep=""))
# 
# clc.map <- raster(paste(dataDir,"consensuslandcover/consensus_full_class_1.tif",sep=""))
# clc.map <- projectRaster(from = clc.map,res = 1000,crs = behrCRS)
# values(clc.map) <- round(values(clc.map))
# 
# toRas <- clc.map
# 
# step <- ncol(clc.map)/5
# 
# nBlocks <- ceiling(nrow(clc.map)/5)
# 
# cat('Splitting map\n')
# 
# for (block in 1:nBlocks){
#   
#   print(block)
#   
#   blockStart <- ((block-1) * step) + 1
#   
#   inds <- rep(blockStart:(blockStart+step - 1),each=5)
#   
#   if (block != nBlocks){
#     inds <- rep(inds,5)
#   } else {
#     inds <- rep(inds,1)
#   }
#   
#   if(block==1){
#     blockVals <- inds
#   } else {
#     blockVals <- c(blockVals,inds)
#   }
# }
# 
# values(toRas) <- blockVals
# 
# rm(blockVals,clc.map,inds)
# gc()
# 
# clc.vals <- array(data = NA,dim = c(12,25,18671963))
# 
# cat('Processing land-cover data\n')
# 
# for(i in 1:length(clc.files)){
#   
#   cat('Processing map',i,'\n')
#   
#   map <- raster(paste(dataDir,"consensuslandcover/",clc.files[i],sep=""))
#   
#   map <- projectRaster(from = map,res = 1000,crs = behrCRS)
#   
#   values(map) <- round(values(map))
#   
#   vals <- do.call('cbind',split(values(map),values(toRas)))
#   
#   clc.vals[i,,] <- vals
#   
#   rm(vals)
#   
#   gc()
# }
# 
# e <- extent(toRas)
# e@ymin <- e@ymin-4000
# 
# newRas <- raster(nrows=nBlocks,ncols=step,ext=e,crs=behrCRS)
# 
# cat('Calculating values\n')
# 
# values(newRas) <- apply(X = clc.vals,MARGIN = 3,FUN = function(x){
#   
#   x <- x/100
#   
#   y <- apply(X = x,MARGIN = 1,FUN = mean,na.rm = TRUE)
#   
#   si <- -sum(y[y>0] * log(y[y>0]))
#   
#   return(si)
# })
# 
# writeRaster(x = newRas,filename = paste(outDir,"HabitatDiversity.tif",sep=""),format="GTiff")

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



t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()
