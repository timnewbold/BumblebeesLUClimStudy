#### Setup ####

# Set input and output directories
dataDir <- "0_data/"
modelDir <- "3_RunModels/"

outDir <- "7_DataChecks/"

# Create output directory, if it doesn't already exist
if(!dir.exists(outDir)) dir.create(outDir)

# Create log file
sink(paste0(outDir,"log.txt"))

# Start timer
t.start <- Sys.time()

print(t.start)

# Load required packages
suppressMessages(suppressWarnings(library(terra))) # Version 1.7-71

# Print session information
sessionInfo()

# Load data used for the main models
modelData <- readRDS(paste0(modelDir,"ModelData.rds"))

## Check numbers of stations used to estimate temperature in the CRU dataset
## for sites in our final model data

# Read estimates of the number of weather stations whose data influenced a given
# grid cell's temperature estimate
# January 1901
stn.1901 <- rast(paste0(dataDir,"cru_ts3.24.01.1901.1910.tmp.stn.nc"))[[1]]
# December 2014
stn.2014 <- rast(paste0(dataDir,"cru_ts3.24.01.2011.2015.tmp.stn.nc"))[[48]]

modelData$stn.1901 <- terra::extract(x = stn.1901,y = modelData[,c('Longitude','Latitude')])$stn_1
modelData$stn.2014 <- terra::extract(x = stn.2014,y = modelData[,c('Longitude','Latitude')])$stn_48

cat(paste0("Minimum number of weather stations supporting temperature estimates in January 1901 is ",min(modelData$stn.1901),"\n"))
cat(paste0("Minimum number of weather stations supporting temperature estimates in December 2014 is ",min(modelData$stn.2014),"\n"))

## Compare maximum temperature estimates at the sampled sites between CRU Version
## 3.24.01 (used in the main analysis) and 4.08 (the most recent version)
tmx.1901.3.24 <- rast(paste0(dataDir,"cru_ts3.24.01.1901.1910.tmx.dat.nc"))[[1]]
tmx.1901.4.08 <- rast(paste0(dataDir,"cru_ts4.08.1901.1910.tmx.dat.nc"))[[1]]

tmx.2014.3.24 <- rast(paste0(dataDir,"cru_ts3.24.01.2011.2015.tmx.dat.nc"))[[48]]
tmx.2014.4.08 <- rast(paste0(dataDir,"cru_ts4.08.2011.2020.tmx.dat.nc"))[[48]]

modelData$tmx.1901.3.24 <- terra::extract(x = tmx.1901.3.24,y = modelData[,c('Longitude','Latitude')])$tmx_1
modelData$tmx.1901.4.08 <- terra::extract(x = tmx.1901.4.08,y = modelData[,c('Longitude','Latitude')])$tmx_1

modelData$tmx.2014.3.24 <- terra::extract(x = tmx.2014.3.24,y = modelData[,c('Longitude','Latitude')])$tmx_48
modelData$tmx.2014.4.08 <- terra::extract(x = tmx.2014.4.08,y = modelData[,c('Longitude','Latitude')])$tmx_48


cat("Correlation of maximum temperatures between CRU versions for January 1901:\n")
print(cor.test(modelData$tmx.1901.3.24,modelData$tmx.1901.4.08))

cat("Correlation of maximum temperatures between CRU versions for December 2014:\n")
print(cor.test(modelData$tmx.2014.3.24,modelData$tmx.2014.4.08))

# End timer
t.end <- Sys.time()

# Display run-time
print(round(t.end - t.start,0))

# Close log file
sink()