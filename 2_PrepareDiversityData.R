#### Setup ####

# Set input and output directories
dataDir <- "0_data/"
mapDir <- "1_PrepareMapData/"

outDir <- "2_PrepareDiversityData/"

################################################################################
########################## OBTAINING DATA ######################################
################################################################################

## All the data below need to be saved in the folder called '0_data'

## PREDICTS DATABASE
## The PREDICTS database, documenting the structure of ecological assemblages
## across different land-use types and intensities, can be downloaded from the 
## data portal of the Natural History Museum
## http://dx.doi.org/10.5519/0066354
## You should download the RDS-format file, and save as database.rds

## ELEVATION
## A Digital Elevation Model at 30 arc-second spatial resolution can be 
## downloaded from WorldClim:
## https://www.worldclim.org/data/worldclim21.html
## Filename should be wc2.1_30s_elev.tif

## BUMBLEBEE NICHE PROPERTIES
## Estimates of the realized climatic niches of bumblebees, and how these have
## changed under recent climate change can be downloaded from the dedicated
## Figshare repository for this paper:
## INSERT LINK
## The creation of these estimates is described in Soroye et al. (2020), DOI: 
## 10.1126/science.aax8591

## LAND-USE HISTORY
## Data on the year in which landscapes were first 30% converted to human land
## uses is also available in the dedicated Figshare repository for this paper:
## INSERT LINK
## For more information on these estimates, see Newbold et al. (2015), DOI:
## 10.1038/nature14324

## UN SUB-REGIONS MAP
## A map of UN sub-regions, for plottingpurposes, can be downladed from the
## dedicated Figshare repository for this paper:
## INSERT LINK

# Create output directory, if it doesn't already exist
if(!dir.exists(outDir)) dir.create(outDir)

# Create log file
sink(paste0(outDir,"log.txt"))

# Start timer
t.start <- Sys.time()

print(t.start)

# Load required packages
suppressMessages(suppressWarnings(library(terra)))
suppressMessages(suppressWarnings(library(predictsFunctions)))
suppressMessages(suppressWarnings(library(raster)))
suppressMessages(suppressWarnings(library(rgdal)))
# To obtain this package, run:
# remotes::install_github("timnewbold/predicts-demo",subdir="predictsFunctions")

# Print session information to log file
sessionInfo()

cat('Loading database extracts\n')
# Read the PREDICTS database
diversity<-readRDS(paste(dataDir,"database.rds",sep=""))

cat('Selecting appropriate data\n')
# Remove sites without geographical coordinates
diversity <- diversity[apply(diversity[,c('Longitude','Latitude')],1,
                             function(r) all(!is.na(r))),]

# Select only data for the genus Bombus
diversity <- droplevels(diversity[(diversity$Genus=="Bombus"),])

# Select data from northern, southern and western Europe, and North America
diversity <- droplevels(diversity[(diversity$UN_subregion %in% 
                                     c("Northern Europe","Southern Europe",
                                       "Western Europe","North America")),])

cat('Correcting for sampling effort\n')
# Correct for any variation in sampling effort within studies
diversity <- CorrectSamplingEffort(diversity)

cat('Merging sites\n')
# Merge sites that have identical coordinates and environmental variables
diversity <- MergeSites(diversity,silent = TRUE)

# Drop data where taxonomy could not be resolved
diversity <- diversity[(diversity$Best_guess_binomial!=""),]
diversity <- diversity[(diversity$Best_guess_binomial!="Unknown bombus"),]

diversity <- droplevels(diversity)

# Read information on climatic position variables for the bumblebees
sp.names <- read.table(paste(dataDir,"bombus_stacknames.txt",sep=""))

# Define WGS84 coordinate system
wgsCRS <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'

TEI_bl <- readRDS(paste(dataDir,"BaselineTEI_Spp.rds",sep=""))
names(TEI_bl) <- paste("Bombus_",sp.names$x,sep="")
TEI_bl <- rast(TEI_bl)
crs(TEI_bl) <- wgsCRS
PEI_bl <- readRDS(paste(dataDir,"BaselinePEI_Spp.rds",sep=""))
names(PEI_bl) <- paste("Bombus_",sp.names$x,sep="")
PEI_bl <- rast(PEI_bl)
crs(PEI_bl) <- wgsCRS
TEI_delta <- readRDS(paste(dataDir,"DeltaTEI_Spp_Period3.rds",sep=""))
names(TEI_delta) <- paste("Bombus_",sp.names$x,sep="")
TEI_delta <- rast(TEI_delta)
crs(TEI_delta) <- wgsCRS
PEI_delta <- readRDS(paste(dataDir,"DeltaPEI_Spp_Period3.rds",sep=""))
names(PEI_delta) <- paste("Bombus_",sp.names$x,sep="")
PEI_delta <- rast(PEI_delta)
crs(PEI_delta) <- wgsCRS

## Check correspondence between measures based on old and new CRU data

# Get species names for layers based on CRU4
sp.names.4 <- read.csv(paste0(dataDir,"ThermalLimits_CRU4/species_names.csv"))

TEI_bl_4 <- readRDS(paste0(dataDir,"ThermalLimits_CRU4/baselineTPI_Spp.RDS"))
names(TEI_bl_4) <- sp.names.4$species
TEI_bl_4 <- rast(TEI_bl_4)
TEI_bl_4 <- terra::project(x = TEI_bl_4,y = TEI_bl)

TEI_delta_4 <- readRDS(paste0(dataDir,"ThermalLimits_CRU4/deltaTPI_Spp.RDS"))
names(TEI_delta_4) <- sp.names.4$species
TEI_delta_4 <- rast(TEI_delta_4)
TEI_delta_4 <- terra::project(x = TEI_delta_4,y = TEI_delta)

cor.tei.bl <- sapply(X = intersect(sp.names$x,sp.names.4$species),FUN = function(sp){
  
  rast.old <- TEI_bl[[paste0("Bombus_",sp)]]
  rast.4 <- TEI_bl_4[[sp]]
  
  return(cor.test(values(rast.old),values(rast.4))$estimate)
  
})

cor.tei.delta <- sapply(X = intersect(sp.names$x,sp.names.4$species),FUN = function(sp){
  
  rast.old <- TEI_delta[[paste0("Bombus_",sp)]]
  rast.4 <- TEI_delta_4[[sp]]
  
  return(cor.test(values(rast.old),values(rast.4))$estimate)
  
})

# Add climatic position variables to the data frame
diversity <- do.call('rbind',lapply(X = split(
  x = diversity,f = diversity$Best_guess_binomial),
  FUN = function(div.sp){
  
  sp <- gsub(" ","_",div.sp$Best_guess_binomial[1])
  sp2 <- strsplit(paste(div.sp$Best_guess_binomial[1])," ")[[1]][2]
  
  if (sp %in% names(TEI_bl)){
    div.sp$TEI_BL <- terra::extract(
      x = TEI_bl[[sp]],y = div.sp[,c('Longitude','Latitude')])[[sp]]
  } else {
    div.sp$TEI_BL <- NA
  }
  
  if (sp %in% names(PEI_bl)){
    div.sp$PEI_BL <- terra::extract(
      x = PEI_bl[[sp]],y = div.sp[,c('Longitude','Latitude')])[[sp]]
  } else {
    div.sp$PEI_BL <- NA
  }
  
  if (sp %in% names(TEI_delta)){
    div.sp$TEI_delta <- terra::extract(
      x = TEI_delta[[sp]],y = div.sp[,c('Longitude','Latitude')])[[sp]]
  } else {
    div.sp$TEI_delta <- NA
  }
  
  if (sp %in% names(PEI_delta)){
    div.sp$PEI_delta <- terra::extract(
      x = PEI_delta[[sp]],y = div.sp[,c('Longitude','Latitude')])[[sp]]
  } else {
    div.sp$PEI_delta <- NA
  }
  
  if (sp2 %in% names(TEI_bl_4)){
    div.sp$TEI_BL_4 <- terra::extract(
      x = TEI_bl_4[[sp2]],y = div.sp[,c('Longitude','Latitude')])[[sp2]]
  } else {
    div.sp$TEI_BL_4 <- NA
  }
  
  if (sp2 %in% names(TEI_delta_4)){
    div.sp$TEI_delta_4 <- terra::extract(
      x = TEI_delta_4[[sp2]],y = div.sp[,c('Longitude','Latitude')])[[sp2]]
  } else {
    div.sp$TEI_delta_4 <- NA
  }
  
  return(div.sp)
}))

# Reclassify land uses as natural or human-dominated
diversity$LandUse <- paste(diversity$Predominant_land_use)
diversity$LandUse[(diversity$LandUse=="Primary vegetation")] <- "Natural"
diversity$LandUse[(diversity$LandUse=="Mature secondary vegetation")] <- "Natural"
diversity$LandUse[(diversity$LandUse=="Intermediate secondary vegetation")] <- "Natural"
diversity$LandUse[(diversity$LandUse=="Young secondary vegetation")] <- "Natural"
diversity$LandUse[(diversity$LandUse=="Secondary vegetation (indeterminate age)")] <- "Natural"
diversity$LandUse[(diversity$LandUse=="Cropland")] <- "Human"
diversity$LandUse[(diversity$LandUse=="Pasture")] <- "Human"
diversity$LandUse[(diversity$LandUse=="Urban")] <- "Human"
diversity$LandUse[(diversity$LandUse=="Cannot decide")] <- NA
diversity$LandUse <- factor(diversity$LandUse,levels=c("Natural","Human"))

# Create a binary occurrence column (1 where recorded abundance was > 0)
diversity$occur <- ifelse(diversity$Measurement>0,1,0)

# Define equal-area (Behrmann) and geographic (WGS84) projections
behrCRS <- '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs'
wgsCRS <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'

# Read elevation and add to data frame
elev <- rast(paste0(dataDir,"wc2.1_30s_elev.tif"))

diversity$Elevation <- terra::extract(
  x = elev,y = diversity[,c('Longitude','Latitude')])$wc2.1_30s_elev

# Read data on landscape percent natural habitat and add to data frame
nathab.1 <- rast(paste0(mapDir,"PercentNaturalHabitat1k.tif"))
nathab.2 <- rast(paste0(mapDir,"PercentNaturalHabitat2k.tif"))
nathab.5 <- rast(paste0(mapDir,"PercentNaturalHabitat5k.tif"))
nathab.10 <- rast(paste0(mapDir,"PercentNaturalHabitat10k.tif"))

# Project to geographic coordiinates, for matching with the PREDICTS sites
nathab.1 <- terra::project(x = nathab.1,y = wgsCRS)
nathab.2 <- terra::project(x = nathab.2,y = wgsCRS)
nathab.5 <- terra::project(x = nathab.5,y = wgsCRS)
nathab.10 <- terra::project(x = nathab.10,y = wgsCRS)

# Extract values to the biodiversity data frame
diversity$NaturalHabitat1k <- terra::extract(
  x = nathab.1,y = diversity[,c('Longitude','Latitude')])$pri_1km_int
diversity$NaturalHabitat2k <- terra::extract(
  x = nathab.2,y = diversity[,c('Longitude','Latitude')])$pri_1km_int
diversity$NaturalHabitat5k <- terra::extract(
  x = nathab.5,y = diversity[,c('Longitude','Latitude')])$pri_1km_int
diversity$NaturalHabitat10k <- terra::extract(
  x = nathab.10,y = diversity[,c('Longitude','Latitude')])$pri_1km_int

# Read estimates of pesticide toxicity and add to data frame
pest.tox.l <- rast(paste0(mapDir,"Pesticide_TotalToxicity_Low.tif"))
pest.tox.h <- rast(paste0(mapDir,"Pesticide_TotalToxicity_High.tif"))

# Replace NA values with zeroes
values(pest.tox.l)[(is.na(values(pest.tox.l)))] <- 0
values(pest.tox.h)[(is.na(values(pest.tox.h)))] <- 0

diversity$PestToxLow <- terra::extract(
  x = pest.tox.l,y = diversity[,c('Longitude','Latitude')])$sum

diversity$PestToxHigh <- terra::extract(
  x = pest.tox.h,y = diversity[,c('Longitude','Latitude')])$sum

# Read data on fertilizer application and add to data frame
fert <- rast(paste0(mapDir,"FertilizerMap.tif"))

# Replace NA values with zeroes
values(fert)[is.na(values(fert))] <- 0

diversity$Fertilizer <- terra::extract(
  x = fert,y = diversity[,c('Longitude','Latitude')])$sum

# Read data on time since 30% landscape conversion to human uses 
YOC30 <- rast(paste0(dataDir,"yocPS15002005_1_0.3.asc"))
# Where land hasn't yet been converted, set year of conversion to most recent year
values(YOC30)[values(YOC30)==0] <- 2005
# Calculate number of years from landscape conversion to most recent year
AgeConv <- 2005 - YOC30

# Add time since conversion to data frame
diversity$AgeConv <- terra::extract(
  x = AgeConv,y = diversity[,c('Longitude','Latitude')])$yocPS15002005_1_0.3

cat('Saving diversity data\n')
# Save data for modelling
saveRDS(object = diversity,file=paste0(outDir,"diversity_data.rds"))

# Create data frame of unique sites
sites <- unique(diversity[,c('SSS','Longitude','Latitude')])

# Make map of sites
occMap <- SpatialPointsDataFrame(coords = data.frame(x=sites$Longitude,
                                                     y=sites$Latitude),
                                 data = sites,
                                 proj4string = CRS(
                                   "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Read map of UN sub-regions, and select relevant regions (Europe and N America)
un_sub <- readOGR(dsn = gsub("/","",dataDir),layer = "UN_subregion",verbose = FALSE)

un_sub <- un_sub[(un_sub$SUBREGION %in% c(21,13,154,155,39)),]

# Define plotting colours for UN sub-regions
cols <- colorRampPalette(colors = c("#000000","#ffffff"))(5)

# Plot map of sites
png(filename = paste(outDir,"MapRecords.png",sep=""),width = 17.5,height = 8,units = "cm",res = 1200)

par(mar=c(0,0,0,0))

plot(un_sub,col=cols,xlim=c(-180,30))

plot(occMap,add=TRUE,pch=16,cex=0.5,col="#cc0000")

invisible(dev.off())

# End process timer
t.end <- Sys.time()

print(round(t.end - t.start,0))

# sink()