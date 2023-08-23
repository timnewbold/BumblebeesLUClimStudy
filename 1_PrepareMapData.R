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

if(!dir.exists(outDir)) dir.create(outDir)

sink(paste0(outDir,"log.txt"))

t.start <- Sys.time()

print(t.start)

suppressMessages(suppressWarnings(library(rgdal)))
suppressMessages(suppressWarnings(library(terra)))

sessionInfo()

behrCRS <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')

cat('% natural habitat\n')

pri <- readGDAL(paste(
  dataDir,"pri_1km_int",sep=""),
  silent = TRUE)

nat <- pri

rm(pri)
gc()

sec <- readGDAL(paste(
  dataDir,"sec_1km_int",sep=""),
  silent = TRUE)

nat$band1 <- nat$band1 + sec$band1

rm(sec)
gc()

nat <- rast(nat)

nat <- raster(nat)

nat <- projectRaster(from = nat,res = 1000,crs = behrCRS)

values(nat) <- round(values(nat))

toRas <- nat

step <- ncol(nat)/5

nBlocks <- ceiling(nrow(nat)/5)

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

values(toRas) <- blockVals

cat('Calculating values\n')

zm <- zonal(nat,toRas)

e <- extent(nat)
e@ymin <- e@ymin-1000

newRas <- raster(nrows=nBlocks,ncols=step,ext=e,crs=behrCRS)
values(newRas) <- zm[,2]

writeRaster(x = newRas,filename = paste(outDir,"PercentNatural.tif",sep=""),format="GTiff")

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

cropDirs <- dir(path = fertilizerDataDir,pattern = "Fertilizer_",full.names = TRUE)

total.application <- sum(stack(lapply(cropDirs,function(cd){
  
  cat(paste0(gsub("FertilizerCropSpecific_Geotiff/Fertilizer_","",cd),"\n"))
  
  dataFiles <- dir(path = cd,pattern = "_Total.tif",full.names = TRUE)
  dataFiles <- dataFiles[!(grepl("aux",dataFiles))]
  dataFiles <- dataFiles[!(grepl("ovr",dataFiles))]
  
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
