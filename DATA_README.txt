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
## DOI: 10.6084/m9.figshare.24081441
## The creation of these estimates is described in Soroye et al. (2020), DOI: 
## 10.1126/science.aax8591
## The files need are named:
## bombus_stacknames.txt
## BaselineTEI_Spp.rds
## BaselinePEI_Spp.rds
## DeltaTEI_Spp_Period3.rds
## DeltaPEI_Spp_Period3.rds

## LAND-USE HISTORY
## Data on the year in which landscapes were first 30% converted to human land
## uses is also available in the dedicated Figshare repository for this paper:
## DOI: 10.6084/m9.figshare.24081441
## For more information on these estimates, see Newbold et al. (2015), DOI:
## 10.1038/nature14324
## The file needed is named yocPS15002005_1_0.3.asc

## UN SUB-REGIONS MAP
## A map of UN sub-regions, for plottingpurposes, can be downladed from the
## dedicated Figshare repository for this paper:
## DOI: 10.6084/m9.figshare.24081441
## The files needed are the 6 shapefile files named UN_subregion
