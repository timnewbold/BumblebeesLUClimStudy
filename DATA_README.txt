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