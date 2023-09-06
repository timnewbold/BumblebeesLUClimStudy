#### Setup ####

# Set input and output directories
dataDir <- "2_PrepareDiversityData/"
modelsDir <- "3_RunModels/"
bestModelDir <- "6_RefitBestModel/"

outDir <- "7_FitMCMCModel/"

# Create output directory, if it doesn't already exist
if(!dir.exists(outDir)) dir.create(outDir)

# Create log file
sink(paste0(outDir,"log.txt"))

# Start timer
t.start <- Sys.time()

print(t.start)

# Load required packages
suppressMessages(suppressWarnings(library(StatisticalModels)))
suppressMessages(suppressWarnings(library(MCMCglmm)))

# Print session information
sessionInfo()

# Load the compiled species-level data
diversity <- readRDS(paste0(dataDir,"diversity_data.rds"))

allData <- readRDS(paste0(modelsDir,"ModelData.Rds"))

sites <- readRDS(paste0(modelsDir,"AnalysisSites.rds"))

# Log-transform elevation, which is very skewed
diversity$LogElevation <- log(diversity$Elevation+2)

# Log-transform fertilizer application, which is also very skewed
diversity$LogFertilizer <- log(diversity$Fertilizer+1)

diversity$UI <- paste(diversity$LandUse,diversity$Use_intensity,sep="_")
diversity$UI[grepl("Natural",diversity$UI)] <- "Natural"
diversity$UI[grepl("NA",diversity$UI)] <- NA
diversity$UI[grepl("Cannot decide",diversity$UI)] <- NA
diversity$UI <- factor(diversity$UI)
diversity$UI <- relevel(diversity$UI,ref="Natural")

# Select relevant columns
modelData <- diversity[,c('occur','LandUse','TEI_BL',
                          'TEI_delta','LogElevation','SS','SSBS',
                          'Taxon_name_entered','Longitude','Latitude',
                          'Best_guess_binomial','Country',
                          'Pesticide','LogFertilizer',
                          'NaturalHabitat','AgeConv')]
# Remove rows with any NA values
modelData <- na.omit(modelData)

# Load the final model
finalModels <- readRDS(paste(modelsDir,"FinalModels.rds",sep=""))

# Load the best model refitted using the GLMER function
bestModel <- readRDS(paste0(bestModelDir,"BestModel.Rds"))

# Fit an MCMC model with the same structure as the best-fitting model
mcmcModel <- MCMCglmm(fixed = occur~LandUse+poly(NaturalHabitat,2)+poly(Pesticide,1)+LandUse:poly(Pesticide,1)+poly(LogFertilizer,2)+poly(AgeConv,2)+poly(TEI_BL,2)+poly(TEI_delta,1)+LandUse:poly(TEI_BL,2)+LandUse:poly(TEI_delta,1)+poly(TEI_BL,2):poly(TEI_delta,1),random = ~SS+SSBS+Taxon_name_entered,data = modelData,family = "categorical",pr=TRUE)

saveRDS(object = mcmcModel,file = paste0(outDir,"MCMCModel.Rds"))

mcmcPreds <- predict.MCMCglmm(object = mcmcModel,newdata = NULL,verbose=TRUE,marginal = NULL)

predsCompare <- data.frame(predsLM4=fitted(bestModel$model),predsMCMC=mcmcPreds)

saveRDS(object = predsCompare,file = paste0(outDir,"MCMCLME4PredictionsCompared.Rds"))

# End timer
t.end <- Sys.time()

# Display run-time
print(round(t.end - t.start,0))

# Close log file
sink()