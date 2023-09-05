#### Setup ####

# Set input and output directories
dataDir <- "2_PrepareDiversityData/"
modelsDir <- "3_RunModels/"

outDir <- "6_RefitBestModel/"

# Create output directory, if it doesn't already exist
if(!dir.exists(outDir)) dir.create(outDir)

# Create log file
sink(paste0(outDir,"log.txt"))

# Start timer
t.start <- Sys.time()

print(t.start)

# Load required packages
suppressMessages(suppressWarnings(library(StatisticalModels)))

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

bestModel <- GLMER(modelData = modelData,responseVar = "occur",
                   fitFamily = "binomial",
                   fixedStruct = "LandUse+poly(NaturalHabitat,2)+poly(Pesticide,1)+LandUse:poly(Pesticide,1)+poly(LogFertilizer,2)+poly(AgeConv,2)+poly(TEI_BL,2)+poly(TEI_delta,1)+LandUse:poly(TEI_BL,2)+LandUse:poly(TEI_delta,1)+poly(TEI_BL,2):poly(TEI_delta,1)",
                   randomStruct = "(1|SS)+(1|SSBS)+(1|Taxon_name_entered)",
                   saveVars = c("Longitude","Latitude"))

saveRDS(object = bestModel,file = paste0(outDir,"BestModel.Rds"))

# End timer
t.end <- Sys.time()

# Display run-time
print(round(t.end - t.start,0))

# Close log file
sink()