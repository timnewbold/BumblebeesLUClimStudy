#### Setup ####

# Set input and output directories
dataDir <- "0_data"
inDir <- "2_PrepareDiversityData/"

outDir <- "3_RunModels/"

# Create output directory, if it doesn't already exist
if(!dir.exists(outDir)) dir.create(outDir)

# Create log file
sink(paste0(outDir,"log.txt"))

# Start timer
t.start <- Sys.time()

print(t.start)

# Load required packages
suppressMessages(suppressWarnings(library(StatisticalModels)))
# To install this package, run the following code:
# remotes::install_github("timnewbold/StatisticalModels")
suppressMessages(suppressWarnings(library(rgdal)))
suppressMessages(suppressWarnings(library(DHARMa)))

# Print session information
sessionInfo()

# Load the compiled species-level data
diversity <- readRDS(paste0(inDir,"diversity_data.rds"))

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

saveRDS(object = modelData,file = paste0(outDir,"ModelData.rds"))

# Create data frame of unique sites
sites <- unique(modelData[,c('SSBS','Longitude','Latitude')])

saveRDS(object = sites,file = paste0(outDir,"AnalysisSites.rds"))

# Candidate fixed-effect structures
candidate.fixefs <- list(
  ## Null model
  "1",
  ## Land-use only
  "LandUse",
  ## Single variable groups (in addition to land use):
  # Landscape natural habitat
  "LandUse+poly(NaturalHabitat,2)",
  # Intensity factors - fertilizer and pesticide
  "LandUse+poly(Pesticide,1)+LandUse:poly(Pesticide,1)+poly(LogFertilizer,2)",
  # Land-use history
  "LandUse+poly(AgeConv,2)",
  # Climatic variables
  "LandUse+poly(TEI_BL,2)+poly(TEI_delta,1)+LandUse:poly(TEI_BL,2)+LandUse:poly(TEI_delta,1)+poly(TEI_BL,2):poly(TEI_delta,1)",
  ## Two groups of variables at a time
  # Natural habitat and intensity
  "LandUse+poly(NaturalHabitat,2)+poly(Pesticide,1)+LandUse:poly(Pesticide,1)+poly(LogFertilizer,2)",
  # Natural habitat and land-use history
  "LandUse+poly(NaturalHabitat,2)+poly(AgeConv,2)",
  # Natural habitat and climate
  "LandUse+poly(NaturalHabitat,2)+poly(TEI_BL,2)+poly(TEI_delta,1)+LandUse:poly(TEI_BL,2)+LandUse:poly(TEI_delta,1)+poly(TEI_BL,2):poly(TEI_delta,1)",
  # Intensity and land-use history
  "LandUse+poly(Pesticide,1)+LandUse:poly(Pesticide,1)+poly(LogFertilizer,2)+poly(AgeConv,2)",
  # Intensity and climate
  "LandUse+poly(Pesticide,1)+LandUse:poly(Pesticide,1)+poly(LogFertilizer,2)+poly(TEI_BL,2)+poly(TEI_delta,1)+LandUse:poly(TEI_BL,2)+LandUse:poly(TEI_delta,1)+poly(TEI_BL,2):poly(TEI_delta,1)",
  # Land-use history and climate
  "LandUse+poly(AgeConv,2)+poly(TEI_BL,2)+poly(TEI_delta,1)+LandUse:poly(TEI_BL,2)+LandUse:poly(TEI_delta,1)+poly(TEI_BL,2):poly(TEI_delta,1)",
  ## Three groups of variables at a time
  # Natural habitat, intensity and land-use history
  "LandUse+poly(NaturalHabitat,2)+poly(Pesticide,1)+LandUse:poly(Pesticide,1)+poly(LogFertilizer,2)+poly(AgeConv,2)",
  # Natural habitat, intensity and climate
  "LandUse+poly(NaturalHabitat,2)+poly(Pesticide,1)+LandUse:poly(Pesticide,1)+poly(LogFertilizer,2)+poly(TEI_BL,2)+poly(TEI_delta,1)+LandUse:poly(TEI_BL,2)+LandUse:poly(TEI_delta,1)+poly(TEI_BL,2):poly(TEI_delta,1)",
  # Natural habitat, land-use history and climate
  "LandUse+poly(NaturalHabitat,2)+poly(AgeConv,2)+poly(TEI_BL,2)+poly(TEI_delta,1)+LandUse:poly(TEI_BL,2)+LandUse:poly(TEI_delta,1)+poly(TEI_BL,2):poly(TEI_delta,1)",
  # Intensity, land-use history and climate
  "LandUse+poly(Pesticide,1)+LandUse:poly(Pesticide,1)+poly(LogFertilizer,2)+poly(AgeConv,2)+poly(TEI_BL,2)+poly(TEI_delta,1)+LandUse:poly(TEI_BL,2)+LandUse:poly(TEI_delta,1)+poly(TEI_BL,2):poly(TEI_delta,1)",
  ## All four variable groups together
  "LandUse+poly(NaturalHabitat,2)+poly(Pesticide,1)+LandUse:poly(Pesticide,1)+poly(LogFertilizer,2)+poly(AgeConv,2)+poly(TEI_BL,2)+poly(TEI_delta,1)+LandUse:poly(TEI_BL,2)+LandUse:poly(TEI_delta,1)+poly(TEI_BL,2):poly(TEI_delta,1)"
)

# Run a set of models with these candidate fixed-effects structures
models <- GLMERCandidates(modelData = modelData,responseVar = "occur",
                          fitFamily = "binomial",
                          candidateFixedStructs = candidate.fixefs,
                          randomStruct = "(1|SS)+(1|SSBS)+(1|Taxon_name_entered)")

# Save the final models
saveRDS(object = models,file = paste0(outDir,"FinalModels.rds"))

# Get the AIC values of each model
aics <- unlist(lapply(models,AIC))

# Identify the best-fitting model (lowest AIC)
bestModel <- which(aics==min(aics))

# Test for over-dispersion in the residuals of the best model
DHARMa::testDispersion(models[[bestModel]],alternative = "greater")

# Run model performance check
simOut <- DHARMa::simulateResiduals(fittedModel = models[[bestModel]])

# Plot results of model performance check
plot(simOut)

# End timer
t.end <- Sys.time()

# Display run-time
print(round(t.end - t.start,0))

# Close log file
sink()