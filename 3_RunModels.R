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
suppressMessages(suppressWarnings(library(brms)))

# Print session information
sessionInfo()

# Load the compiled species-level data
diversity <- readRDS(paste0(inDir,"diversity_data.rds"))

# Log-transform elevation, which is very skewed
diversity$LogElevation <- log(diversity$Elevation+2)

# Log-transform fertilizer application, which is also very skewed
diversity$LogFertilizer <- log(diversity$Fertilizer+1)

# Log-transform pesticide application, which is also very skewed
diversity$LogPestToxLow <- log(diversity$PestToxLow+1)

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
                          'LogPestToxLow',
                          'NaturalHabitat2k','AgeConv')]
# Remove rows with any NA values
modelData <- na.omit(modelData)

saveRDS(object = modelData,file = paste0(outDir,"ModelData.rds"))

# Create data frame of unique sites
sites <- unique(modelData[,c('SSBS','Longitude','Latitude')])

saveRDS(object = sites,file = paste0(outDir,"AnalysisSites.rds"))

# Run model using BRMS
model <- brm(formula = occur~
               # Control for elevational effects
               poly(LogElevation,1)+
               # Main effects of land use, and landscape pesticide toxicity and
               # conversion age
               LandUse+
               poly(NaturalHabitat2k,1)+
               poly(LogPestToxLow,1)+
               poly(AgeConv,1)+
               # Interactions between land use and other variables
               LandUse:poly(NaturalHabitat2k,1)+
               LandUse:poly(LogPestToxLow,1)+
               LandUse:poly(AgeConv,1)+
               # Climate niche variables and interactions
               poly(TEI_BL,2)+poly(TEI_delta,1)+
               LandUse:poly(TEI_BL,2)+
               LandUse:poly(TEI_delta,1)+
               # Random effects
               (1|SS)+(1|SSBS)+(1|Taxon_name_entered),
             data=modelData,family='bernoulli',
             iter=2000,chains=4,cores=4)

saveRDS(model,paste0(outDir,"FinalModelBRMS.rds"))

pdf(file = paste0(outDir,"ModelResultsDiagnostics.pdf"),width = 16,height = 24)

plot(model,ask = FALSE)

invisible(dev.off())

# End timer
t.end <- Sys.time()

# Display run-time
print(round(t.end - t.start,0))

# Close log file
sink()