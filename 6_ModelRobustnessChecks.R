#### Setup ####

# Set input and output directories
dataDir <- "0_data"
inDir <- "2_PrepareDiversityData/"

outDir <- "6_ModelRobustnessChecks/"

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
diversity$LogPestToxHigh <- log(diversity$PestToxHigh+1)

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
                          'LogPestToxLow','LogPestToxHigh',
                          'NaturalHabitat1k','NaturalHabitat2k','NaturalHabitat5k',
                          'AgeConv','AgeConv10','AgeConv50')]

# Remove rows with any NA values
modelData <- na.omit(modelData)

# Run final model as in main text

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

saveRDS(model,paste0(outDir,"ModelOriginal.rds"))

# Run a version of the model with natural habitat estimated at different grains

model2 <- brm(formula = occur~
               # Control for elevational effects
               poly(LogElevation,1)+
               # Main effects of land use, and landscape pesticide toxicity and
               # conversion age
               LandUse+
               poly(NaturalHabitat1k,1)+
               poly(LogPestToxLow,1)+
               poly(AgeConv,1)+
               # Interactions between land use and other variables
               LandUse:poly(NaturalHabitat1k,1)+
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

saveRDS(model2,paste0(outDir,"ModelNatHab1k.rds"))

model3 <- brm(formula = occur~
                # Control for elevational effects
                poly(LogElevation,1)+
                # Main effects of land use, and landscape pesticide toxicity and
                # conversion age
                LandUse+
                poly(NaturalHabitat5k,1)+
                poly(LogPestToxLow,1)+
                poly(AgeConv,1)+
                # Interactions between land use and other variables
                LandUse:poly(NaturalHabitat5k,1)+
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

saveRDS(model3,paste0(outDir,"ModelNatHab5k.rds"))

# Run models with 10% and 50% thresholds of conversion for calculating duration
# since substantial landscape modification

model4 <- brm(formula = occur~
                # Control for elevational effects
                poly(LogElevation,1)+
                # Main effects of land use, and landscape pesticide toxicity and
                # conversion age
                LandUse+
                poly(NaturalHabitat2k,1)+
                poly(LogPestToxLow,1)+
                poly(AgeConv10,1)+
                # Interactions between land use and other variables
                LandUse:poly(NaturalHabitat2k,1)+
                LandUse:poly(LogPestToxLow,1)+
                LandUse:poly(AgeConv10,1)+
                # Climate niche variables and interactions
                poly(TEI_BL,2)+poly(TEI_delta,1)+
                LandUse:poly(TEI_BL,2)+
                LandUse:poly(TEI_delta,1)+
                # Random effects
                (1|SS)+(1|SSBS)+(1|Taxon_name_entered),
              data=modelData,family='bernoulli',
              iter=2000,chains=4,cores=4)

saveRDS(model4,paste0(outDir,"ModelAge10.rds"))

model5 <- brm(formula = occur~
                # Control for elevational effects
                poly(LogElevation,1)+
                # Main effects of land use, and landscape pesticide toxicity and
                # conversion age
                LandUse+
                poly(NaturalHabitat2k,1)+
                poly(LogPestToxLow,1)+
                poly(AgeConv50,1)+
                # Interactions between land use and other variables
                LandUse:poly(NaturalHabitat2k,1)+
                LandUse:poly(LogPestToxLow,1)+
                LandUse:poly(AgeConv50,1)+
                # Climate niche variables and interactions
                poly(TEI_BL,2)+poly(TEI_delta,1)+
                LandUse:poly(TEI_BL,2)+
                LandUse:poly(TEI_delta,1)+
                # Random effects
                (1|SS)+(1|SSBS)+(1|Taxon_name_entered),
              data=modelData,family='bernoulli',
              iter=2000,chains=4,cores=4)

saveRDS(model5,paste0(outDir,"ModelAge50.rds"))

# Finally, run a model with high instead of low estimates of pesticide toxicity

model6 <- brm(formula = occur~
                # Control for elevational effects
                poly(LogElevation,1)+
                # Main effects of land use, and landscape pesticide toxicity and
                # conversion age
                LandUse+
                poly(NaturalHabitat2k,1)+
                poly(LogPestToxHigh,1)+
                poly(AgeConv,1)+
                # Interactions between land use and other variables
                LandUse:poly(NaturalHabitat2k,1)+
                LandUse:poly(LogPestToxHigh,1)+
                LandUse:poly(AgeConv,1)+
                # Climate niche variables and interactions
                poly(TEI_BL,2)+poly(TEI_delta,1)+
                LandUse:poly(TEI_BL,2)+
                LandUse:poly(TEI_delta,1)+
                # Random effects
                (1|SS)+(1|SSBS)+(1|Taxon_name_entered),
              data=modelData,family='bernoulli',
              iter=2000,chains=4,cores=4)

saveRDS(model6,paste0(outDir,"ModelPestHigh.rds"))