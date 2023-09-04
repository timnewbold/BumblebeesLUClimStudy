#### Setup ####

# Set input and output directories
dataDir <- "2_PrepareDiversityData/"
inDir <- "3_RunModels/"

outDir <- "5_GetEstimatedBiodiversityChanges/"

# Create output directory, if it doesn't already exist
if(!dir.exists(outDir)) dir.create(outDir)

# Create log file
sink(paste0(outDir,"log.txt"))

# Start the process timer
t.start <- Sys.time()

print(t.start)

# Load required packages
suppressMessages(suppressWarnings(library(StatisticalModels)))

# Print session information
sessionInfo()

# Load the compiled species-level data
diversity <- readRDS(paste0(dataDir,"diversity_data.rds"))

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
finalModels <- readRDS(paste(inDir,"FinalModels.rds",sep=""))

#### Predictions for sampled locations ####

# Create baseline data-frame with all natural land use, 100% surrounding natural
# habitat (note that 100% is stored as 1,000), no pesticide and no climate change
nd.bl <- modelData
nd.bl$LandUse <- factor("Natural",levels=c("Natural","Human"))
nd.bl$NaturalHabitat <- 1000
nd.bl$Pesticide <- 0
nd.bl$LogFertilizer <- 0
nd.bl$AgeConv <- 0
nd.bl$TEI_delta <- 0

# Extract observed data-frame
nd.obs <- modelData

# Get AIC values of final set of models
aics <- unlist(lapply(finalModels,AIC))

# Get differences in AIC values from best-fitting model
aic.diffs <- aics - min(aics)

# Calculate AIC weights
aic.weights <- exp(-0.5 * aic.diffs)/sum(exp(-0.5 * aic.diffs))

preds <- sapply(X = 1:1000, FUN = function(i) {
  
  ## Select model at random, weighted by AIC weight
  model <- finalModels[[sample(x = 1:length(finalModels),size = 1,
                  replace = FALSE,prob = aic.weights)]]
  
  ## Draw model coefficients
  fe.draw <- mvrnorm(n = 1, mu = fixef(model), 
                     Sigma = vcov(model))
  
  ## Make predictions using the baseline data
  # Create the model matrix
  mm.bl <- model.matrix(terms(model), nd.bl)
  # Remove any redundant columns from the model matrix
  if (ncol(mm.bl) > length(lme4::fixef(model))) {
    mm.bl <- mm.bl[, -which(!(
      names(mm.bl[1, ]) %in% names(
        lme4::fixef(model))))]
  }
  # Calculate predicted 
  y.bl <- mm.bl %*% fe.draw
  y.bl <- 1/(1+exp(-(y.bl)))
  
  mm.obs <- model.matrix(terms(model), nd.obs)
  if (ncol(mm.obs) > length(lme4::fixef(model))) {
    mm.obs <- mm.obs[, -which(!(
      names(mm.obs[1, ]) %in% names(
        lme4::fixef(model))))]
  }
  
  y.obs <- mm.obs %*% fe.draw
  y.obs <- 1/(1+exp(-(y.obs)))
  
  y <- ((y.obs/y.bl)*100)-100
  
  return(y)
})

y.median <- apply(X = preds,MARGIN = 1,FUN = median)
y.lower <- apply(X = preds,MARGIN = 1,FUN = quantile,probs=0.025)
y.upper <- apply(X = preds,MARGIN = 1,FUN = quantile,probs=0.975)

lu.means.median <- tapply(X = y.median,INDEX = nd.obs$LandUse,FUN = median)
lu.means.lower <- tapply(X = y.lower,INDEX = nd.obs$LandUse,FUN = median)
lu.means.upper <- tapply(X = y.upper,INDEX = nd.obs$LandUse,FUN = median)

cat("\n\n")

cat(paste0("Average decline in natural habitat = ",
    round(lu.means.median['Natural'],1),
    "% (95% CI: ",
    round(lu.means.lower['Natural'],1),
    "% - ",
    round(lu.means.upper['Natural'],1),
    "%)"))

cat("\n\n")

cat(paste0("Average decline in human-dominated land use = ",
           round(lu.means.median['Human'],1),
           "% (95% CI: ",
           round(lu.means.lower['Human'],1),
           "% - ",
           round(lu.means.upper['Human'],1),
           "%)"))

t.end <- Sys.time()

cat("\n\n")

print(round(t.end - t.start,0))

sink()