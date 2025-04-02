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
suppressMessages(suppressWarnings(library(brms))) # Version 2.21.0
suppressMessages(suppressWarnings(library(parallel))) # Version 4.3.1
suppressMessages(suppressWarnings(library(snow))) # Version 0.4-4
suppressMessages(suppressWarnings(library(terra))) # Version 1.7-71
suppressMessages(suppressWarnings(library(ggplot2))) # Version 3.5.1
suppressMessages(suppressWarnings(library(ggmap))) # Version 4.0.0
suppressMessages(suppressWarnings(library(cowplot))) # Version 1.1.3
suppressMessages(suppressWarnings(library(sf))) # Version 1.0-19

# Print session information
sessionInfo()

# Load the compiled species-level data
diversity <- readRDS(paste0(dataDir,"diversity_data.rds"))

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

# Load the final model
finalModel <- readRDS(paste0(inDir,"FinalModelBRMS.rds"))

#### Predictions for sampled locations ####

# Create baseline data-frame with all natural land use, 100% surrounding natural
# habitat (note that 100% is stored as 1,000), no pesticide and no climate change
nd.bl <- modelData
nd.bl$LandUse <- factor("Natural",levels=c("Natural","Human"))
nd.bl$NaturalHabitat2k <- 1000
nd.bl$LogPestToxLow <- 0
nd.bl$AgeConv <- 0
nd.bl$TEI_delta <- 0

# Extract observed data-frame
nd.obs <- modelData

# cl <- makeCluster(detectCores()-1)
cl <- makeCluster(2)
clusterExport(cl = cl,list = c("finalModel","nd.bl","nd.obs"))

# Make 1,000 predictions of % difference between observed
# and baseline conditions
preds <- parSapply(cl = cl,X = 1:10,FUN = function(i){
  
  # Draw one of the posterior parameter estimates at random
  i <- sample(x = 1:4000,size = 1,replace = FALSE)
  
  # Use this set of posterior parameter estimates to get predicted
  # probability of occurrence under baseline and observed environmental
  # conditions
  y.bl <- brms::posterior_epred(
    object = finalModel,newdata = nd.bl,re_formula = NA,
    draw_ids = i)[1,]
  
  y.obs <- brms::posterior_epred(
    object = finalModel,newdata = nd.obs,re_formula = NA,
    draw_ids = i)[1,]
  
  y <- ((y.obs/y.bl)*100)-100
  
  return(y)
})

stopCluster(cl)

# Calculate median and credible intervals across replicate predictions
y.median <- apply(X = preds,MARGIN = 1,FUN = median)
y.lower.95 <- apply(X = preds,MARGIN = 1,FUN = quantile,probs=0.025)
y.upper.95 <- apply(X = preds,MARGIN = 1,FUN = quantile,probs=0.975)
y.lower.67 <- apply(X = preds,MARGIN = 1,FUN = quantile,probs=0.16667)
y.upper.67 <- apply(X = preds,MARGIN = 1,FUN = quantile,probs=0.83333)

# Summarize median, and credible intervals of predictions
# across the two land-use types
lu <- factor(c("Natural","Human"),levels=c("Natural","Human"))
lu.means.median <- tapply(X = y.median,INDEX = nd.obs$LandUse,FUN = median)
lu.means.lower.95 <- tapply(X = y.lower.95,INDEX = nd.obs$LandUse,FUN = median)
lu.means.upper.95 <- tapply(X = y.upper.95,INDEX = nd.obs$LandUse,FUN = median)
lu.means.lower.67 <- tapply(X = y.lower.67,INDEX = nd.obs$LandUse,FUN = median)
lu.means.upper.67 <- tapply(X = y.upper.67,INDEX = nd.obs$LandUse,FUN = median)

# Print land-use summaries
cat("\n\n")

cat(paste0("Average decline in natural habitat = ",
    round(lu.means.median['Natural'],1),
    "% (95% CI: ",
    round(lu.means.lower.95['Natural'],1),
    "% - ",
    round(lu.means.upper.95['Natural'],1),
    "%)"))

cat("\n\n")

cat(paste0("Average decline in human-dominated land use = ",
           round(lu.means.median['Human'],1),
           "% (95% CI: ",
           round(lu.means.lower.95['Human'],1),
           "% - ",
           round(lu.means.upper.95['Human'],1),
           "%)"))

cat("\n\n")

# Average predicted differences in occurrence probability across sampled sites
y.median.sites <- tapply(X = y.median,INDEX = nd.obs$SSBS,FUN = mean)
y.upper.sites <- tapply(X = y.upper.95,INDEX = nd.obs$SSBS,FUN = mean)
y.lower.sites <- tapply(X = y.lower.95,INDEX = nd.obs$SSBS,FUN = mean)

# Create data frame for site-level averages
sites <- data.frame(y=y.median.sites,yplus=y.upper.sites,yminus=y.lower.sites,
                    SSBS=names(y.median.sites))
sites$Longitude <- nd.obs$Longitude[match(sites$SSBS,nd.obs$SSBS)]
sites$Latitude <- nd.obs$Latitude[match(sites$SSBS,nd.obs$SSBS)]

# Create basemap for North American and Western European sub-regions
baseMap <- vect("0_data/UN_subregion.shp")
baseMap <- terra::subset(baseMap,baseMap$SUBREGION %in% c(21,13,154,155,39))
baseMap <- st_as_sf(baseMap)

# Create map of predictions for sites sampled in PREDICTS
predsMap <- vect(x = sites,geom=c("Longitude","Latitude"),
                 crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
predsMap <- st_as_sf(predsMap)

# Plot base map and predictions for sampled locations
mapPreds <- ggplot() + geom_sf(data = baseMap,fill = "#ffffff") + 
  geom_sf(data=predsMap,mapping = aes(colour=y)) + 
  scale_color_continuous("\u0394 Probability of\noccurrence (%)",type = "viridis") + 
  scale_x_continuous(limits = c(-180,60)) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size=10))

# Create error-bar plot showing predictions at sampled locations
# by land-use type
plotLU <- ggplot() + 
  geom_linerange(mapping = aes(x = lu,
                               ymin = lu.means.lower.95,
                               ymax = lu.means.upper.95,
                               col = lu)) +
  geom_pointrange(mapping = aes(x = lu,
                                y = lu.means.median,
                                ymin = lu.means.lower.67,
                                ymax = lu.means.upper.67,
                                col = lu),
                  linewidth = 1, size = 0.8) +
  scale_colour_manual(values = c("#009F81","#A40122")) +
  scale_y_continuous(name = "\u0394 Probability of occurrence (%)",
                     limits = c(-100,0)) + 
  scale_x_discrete(name = "Land use") + 
  geom_hline(mapping = aes(yintercept = 0),
             alpha = 0.3, linetype = 'dashed') +
  theme_classic() + 
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

# Create final figure
plots <- cowplot::plot_grid(mapPreds,plotLU,nrow = 2,ncol = 1,
                   rel_heights = c(0.4,0.6),rel_widths = c(1,0.1))
save_plot(filename = paste0(outDir,"EstimatedBiodiversityDifferences.png"),
          plot = plots,ncol = 1,nrow = 2,base_width = 12.5/2.54,
          base_height = 6.5/2.54)

t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()