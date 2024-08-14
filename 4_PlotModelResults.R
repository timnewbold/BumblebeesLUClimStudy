#### Setup ####

# Set input and output directories
dataDir <- "2_PrepareDiversityData/"
inDir <- "3_RunModels/"

outDir <- "4_PlotModelResults/"

# Create output directory, if it doesn't already exist
if(!dir.exists(outDir)) dir.create(outDir)

# Create log file
sink(paste0(outDir,"log.txt"))

# Start timer
t.start <- Sys.time()

print(t.start)

# Load required packages
suppressMessages(suppressWarnings(library(brms)))
suppressMessages(suppressWarnings(library(Hmisc)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(ggstance)))
suppressMessages(suppressWarnings(library(cowplot)))

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

draws <- as_draws_df(finalModel)[,2:15]

draw.means <- apply(X = draws,MARGIN = 2,FUN = mean)
draw.lower.95 <- apply(X = draws,MARGIN = 2,FUN = quantile,probs = 0.025)
draw.upper.95 <- apply(X = draws,MARGIN = 2,FUN = quantile,probs = 0.975)
draw.lower.67 <- apply(X = draws,MARGIN = 2,FUN = quantile,probs = 0.16667)
draw.upper.67 <- apply(X = draws,MARGIN = 2,FUN = quantile,probs = 0.83333)

draw.summary <- data.frame(Variable = names(draw.means),
                           Mean = draw.means,
                           Lower95 = draw.lower.95,
                           Upper95 = draw.upper.95,
                           Lower67 = draw.lower.67,
                           Upper67 = draw.upper.67)

draw.summary$Variable <- gsub("b_","",draw.summary$Variable)
draw.summary$Variable <- gsub("poly","",draw.summary$Variable)
draw.summary$Variable <- gsub("2k1","2k",draw.summary$Variable)
draw.summary$Variable <- gsub("LogPestToxLow1","PesticideToxicity",draw.summary$Variable)
draw.summary$Variable <- gsub("AgeConv1","ConversionDuration",draw.summary$Variable)
draw.summary$Variable <- gsub("TEI_BL21","ThermalPosition_Baseline",draw.summary$Variable)
draw.summary$Variable <- gsub("TEI_BL22","ThermalPosition_Baseline2",draw.summary$Variable)
draw.summary$Variable <- gsub("TEI_delta1","ThermalPosition_Delta",draw.summary$Variable)
draw.summary$Variable <- gsub("LogElevation1","Elevation",draw.summary$Variable)

draw.summary$Variable <- factor(draw.summary$Variable,
                                levels = c("Elevation",
                                           "LandUseHuman:ThermalPosition_Delta",
                                           "ThermalPosition_Delta",
                                           "LandUseHuman:ThermalPosition_Baseline2",
                                           "LandUseHuman:ThermalPosition_Baseline",
                                           "ThermalPosition_Baseline2",
                                           "ThermalPosition_Baseline",
                                           "LandUseHuman:ConversionDuration",
                                           "ConversionDuration",
                                           "LandUseHuman:PesticideToxicity",
                                           "PesticideToxicity",
                                           "LandUseHuman:NaturalHabitat2k",
                                           "NaturalHabitat2k",
                                           "LandUseHuman"))

plotCols <- c("#0072B2","#D55E00","#D55E00",
              "#D55E00","#D55E00","#D55E00","#D55E00",
              "#F0E442","#F0E442","#E69F00","#E69F00",
              "#009E73","#009E73","#CC79A7")

forestPlot <- ggplot() + geom_linerangeh(mapping = aes(y=draw.summary$Variable,
                                        xmin = draw.summary$Lower95, 
                                        xmax = draw.summary$Upper95,
                                        colour = draw.summary$Variable)) + 
  geom_pointrangeh(mapping = aes(y = draw.summary$Variable,
                                 x = draw.summary$Mean,
                                 xmin = draw.summary$Lower67,
                                 xmax = draw.summary$Upper67,
                                 colour = draw.summary$Variable),
                   fatten=1.8,size=1) +
  scale_colour_manual(values = plotCols) +
  geom_vline(mapping = aes(xintercept = 0)) +
  labs(x = "Effect size",y = element_blank()) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
  
save_plot(filename = paste0(outDir,"CoefficientForestPlot.png"),plot = forestPlot,
          base_height = 12/2.54,base_width = 17.5/2.54)

#### Plotting (natural habitat) ####

# Create data frame with variable values for which to predict occurrence probability
# Set elevation and landscape age variables at their median
# values (to reflect a 'typical' landscape)
# Ignore pesticide toxicity for now (i.e., set to zero)
# Set baseline TEI to the median value
# Set delta TEI to zero (i.e. ignore climate change effects)
nd <- data.frame(LogPestToxLow=0,
                 AgeConv=median(modelData$AgeConv),
                 LogElevation=median(modelData$LogElevation),
                 TEI_BL=median(modelData$TEI_BL),
                 TEI_delta=0)

# Get predicted values with both 67% (~ 1 SE) and 95% credible intervals
preds.67 <- conditional_effects(x = finalModel,effects = "NaturalHabitat2k:LandUse",
                                resolution = 1000,
                                conditions = nd,plot = FALSE,prob=0.67)
preds.95 <- conditional_effects(x = finalModel,effects = "NaturalHabitat2k:LandUse",
                                resolution = 1000,
                                conditions = nd,plot = FALSE,prob=0.95)

plotHab <- ggplot() + 
  geom_ribbon(data = preds.95$`NaturalHabitat2k:LandUse`, 
              mapping = aes(x = NaturalHabitat2k/10,ymin = lower__,
                            ymax = upper__,fill = LandUse),alpha = 0.2) + 
  geom_ribbon(data = preds.67$`NaturalHabitat2k:LandUse`,
              mapping = aes(x = NaturalHabitat2k/10,ymin = lower__,
                            ymax = upper__,fill = LandUse),alpha = 0.4) + 
  geom_line(data = preds.67$`NaturalHabitat2k:LandUse`,
            mapping = aes(x = NaturalHabitat2k/10,y = estimate__,col = LandUse),
            linewidth = 1) + 
  scale_fill_manual(values = c("#009F81","#A40122")) + 
  scale_color_manual(values = c("#009F81","#A40122")) + 
  scale_x_continuous(name = "Natural habitat (%)",limits = c(0,100)) + 
  scale_y_continuous(name = "Probability of occurrence",limits = c(0,1)) + 
  theme_classic() + 
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

#### Plotting (landscape conversion age) ####

# Create data frame with variable values for which to predict occurrence probability
# Set elevation and natural habitat variables at their median
# values (to reflect a 'typical' landscape)
# Ignore pesticide toxicity for now (i.e., set to zero)
# Set baseline TEI to the median value
# Set delta TEI to zero (i.e. ignore climate change effects)
nd <- data.frame(LogPestToxLow=0,
                 NaturalHabitat2k=median(modelData$NaturalHabitat2k),
                 LogElevation=median(modelData$LogElevation),
                 TEI_BL=median(modelData$TEI_BL),
                 TEI_delta=0)

# Get predicted values with both 67% (~ 1 SE) and 95% credible intervals
preds.67 <- conditional_effects(x = finalModel,effects = "AgeConv:LandUse",
                                resolution = 1000,
                                conditions = nd,plot = FALSE,prob=0.67)
preds.95 <- conditional_effects(x = finalModel,effects = "AgeConv:LandUse",
                                resolution = 1000,
                                conditions = nd,plot = FALSE,prob=0.95)

plotAge <- ggplot() + 
  geom_ribbon(data = preds.95$`AgeConv:LandUse`, 
              mapping = aes(x = AgeConv,ymin = lower__,
                            ymax = upper__,fill = LandUse),alpha = 0.2) + 
  geom_ribbon(data = preds.67$`AgeConv:LandUse`,
              mapping = aes(x = AgeConv,ymin = lower__,
                            ymax = upper__,fill = LandUse),alpha = 0.4) + 
  geom_line(data = preds.67$`AgeConv:LandUse`,
            mapping = aes(x = AgeConv,y = estimate__,col = LandUse),
            linewidth = 1) + 
  scale_fill_manual(values = c("#009F81","#A40122")) + 
  scale_color_manual(values = c("#009F81","#A40122")) + 
  scale_x_continuous(name = "Duration of modification (yrs)") + 
  scale_y_continuous(name = "Probability of occurrence",limits = c(0,1)) + 
  theme_classic() + 
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))


#### Plotting (pesticide toxicity) ####

# Create data frame with variable values for which to predict occurrence probability
# Set elevation, natural habitat and landscape age variables at their median
# values (to reflect a 'typical' landscape)
# Set baseline TEI to the median value
# Set delta TEI to zero (i.e. ignore climate change effects)
nd <- data.frame(NaturalHabitat2k=median(modelData$NaturalHabitat2k),
                 AgeConv=median(modelData$AgeConv),
                 LogElevation=median(modelData$LogElevation),
                 TEI_BL=median(modelData$TEI_BL),
                 TEI_delta=0)

# Get predicted values with both 67% (~ 1 SE) and 95% credible intervals
preds.67 <- conditional_effects(x = finalModel,effects = "LogPestToxLow:LandUse",
                                resolution = 1000,
                                conditions = nd,plot = FALSE,prob=0.67)
preds.95 <- conditional_effects(x = finalModel,effects = "LogPestToxLow:LandUse",
                                resolution = 1000,
                                conditions = nd,plot = FALSE,prob=0.95)

plotPest <- ggplot() + 
  geom_ribbon(data = preds.95$`LogPestToxLow:LandUse`, 
              mapping = aes(x = exp(LogPestToxLow)-1,ymin = lower__,
                            ymax = upper__,fill = LandUse),alpha = 0.2) + 
  geom_ribbon(data = preds.67$`LogPestToxLow:LandUse`,
              mapping = aes(x = exp(LogPestToxLow)-1,ymin = lower__,
                            ymax = upper__,fill = LandUse),alpha = 0.4) + 
  geom_line(data = preds.67$`LogPestToxLow:LandUse`,
            mapping = aes(x = exp(LogPestToxLow)-1,y = estimate__,col = LandUse),
            linewidth = 1) + 
  scale_fill_manual(values = c("#009F81","#A40122")) + 
  scale_color_manual(values = c("#009F81","#A40122")) + 
  scale_x_continuous(name = "Pesticide toxicity") + 
  scale_y_continuous(name = "Probability of occurrence",limits = c(0,1)) + 
  theme_classic() + 
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))


#### Plotting (Baseline TEI) ####

# Create data frame with variable values for which to predict occurrence probability
# Set elevation, natural habitat, and landscape age variables at their median
# values (to reflect a 'typical' landscape)
# Ignore pesticide toxicity for now (i.e., set to zero)
# Set delta TEI to zero (i.e. no climate change)
nd <- data.frame(NaturalHabitat2k=median(modelData$NaturalHabitat2k),
                 LogPestToxLow=0,
                 AgeConv=median(modelData$AgeConv),
                 LogElevation=median(modelData$LogElevation),
                 TEI_delta=0)

# Get predicted values with both 67% (~ 1 SE) and 95% credible intervals
preds.67 <- conditional_effects(x = finalModel,effects = "TEI_BL:LandUse",
                                resolution = 1000,
                                conditions = nd,plot = FALSE,prob=0.67)
preds.95 <- conditional_effects(x = finalModel,effects = "TEI_BL:LandUse",
                                resolution = 1000,
                                conditions = nd,plot = FALSE,prob=0.95)

plotTEI_bl <- ggplot() + 
  geom_ribbon(data = preds.95$`TEI_BL:LandUse`, 
              mapping = aes(x = TEI_BL,ymin = lower__,
                            ymax = upper__,fill = LandUse),alpha = 0.2) + 
  geom_ribbon(data = preds.67$`TEI_BL:LandUse`,
              mapping = aes(x = TEI_BL,ymin = lower__,
                            ymax = upper__,fill = LandUse),alpha = 0.4) + 
  geom_line(data = preds.67$`TEI_BL:LandUse`,
            mapping = aes(x = TEI_BL,y = estimate__,col = LandUse),
            linewidth = 1) + 
  scale_fill_manual(values = c("#009F81","#A40122")) + 
  scale_color_manual(values = c("#009F81","#A40122")) + 
  scale_x_continuous(name = "Baseline thermal position") + 
  scale_y_continuous(name = "Probability of occurrence",limits = c(0,1)) + 
  theme_classic() + 
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))


#### Plotting (delta TEI) ####

# Create data frame with variable values for which to predict occurrence probability
# Set elevation, natural habitat, and landscape age variables at their median
# values (to reflect a 'typical' landscape)
# Ignore pesticide toxicity for now (i.e., set to zero)
# Set baseline TEI to the median value
nd <- data.frame(NaturalHabitat2k=median(modelData$NaturalHabitat2k),
                 LogPestToxLow=0,
                 AgeConv=median(modelData$AgeConv),
                 LogElevation=median(modelData$LogElevation),
                 TEI_BL=median(modelData$TEI_BL))

# Get predicted values with both 67% (~ 1 SE) and 95% credible intervals
preds.67 <- conditional_effects(x = finalModel,effects = "TEI_delta:LandUse",
                                resolution = 1000,
                                conditions = nd,plot = FALSE,prob=0.67)
preds.95 <- conditional_effects(x = finalModel,effects = "TEI_delta:LandUse",
                                resolution = 1000,
                                conditions = nd,plot = FALSE,prob=0.95)

plotTEI_delta <- ggplot() + 
  geom_ribbon(data = preds.95$`TEI_delta:LandUse`, 
              mapping = aes(x = TEI_delta,ymin = lower__,
                            ymax = upper__,fill = LandUse),alpha = 0.2) + 
  geom_ribbon(data = preds.67$`TEI_delta:LandUse`,
              mapping = aes(x = TEI_delta,ymin = lower__,
                            ymax = upper__,fill = LandUse),alpha = 0.4) + 
  geom_line(data = preds.67$`TEI_delta:LandUse`,
            mapping = aes(x = TEI_delta,y = estimate__,col = LandUse),
            linewidth = 1) + 
  scale_fill_manual(values = c("#009F81","#A40122")) + 
  scale_color_manual(values = c("#009F81","#A40122")) + 
  scale_x_continuous(name = "\u0394 thermal position") + 
  scale_y_continuous(name = "Probability of occurrence",limits = c(0,1)) + 
  theme_classic() + 
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

plots <- plot_grid(plotHab,plotAge,plotPest,plotTEI_bl,plotTEI_delta,
                   nrow = 3,ncol = 2,labels =c("a)","b)","c)","d)","e)"),
                   label_x = -0.02,label_size = 12,label_fontface = "plain")

save_plot(filename = paste0(outDir,"ConditionalEffects.png"),
          plot = plots,ncol = 2,nrow = 3,
          base_height = 7.3333/2.54,base_width = 8.75/2.54)

t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()