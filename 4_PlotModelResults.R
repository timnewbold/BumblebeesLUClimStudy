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

#### Plotting (land use) ####

# Plot results
pdf(file = paste0(outDir,"ResultsPlotLandUse.pdf"),
    width = 8.5/2.54,height = 6/2.54)

par(tck=-0.01)
par(mar=c(2.6,2.6,0.2,0.2))
par(mgp=c(1.6,0.2,0))
par(las=1)

### Land use
# Create data frame with variable values for which to predict occurrence probability
# Start with rows for natural and human land use
nd <- data.frame(LandUse=factor(c("Natural","Human")))
nd$LandUse <- relevel(x = nd$LandUse,ref = "Natural")
# Set natural habitat, pesticide and landscape age variables at their median
# values for each land-use type
nd$NaturalHabitat <- c(
  median(modelData[modelData$LandUse=="Natural",]$NaturalHabitat),
  median(modelData[modelData$LandUse=="Human",]$NaturalHabitat))
nd$Pesticide <- c(
  median(modelData[modelData$LandUse=="Natural",]$Pesticide),
  median(modelData[modelData$LandUse=="Human",]$Pesticide))
nd$LogFertilizer <- c(
  median(modelData[modelData$LandUse=="Natural",]$LogFertilizer),
  median(modelData[modelData$LandUse=="Human",]$LogFertilizer))
nd$AgeConv <- c(
  median(modelData[modelData$LandUse=="Natural",]$AgeConv),
  median(modelData[modelData$LandUse=="Human",]$AgeConv))

# Set baseline TEI to be near range centre, but within 
# sampled range (10th percentile of sampled values, i.e. = 0.59)
nd$TEI_BL <- 0.59
# Set delta TEI to zero (i.e. no climate change)
nd$TEI_delta <- 0
# Add dummy variable representing the response variable
nd$occur <- 0

# Make predictions from model
preds <- PredictGLMERMultiModel(model = finalModels,data = nd,nIters = 1000)

# Back-transform to occurrence probabilities
preds <- 1/(1+exp(-preds))

# Rescale predictions to reference level (i.e. all probabilities = 1 in 
# natural habitat)
preds <- sweep(x = preds,MARGIN = 2,STATS = preds[1,],FUN = '/')

# Calculate summary statistics across predictions
preds.median <- ((apply(X = preds,MARGIN = 1,FUN = median,na.rm=TRUE))*100)-100
preds.lower.67 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=(1/6),
                          na.rm=TRUE))*100)-100
preds.upper.67 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=(5/6),
                          na.rm=TRUE))*100)-100
preds.lower.95 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=0.025,
                       na.rm=TRUE))*100)-100
preds.upper.95 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=0.975,
                       na.rm=TRUE))*100)-100

# Plot error bars for 67% and 95% confidence intervals
errbar(x = 1:2,y = preds.median,yplus = preds.upper.95,yminus = preds.lower.95,
       pch=16,cap = 0,xaxt="n",xlim=c(0.6,2.4),lwd=1.5,cex = 1.6,bty="l",
       errbar.col = c("#009F81","#A40122"),col = c("#009F81","#A40122"),
       xlab=NA,ylab="P.occur difference (%)")
errbar(x = 1:2,y = preds.median,yplus = preds.upper.67,yminus = preds.lower.67,
       pch=16,cap = 0,xaxt="n",xlim=c(0.6,2.4),lwd=3.5,cex = 1.6,add = TRUE,
       errbar.col = c("#009F81","#A40122"),col = c("#009F81","#A40122"),
       xlab=NA,ylab="P.occur difference (%)")

axis(side = 1,at = c(1,2),labels = c("Natural","Human-used"))
title(xlab="Land use")
abline(h=0,lty=2,col="#00000055")

invisible(dev.off())

#### Plotting (landscape) ####

# Plot results
pdf(file = paste0(outDir,"ResultsPlotLandscape.pdf"),
    width = 17.5/2.54,height = 12/2.54)

par(mfrow=c(2,2))
par(tck=-0.01)
par(mar=c(2.6,2.6,0.2,0.2))
par(mgp=c(1.6,0.2,0))
par(las=1)


### Natural habitat

# Get the 95% ranges of values for natural habitat for each land-use type
natHabRangeNatural <- quantile(
  modelData$NaturalHabitat[modelData$LandUse=="Natural"],
  probs = c(0.025,0.975))
natHabRangeHuman <- quantile(
  modelData$NaturalHabitat[modelData$LandUse=="Human"],
  probs = c(0.025,0.975))

# Create data-frame with first row as reference (100% natural habitat), and 
# then plotting ranges for natural and human land use
# Note 100% natural habitat = a value of 1,000
nd <- data.frame(NaturalHabitat=c(1000,
  seq(from = natHabRangeNatural[1],to = natHabRangeNatural[2],
      length.out = 1000),
  seq(from = natHabRangeHuman[1],to = natHabRangeHuman[2],
      length.out = 1000)))
# Reference row takes natural habitat, then 1000 values each of natural
# and human land use
nd$LandUse <- factor(c(rep("Natural",1001),rep("Human",1000)))
nd$LandUse <- relevel(x = nd$LandUse,ref = "Natural")
# For reference row, set pesticide to 0, then medians for each land-use type
nd$Pesticide <- c(0,
  rep(median(modelData[modelData$LandUse=="Natural",]$Pesticide),
      1000),
  rep(median(modelData[modelData$LandUse=="Human",]$Pesticide),
      1000))
# For reference row, set fertilizer to 0, then medians for each land-use type
nd$LogFertilizer <- c(0,
                  rep(median(modelData[modelData$LandUse=="Natural",]$LogFertilizer),
                      1000),
                  rep(median(modelData[modelData$LandUse=="Human",]$LogFertilizer),
                      1000))
# For reference row, set landscape age to median across whole dataset, then
# median values for each land-use type
nd$AgeConv <- c(median(modelData$AgeConv),
  rep(median(modelData[modelData$LandUse=="Natural",]$AgeConv),
      1000),
  rep(median(modelData[modelData$LandUse=="Human",]$AgeConv),1000))
# Set baseline TEI to be near range centre, but within 
# sampled range (10th percentile of sampled values, i.e. = 0.59)
nd$TEI_BL <- 0.59
# Set delta TEI to 0 (i.e. no climate change)
nd$TEI_delta <- 0
# Add dummy variable representing the response variable
nd$occur <- 0

# Make predictions from model
preds <- PredictGLMERMultiModel(models = finalModels,data = nd,nIters = 1000)

# Back-transform to occurrence probabilities
preds <- 1/(1+exp(-preds))

# Rescale predictions to reference level (i.e. all probabilities = 1 in 
# natural habitat)
preds <- sweep(x = preds,MARGIN = 2,STATS = preds[1,],FUN = '/')

# Calculate summary statistics across predictions
preds.median <- ((apply(X = preds,MARGIN = 1,FUN = median,na.rm=TRUE))*100)-100
preds.lower.67 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=(1/6),
                          na.rm=TRUE))*100)-100
preds.upper.67 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=(5/6),
                          na.rm=TRUE))*100)-100
preds.lower.95 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=0.025,
                          na.rm=TRUE))*100)-100
preds.upper.95 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=0.975,
                          na.rm=TRUE))*100)-100

plot(-9e99,-9e99,xlim=c(0,1000),
     ylim=c(min(preds.lower.95),max(preds.upper.95)),xaxt="n",
     xlab=NA,ylab="P.occur difference (%)",bty="l")
axis(side = 1,at = c(0,250,500,750,1000),labels = c(0,25,50,75,100))
title(xlab="Landscape natural habitat (%)")

X.Vec <- c(nd$NaturalHabitat[2:1001], max(nd$NaturalHabitat[2:1001]), 
           rev(nd$NaturalHabitat[2:1001]), min(nd$NaturalHabitat[2:1001]))
Y.Vec <- c(preds.lower.95[2:1001], tail(preds.upper.95[2:1001], 1), 
           rev(preds.upper.95[2:1001]), (preds.lower.95[2:1001])[1])
Y.Vec2 <- c(preds.lower.67[2:1001], tail(preds.upper.67[2:1001], 1), 
           rev(preds.upper.67[2:1001]), (preds.lower.67[2:1001])[1])

polygon(x = X.Vec,y = Y.Vec,col="#009F8133",border=NA)
polygon(x = X.Vec,y = Y.Vec2,col="#009F8155",border=NA)
points(x = nd$NaturalHabitat[2:1001],y = preds.median[2:1001],type="l",lwd=2,col="#009F81")

abline(h=0,lty=2,col="#00000055")

plot(-9e99,-9e99,xlim=c(0,1000),
     ylim=c(min(preds.lower.95),max(preds.upper.95)),xaxt="n",
     xlab=NA,ylab="P.occur difference (%)",bty="l")
axis(side = 1,at = c(0,250,500,750,1000),labels = c(0,25,50,75,100))
title(xlab="Landscape natural habitat (%)")

X.Vec <- c(nd$NaturalHabitat[1002:2001], max(nd$NaturalHabitat[1002:2001]), 
           rev(nd$NaturalHabitat[1002:2001]), min(nd$NaturalHabitat[1002:2001]))
Y.Vec <- c(preds.lower.95[1002:2001], tail(preds.upper.95[1002:2001], 1), 
           rev(preds.upper.95[1002:2001]), (preds.lower.95[1002:2001])[1])
Y.Vec2 <- c(preds.lower.67[1002:2001], tail(preds.upper.67[1002:2001], 1), 
            rev(preds.upper.67[1002:2001]), (preds.lower.67[1002:2001])[1])

polygon(x = X.Vec,y = Y.Vec,col="#A4012233",border=NA)
polygon(x = X.Vec,y = Y.Vec2,col="#A4012255",border=NA)
points(x = nd$NaturalHabitat[1002:2001],y = preds.median[1002:2001],
       type="l",lwd=2,col="#A40122")

abline(h=0,lty=2,col="#00000055")

### Landscape conversion age

# Get the 95% ranges of values for landscape conversion age for each land-use type
ageConvRangeNatural <- quantile(
  modelData$AgeConv[modelData$LandUse=="Natural"],
  probs = c(0.025,0.975))
ageConvRangeHuman <- quantile(
  modelData$AgeConv[modelData$LandUse=="Human"],
  probs = c(0.025,0.975))

# Create data-frame with first row as reference (the median landscape age 
# across the whole dataset), and  then plotting ranges for natural and 
# human land use
nd <- data.frame(AgeConv=c(median(modelData$AgeConv),
                           seq(from = ageConvRangeNatural[1],
                               to = ageConvRangeNatural[2],
                               length.out = 1000),
                           seq(from = ageConvRangeHuman[1],
                               to = ageConvRangeHuman[2],
                               length.out = 1000)))
# Reference row takes natural habitat, then 1000 values each of natural
# and human land use
nd$LandUse <- factor(c(rep("Natural",1001),rep("Human",1000)))
nd$LandUse <- relevel(x = nd$LandUse,ref = "Natural")
# For reference row, set pesticide to 0, then medians for each land-use type
nd$Pesticide <- c(0,
                  rep(median(modelData[modelData$LandUse=="Natural",]$Pesticide),
                      1000),
                  rep(median(modelData[modelData$LandUse=="Human",]$Pesticide),
                      1000))
# For reference row, set fertilizer to 0, then medians for each land-use type
nd$LogFertilizer <- c(0,
                  rep(median(modelData[modelData$LandUse=="Natural",]$LogFertilizer),
                      1000),
                  rep(median(modelData[modelData$LandUse=="Human",]$LogFertilizer),
                      1000))
# For reference row, set natural habitat to 100%, then
# median values for each land-use type
# Note 100% natural habitat = a value of 1,000
nd$NaturalHabitat <- c(1000,
                       rep(median(modelData[
                         modelData$LandUse=="Natural",]$NaturalHabitat),
                         1000),
                       rep(median(modelData[
                         modelData$LandUse=="Human",]$NaturalHabitat),
                         1000))
# Set baseline TEI to be near range centre, but within 
# sampled range (10th percentile of sampled values, i.e. = 0.59)
nd$TEI_BL <- 0.59
# Set delta TEI to 0 (i.e. no climate change)
nd$TEI_delta <- 0
# Add dummy variable representing the response variable
nd$occur <- 0

# Make predictions from model
preds <- PredictGLMERMultiModel(models = finalModels,data = nd,nIters = 1000)

# Back-transform to occurrence probabilities
preds <- 1/(1+exp(-preds))

# Rescale predictions to reference level (i.e. all probabilities = 1 in 
# natural habitat)
preds <- sweep(x = preds,MARGIN = 2,STATS = preds[1,],FUN = '/')

# Calculate summary statistics across predictions
preds.median <- ((apply(X = preds,MARGIN = 1,FUN = median,na.rm=TRUE))*100)-100
preds.lower.67 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=(1/6),
                          na.rm=TRUE))*100)-100
preds.upper.67 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=(5/6),
                          na.rm=TRUE))*100)-100
preds.lower.95 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=0.025,
                          na.rm=TRUE))*100)-100
preds.upper.95 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=0.975,
                          na.rm=TRUE))*100)-100

plot(-9e99,-9e99,xlim=c(0,505),
     ylim=c(min(preds.lower.95),max(preds.upper.95)),xaxt="n",
     xlab=NA,ylab="P.occur difference (%)",bty="l")
axis(side = 1,at = c(0,125,250,375,500))
title(xlab="Duration of landscape modification (yrs)")

X.Vec <- c(nd$AgeConv[2:1001], max(nd$AgeConv[2:1001]), 
           rev(nd$AgeConv[2:1001]), min(nd$AgeConv[2:1001]))
Y.Vec <- c(preds.lower.95[2:1001], tail(preds.upper.95[2:1001], 1), 
           rev(preds.upper.95[2:1001]), (preds.lower.95[2:1001])[1])
Y.Vec2 <- c(preds.lower.67[2:1001], tail(preds.upper.67[2:1001], 1), 
            rev(preds.upper.67[2:1001]), (preds.lower.67[2:1001])[1])

polygon(x = X.Vec,y = Y.Vec,col="#009F8133",border=NA)
polygon(x = X.Vec,y = Y.Vec2,col="#009F8155",border=NA)
points(x = nd$AgeConv[2:1001],y = preds.median[2:1001],type="l",lwd=2,col="#009F81")

abline(h=0,lty=2,col="#00000055")

plot(-9e99,-9e99,xlim=c(0,505),
     ylim=c(min(preds.lower.95),max(preds.upper.95)),xaxt="n",
     xlab=NA,ylab="P.occur difference (%)",bty="l")
axis(side = 1,at = c(0,125,250,375,500))
title(xlab="Duration of landscape modification (yrs)")

X.Vec <- c(nd$AgeConv[1002:2001], max(nd$AgeConv[1002:2001]), 
           rev(nd$AgeConv[1002:2001]), min(nd$AgeConv[1002:2001]))
Y.Vec <- c(preds.lower.95[1002:2001], tail(preds.upper.95[1002:2001], 1), 
           rev(preds.upper.95[1002:2001]), (preds.lower.95[1002:2001])[1])
Y.Vec2 <- c(preds.lower.67[1002:2001], tail(preds.upper.67[1002:2001], 1), 
            rev(preds.upper.67[1002:2001]), (preds.lower.67[1002:2001])[1])

polygon(x = X.Vec,y = Y.Vec,col="#A4012233",border=NA)
polygon(x = X.Vec,y = Y.Vec2,col="#A4012255",border=NA)
points(x = nd$AgeConv[1002:2001],y = preds.median[1002:2001],
       type="l",lwd=2,col="#A40122")

abline(h=0,lty=2,col="#00000055")

invisible(dev.off())

#### Plotting (land-use intensity) ####

# Plot results
pdf(file = paste0(outDir,"ResultsPlotLandUseIntensity.pdf"),
    width = 17.5/2.54,height = 12/2.54)

par(mfrow=c(2,2))
par(tck=-0.01)
par(mar=c(2.6,2.6,0.2,0.2))
par(mgp=c(1.6,0.2,0))
par(las=1)


### Pesticide application

# Get the 95% ranges of values for pesticide application for each land-use type
pesticideRangeNatural <- quantile(
  modelData$Pesticide[modelData$LandUse=="Natural"],
  probs = c(0.025,0.975))
pesticideRangeHuman <- quantile(
  modelData$Pesticide[modelData$LandUse=="Human"],
  probs = c(0.025,0.975))

# Create data-frame with first row as reference (zero pesticide application), 
# and  then plotting ranges for natural and human land use
nd <- data.frame(Pesticide=c(0,
                             seq(from = pesticideRangeNatural[1],
                                   to = pesticideRangeNatural[2],
                                 length.out = 1000),
                             seq(from = pesticideRangeHuman[1],
                                 to = pesticideRangeHuman[2],
                                 length.out = 1000)))
# For reference row, set fertilizer to 0, then medians for each land-use type
nd$LogFertilizer <- c(0,
                      rep(median(modelData[modelData$LandUse=="Natural",]$LogFertilizer),
                          1000),
                      rep(median(modelData[modelData$LandUse=="Human",]$LogFertilizer),
                          1000))
# Reference row takes natural habitat, then 1000 values each of natural
# and human land use
nd$LandUse <- factor(c(rep("Natural",1001),rep("Human",1000)))
nd$LandUse <- relevel(x = nd$LandUse,ref = "Natural")
# For reference row, set landscape conversion age to median across whole
# dataset, then medians for each land-use type
nd$AgeConv <- c(median(modelData$AgeConv),
                  rep(median(modelData[
                    modelData$LandUse=="Natural",]$AgeConv),1000),
                  rep(median(modelData[
                    modelData$LandUse=="Human",]$AgeConv),1000))
# For reference row, set natural habitat to 100%, then
# median values for each land-use type
# Note 100% natural habitat = a value of 1,000
nd$NaturalHabitat <- c(1000,
                       rep(median(modelData[
                         modelData$LandUse=="Natural",]$NaturalHabitat),
                         1000),
                       rep(median(modelData[
                         modelData$LandUse=="Human",]$NaturalHabitat),
                         1000))
# Set baseline TEI to be near range centre, but within 
# sampled range (10th percentile of sampled values, i.e. = 0.59)
nd$TEI_BL <- 0.59
# Set delta TEI to 0 (i.e. no climate change)
nd$TEI_delta <- 0
# Add dummy variable representing the response variable
nd$occur <- 0

# Make predictions from model
preds <- PredictGLMERMultiModel(models = finalModels,data = nd,nIters = 1000)

# Back-transform to occurrence probabilities
preds <- 1/(1+exp(-preds))

# Rescale predictions to reference level (i.e. all probabilities = 1 in 
# natural habitat)
preds <- sweep(x = preds,MARGIN = 2,STATS = preds[1,],FUN = '/')

# Calculate summary statistics across predictions
preds.median <- ((apply(X = preds,MARGIN = 1,FUN = median,na.rm=TRUE))*100)-100
preds.lower.67 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=(1/6),
                          na.rm=TRUE))*100)-100
preds.upper.67 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=(5/6),
                          na.rm=TRUE))*100)-100
preds.lower.95 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=0.025,
                          na.rm=TRUE))*100)-100
preds.upper.95 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=0.975,
                          na.rm=TRUE))*100)-100

plot(-9e99,-9e99,xlim=c(0,40),
     ylim=c(min(preds.lower.95),max(preds.upper.95)),xaxt="n",
     xlab=NA,ylab="P.occur difference (%)",bty="l")
axis(side = 1,at = c(0,10,20,30,40))
title(xlab="Pesticide application (check units)")

X.Vec <- c(nd$Pesticide[2:1001], max(nd$Pesticide[2:1001]), 
           rev(nd$Pesticide[2:1001]), min(nd$Pesticide[2:1001]))
Y.Vec <- c(preds.lower.95[2:1001], tail(preds.upper.95[2:1001], 1), 
           rev(preds.upper.95[2:1001]), (preds.lower.95[2:1001])[1])
Y.Vec2 <- c(preds.lower.67[2:1001], tail(preds.upper.67[2:1001], 1), 
            rev(preds.upper.67[2:1001]), (preds.lower.67[2:1001])[1])

polygon(x = X.Vec,y = Y.Vec,col="#009F8133",border=NA)
polygon(x = X.Vec,y = Y.Vec2,col="#009F8155",border=NA)
points(x = nd$Pesticide[2:1001],y = preds.median[2:1001],type="l",lwd=2,col="#009F81")

abline(h=0,lty=2,col="#00000055")

plot(-9e99,-9e99,xlim=c(0,40),
     ylim=c(min(preds.lower.95),max(preds.upper.95)),xaxt="n",
     xlab=NA,ylab="P.occur difference (%)",bty="l")
axis(side = 1,at = c(0,10,20,30,40))
title(xlab="Pesticide application (check units)")

X.Vec <- c(nd$Pesticide[1002:2001], max(nd$Pesticide[1002:2001]), 
           rev(nd$Pesticide[1002:2001]), min(nd$Pesticide[1002:2001]))
Y.Vec <- c(preds.lower.95[1002:2001], tail(preds.upper.95[1002:2001], 1), 
           rev(preds.upper.95[1002:2001]), (preds.lower.95[1002:2001])[1])
Y.Vec2 <- c(preds.lower.67[1002:2001], tail(preds.upper.67[1002:2001], 1), 
            rev(preds.upper.67[1002:2001]), (preds.lower.67[1002:2001])[1])

polygon(x = X.Vec,y = Y.Vec,col="#A4012233",border=NA)
polygon(x = X.Vec,y = Y.Vec2,col="#A4012255",border=NA)
points(x = nd$Pesticide[1002:2001],y = preds.median[1002:2001],
       type="l",lwd=2,col="#A40122")

abline(h=0,lty=2,col="#00000055")



### Fertilizer application

# Get the 95% ranges of values for pesticide application for each land-use type
fertilizerRangeNatural <- quantile(
  modelData$LogFertilizer[modelData$LandUse=="Natural"],
  probs = c(0.025,0.975))
fertilizerRangeHuman <- quantile(
  modelData$LogFertilizer[modelData$LandUse=="Human"],
  probs = c(0.025,0.975))

# Create data-frame with first row as reference (zero pesticide application), 
# and  then plotting ranges for natural and human land use
nd <- data.frame(LogFertilizer=c(0,
                             seq(from = fertilizerRangeNatural[1],
                                 to = fertilizerRangeNatural[2],
                                 length.out = 1000),
                             seq(from = fertilizerRangeHuman[1],
                                 to = fertilizerRangeHuman[2],
                                 length.out = 1000)))
# For reference row, set pesticide to 0, then medians for each land-use type
nd$Pesticide <- c(0,
                      rep(median(modelData[modelData$LandUse=="Natural",]$Pesticide),
                          1000),
                      rep(median(modelData[modelData$LandUse=="Human",]$Pesticide),
                          1000))
# Reference row takes natural habitat, then 1000 values each of natural
# and human land use
nd$LandUse <- factor(c(rep("Natural",1001),rep("Human",1000)))
nd$LandUse <- relevel(x = nd$LandUse,ref = "Natural")
# For reference row, set landscape conversion age to median across whole
# dataset, then medians for each land-use type
nd$AgeConv <- c(median(modelData$AgeConv),
                rep(median(modelData[
                  modelData$LandUse=="Natural",]$AgeConv),1000),
                rep(median(modelData[
                  modelData$LandUse=="Human",]$AgeConv),1000))
# For reference row, set natural habitat to 100%, then
# median values for each land-use type
# Note 100% natural habitat = a value of 1,000
nd$NaturalHabitat <- c(1000,
                       rep(median(modelData[
                         modelData$LandUse=="Natural",]$NaturalHabitat),
                         1000),
                       rep(median(modelData[
                         modelData$LandUse=="Human",]$NaturalHabitat),
                         1000))
# Set baseline TEI to be near range centre, but within 
# sampled range (10th percentile of sampled values, i.e. = 0.59)
nd$TEI_BL <- 0.59
# Set delta TEI to 0 (i.e. no climate change)
nd$TEI_delta <- 0
# Add dummy variable representing the response variable
nd$occur <- 0

# Make predictions from model
preds <- PredictGLMERMultiModel(models = finalModels,data = nd,nIters = 1000)

# Back-transform to occurrence probabilities
preds <- 1/(1+exp(-preds))

# Rescale predictions to reference level (i.e. all probabilities = 1 in 
# natural habitat)
preds <- sweep(x = preds,MARGIN = 2,STATS = preds[1,],FUN = '/')

# Calculate summary statistics across predictions
preds.median <- ((apply(X = preds,MARGIN = 1,FUN = median,na.rm=TRUE))*100)-100
preds.lower.67 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=(1/6),
                          na.rm=TRUE))*100)-100
preds.upper.67 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=(5/6),
                          na.rm=TRUE))*100)-100
preds.lower.95 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=0.025,
                          na.rm=TRUE))*100)-100
preds.upper.95 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=0.975,
                          na.rm=TRUE))*100)-100

plot(-9e99,-9e99,xlim=c(0,13.5),
     ylim=c(min(preds.lower.95),max(preds.upper.95)),xaxt="n",
     xlab=NA,ylab="P.occur difference (%)",bty="l")
axis(side = 1,
     at = log(c(0,20,400,8000,700000)+1),
     labels=c(0,20,400,8000,700000))
title(xlab="Fertilizer application (check units)")

X.Vec <- c(nd$LogFertilizer[2:1001], max(nd$LogFertilizer[2:1001]), 
           rev(nd$LogFertilizer[2:1001]), min(nd$LogFertilizer[2:1001]))
Y.Vec <- c(preds.lower.95[2:1001], tail(preds.upper.95[2:1001], 1), 
           rev(preds.upper.95[2:1001]), (preds.lower.95[2:1001])[1])
Y.Vec2 <- c(preds.lower.67[2:1001], tail(preds.upper.67[2:1001], 1), 
            rev(preds.upper.67[2:1001]), (preds.lower.67[2:1001])[1])

polygon(x = X.Vec,y = Y.Vec,col="#009F8133",border=NA)
polygon(x = X.Vec,y = Y.Vec2,col="#009F8155",border=NA)
points(x = nd$LogFertilizer[2:1001],y = preds.median[2:1001],type="l",lwd=2,col="#009F81")

abline(h=0,lty=2,col="#00000055")

plot(-9e99,-9e99,xlim=c(0,13.5),
     ylim=c(min(preds.lower.95),max(preds.upper.95)),xaxt="n",
     xlab=NA,ylab="P.occur difference (%)",bty="l")
axis(side = 1,
     at = log(c(0,20,400,8000,700000)+1),
     labels=c(0,20,400,8000,700000))
title(xlab="Fertilizer application (check units)")

X.Vec <- c(nd$LogFertilizer[1002:2001], max(nd$LogFertilizer[1002:2001]), 
           rev(nd$LogFertilizer[1002:2001]), min(nd$LogFertilizer[1002:2001]))
Y.Vec <- c(preds.lower.95[1002:2001], tail(preds.upper.95[1002:2001], 1), 
           rev(preds.upper.95[1002:2001]), (preds.lower.95[1002:2001])[1])
Y.Vec2 <- c(preds.lower.67[1002:2001], tail(preds.upper.67[1002:2001], 1), 
            rev(preds.upper.67[1002:2001]), (preds.lower.67[1002:2001])[1])

polygon(x = X.Vec,y = Y.Vec,col="#A4012233",border=NA)
polygon(x = X.Vec,y = Y.Vec2,col="#A4012255",border=NA)
points(x = nd$LogFertilizer[1002:2001],y = preds.median[1002:2001],
       type="l",lwd=2,col="#A40122")

abline(h=0,lty=2,col="#00000055")

invisible(dev.off())

#### Plotting (climate interaction) ####

# Plot results
pdf(file = paste0(outDir,"ResultsPlotClimateInteraction.pdf"),
    width = 17.5/2.54,height = 12/2.54)

par(mfrow=c(2,2))
par(tck=-0.01)
par(mar=c(2.6,2.6,0.2,0.2))
par(mgp=c(1.6,0.2,0))
par(las=1)

### Baseline niche position

# Get the 95% ranges of values for pesticide application for each land-use type
tei_blRangeNatural <- quantile(
  modelData$TEI_BL[modelData$LandUse=="Natural"],
  probs = c(0.025,0.975))
tei_blRangeHuman <- quantile(
  modelData$TEI_BL[modelData$LandUse=="Human"],
  probs = c(0.025,0.975))

# Create data-frame with first row as reference (near niche centre, 
# 10th percentile of sampled values in the data, i.e. 0.59), 
# and  then plotting ranges for natural and human land use
nd <- data.frame(TEI_BL=c(0.59,
                             seq(from = tei_blRangeNatural[1],
                                 to = tei_blRangeNatural[2],
                                 length.out = 1000),
                             seq(from = tei_blRangeHuman[1],
                                 to = tei_blRangeHuman[2],
                                 length.out = 1000)))
# Reference row takes natural habitat, then 1000 values each of natural
# and human land use
nd$LandUse <- factor(c(rep("Natural",1001),rep("Human",1000)))
nd$LandUse <- relevel(x = nd$LandUse,ref = "Natural")
# For reference row, set landscape conversion age to median across whole
# dataset, then medians for each land-use type
nd$AgeConv <- c(median(modelData$AgeConv),
                rep(median(modelData[
                  modelData$LandUse=="Natural",]$AgeConv),1000),
                rep(median(modelData[
                  modelData$LandUse=="Human",]$AgeConv),1000))
# For reference row, set natural habitat to 100%, then
# median values for each land-use type
# Note 100% natural habitat = a value of 1,000
nd$NaturalHabitat <- c(1000,
                       rep(median(modelData[
                         modelData$LandUse=="Natural",]$NaturalHabitat),
                         1000),
                       rep(median(modelData[
                         modelData$LandUse=="Human",]$NaturalHabitat),
                         1000))
# For reference row, set pesticide to 0, then medians for each land-use type
nd$Pesticide <- c(0,
                  rep(median(modelData[modelData$LandUse=="Natural",]$Pesticide),
                      1000),
                  rep(median(modelData[modelData$LandUse=="Human",]$Pesticide),
                      1000))
# For reference row, set fertilizer to 0, then medians for each land-use type
nd$LogFertilizer <- c(0,
                  rep(median(modelData[modelData$LandUse=="Natural",]$LogFertilizer),
                      1000),
                  rep(median(modelData[modelData$LandUse=="Human",]$LogFertilizer),
                      1000))
# Set delta TEI to 0 (i.e. no climate change)
nd$TEI_delta <- 0
# Add dummy variable representing the response variable
nd$occur <- 0

# Make predictions from model
preds <- PredictGLMERMultiModel(models = finalModels,data = nd,nIters = 1000)

# Back-transform to occurrence probabilities
preds <- 1/(1+exp(-preds))

# Rescale predictions to reference level (i.e. all probabilities = 1 in 
# natural habitat)
preds <- sweep(x = preds,MARGIN = 2,STATS = preds[1,],FUN = '/')

# Calculate summary statistics across predictions
preds.median <- ((apply(X = preds,MARGIN = 1,FUN = median,na.rm=TRUE))*100)-100
preds.lower.67 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=(1/6),
                          na.rm=TRUE))*100)-100
preds.upper.67 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=(5/6),
                          na.rm=TRUE))*100)-100
preds.lower.95 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=0.025,
                          na.rm=TRUE))*100)-100
preds.upper.95 <- ((apply(X = preds,MARGIN = 1,FUN = quantile,probs=0.975,
                          na.rm=TRUE))*100)-100

plot(-9e99,-9e99,xlim=c(0.5,0.8),
     ylim=c(min(preds.lower.95),max(preds.upper.95)),xaxt="n",
     xlab=NA,ylab="P.occur difference (%)",bty="l")
axis(side = 1,at = c(0.5,0.6,0.7,0.8))
title(xlab="Baseline temperature niche position")

X.Vec <- c(nd$TEI_BL[2:1001], max(nd$TEI_BL[2:1001]), 
           rev(nd$TEI_BL[2:1001]), min(nd$TEI_BL[2:1001]))
Y.Vec <- c(preds.lower.95[2:1001], tail(preds.upper.95[2:1001], 1), 
           rev(preds.upper.95[2:1001]), (preds.lower.95[2:1001])[1])
Y.Vec2 <- c(preds.lower.67[2:1001], tail(preds.upper.67[2:1001], 1), 
            rev(preds.upper.67[2:1001]), (preds.lower.67[2:1001])[1])

polygon(x = X.Vec,y = Y.Vec,col="#009F8133",border=NA)
polygon(x = X.Vec,y = Y.Vec2,col="#009F8155",border=NA)
points(x = nd$TEI_BL[2:1001],y = preds.median[2:1001],type="l",lwd=2,col="#009F81")

abline(h=0,lty=2,col="#00000055")

plot(-9e99,-9e99,xlim=c(0.5,0.8),
     ylim=c(min(preds.lower.95),max(preds.upper.95)),xaxt="n",
     xlab=NA,ylab="P.occur difference (%)",bty="l")
axis(side = 1,at = c(0.5,0.6,0.7,0.8))
title(xlab="Baseline temperature niche position")

X.Vec <- c(nd$TEI_BL[1002:2001], max(nd$TEI_BL[1002:2001]), 
           rev(nd$TEI_BL[1002:2001]), min(nd$TEI_BL[1002:2001]))
Y.Vec <- c(preds.lower.95[1002:2001], tail(preds.upper.95[1002:2001], 1), 
           rev(preds.upper.95[1002:2001]), (preds.lower.95[1002:2001])[1])
Y.Vec2 <- c(preds.lower.67[1002:2001], tail(preds.upper.67[1002:2001], 1), 
            rev(preds.upper.67[1002:2001]), (preds.lower.67[1002:2001])[1])

polygon(x = X.Vec,y = Y.Vec,col="#A4012233",border=NA)
polygon(x = X.Vec,y = Y.Vec2,col="#A4012255",border=NA)
points(x = nd$TEI_BL[1002:2001],y = preds.median[1002:2001],
       type="l",lwd=2,col="#A40122")

abline(h=0,lty=2,col="#00000055")


### Niche position change

# Get the 95% ranges of values for pesticide application for each land-use type
tei_deltaRangeNatural <- quantile(
  modelData$TEI_delta[modelData$LandUse=="Natural"],
  probs = c(0.025,0.975))
tei_deltaRangeHuman <- quantile(
  modelData$TEI_delta[modelData$LandUse=="Human"],
  probs = c(0.025,0.975))

# Create data-frame for niche-centre predictions, with first row as 
# reference (no climate change), and  then plotting ranges for natural 
# and human land use
nd.centre <- data.frame(TEI_delta=c(0,
                          seq(from = tei_deltaRangeNatural[1],
                              to = tei_deltaRangeNatural[2],
                              length.out = 1000),
                          seq(from = tei_deltaRangeHuman[1],
                              to = tei_deltaRangeHuman[2],
                              length.out = 1000)))
# Reference row takes natural habitat, then 1000 values each of natural
# and human land use
nd.centre$LandUse <- factor(c(rep("Natural",1001),rep("Human",1000)))
nd.centre$LandUse <- relevel(x = nd.centre$LandUse,ref = "Natural")
# For reference row, set landscape conversion age to median across whole
# dataset, then medians for each land-use type
nd.centre$AgeConv <- c(median(modelData$AgeConv),
                rep(median(modelData[
                  modelData$LandUse=="Natural",]$AgeConv),1000),
                rep(median(modelData[
                  modelData$LandUse=="Human",]$AgeConv),1000))
# For reference row, set natural habitat to 100%, then
# median values for each land-use type
# Note 100% natural habitat = a value of 1,000
nd.centre$NaturalHabitat <- c(1000,
                       rep(median(modelData[
                         modelData$LandUse=="Natural",]$NaturalHabitat),
                         1000),
                       rep(median(modelData[
                         modelData$LandUse=="Human",]$NaturalHabitat),
                         1000))
# For reference row, set pesticide to 0, then medians for each land-use type
nd.centre$Pesticide <- c(0,
                  rep(median(modelData[modelData$LandUse=="Natural",]$Pesticide),
                      1000),
                  rep(median(modelData[modelData$LandUse=="Human",]$Pesticide),
                      1000))
# For reference row, set pesticide to 0, then medians for each land-use type
nd.centre$LogFertilizer <- c(0,
                         rep(median(modelData[modelData$LandUse=="Natural",]$LogFertilizer),
                             1000),
                         rep(median(modelData[modelData$LandUse=="Human",]$LogFertilizer),
                             1000))
# Set baseline TEI to be near range centre, but within 
# sampled range (10th percentile of sampled values, i.e. = 0.59)
nd.centre$TEI_BL <- 0.59
# Add dummy variable representing the response variable
nd.centre$occur <- 0

# Create data-frame for niche-edge predictions, with first row as 
# reference (no climate change), and  then plotting ranges for natural 
# and human land use
nd.edge <- data.frame(TEI_delta=c(0,
                                    seq(from = tei_deltaRangeNatural[1],
                                        to = tei_deltaRangeNatural[2],
                                        length.out = 1000),
                                    seq(from = tei_deltaRangeHuman[1],
                                        to = tei_deltaRangeHuman[2],
                                        length.out = 1000)))
# Reference row takes natural habitat, then 1000 values each of natural
# and human land use
nd.edge$LandUse <- factor(c(rep("Natural",1001),rep("Human",1000)))
nd.edge$LandUse <- relevel(x = nd.edge$LandUse,ref = "Natural")
# For reference row, set landscape conversion age to median across whole
# dataset, then medians for each land-use type
nd.edge$AgeConv <- c(median(modelData$AgeConv),
                       rep(median(modelData[
                         modelData$LandUse=="Natural",]$AgeConv),1000),
                       rep(median(modelData[
                         modelData$LandUse=="Human",]$AgeConv),1000))
# For reference row, set natural habitat to 100%, then
# median values for each land-use type
# Note 100% natural habitat = a value of 1,000
nd.edge$NaturalHabitat <- c(1000,
                              rep(median(modelData[
                                modelData$LandUse=="Natural",]$NaturalHabitat),
                                1000),
                              rep(median(modelData[
                                modelData$LandUse=="Human",]$NaturalHabitat),
                                1000))
# For reference row, set pesticide to 0, then medians for each land-use type
nd.edge$Pesticide <- c(0,
                         rep(median(modelData[
                           modelData$LandUse=="Natural",]$Pesticide),
                             1000),
                         rep(median(modelData[
                           modelData$LandUse=="Human",]$Pesticide),
                             1000))
# For reference row, set pesticide to 0, then medians for each land-use type
nd.edge$LogFertilizer <- c(0,
                       rep(median(modelData[
                         modelData$LandUse=="Natural",]$LogFertilizer),
                         1000),
                       rep(median(modelData[
                         modelData$LandUse=="Human",]$LogFertilizer),
                         1000))
# Set baseline niche position to near niche edge (position = 0.7)
# Represents third quartile of sampled values
nd.edge$TEI_BL <- 0.75
# Add dummy variable representing the response variable
nd.edge$occur <- 0

# Make predictions from model, for both niche-centre and niche-edge
preds.centre <- PredictGLMERMultiModel(
  models = finalModels,data = nd.centre,nIters = 1000)
preds.edge <- PredictGLMERMultiModel(
  models = finalModels,data = nd.edge,nIters = 1000)

# Back-transform to occurrence probabilities
preds.centre <- 1/(1+exp(-preds.centre))
preds.edge <- 1/(1+exp(-preds.edge))

# Rescale predictions to reference level (i.e. all probabilities = 1 in 
# natural habitat)
preds.centre <- sweep(x = preds.centre,MARGIN = 2,STATS = preds.centre[1,],FUN = '/')
preds.edge <- sweep(x = preds.edge,MARGIN = 2,STATS = preds.edge[1,],FUN = '/')

# Calculate summary statistics across predictions
preds.median.centre <- ((apply(X = preds.centre,MARGIN = 1,FUN = median,na.rm=TRUE))*100)-100
preds.lower.67.centre <- ((apply(X = preds.centre,MARGIN = 1,FUN = quantile,probs=(1/6),
                          na.rm=TRUE))*100)-100
preds.upper.67.centre <- ((apply(X = preds.centre,MARGIN = 1,FUN = quantile,probs=(5/6),
                          na.rm=TRUE))*100)-100
preds.lower.95.centre <- ((apply(X = preds.centre,MARGIN = 1,FUN = quantile,probs=0.025,
                          na.rm=TRUE))*100)-100
preds.upper.95.centre <- ((apply(X = preds.centre,MARGIN = 1,FUN = quantile,probs=0.975,
                          na.rm=TRUE))*100)-100

preds.median.edge <- ((apply(X = preds.edge,MARGIN = 1,FUN = median,na.rm=TRUE))*100)-100
preds.lower.67.edge <- ((apply(X = preds.edge,MARGIN = 1,FUN = quantile,probs=(1/6),
                                 na.rm=TRUE))*100)-100
preds.upper.67.edge <- ((apply(X = preds.edge,MARGIN = 1,FUN = quantile,probs=(5/6),
                                 na.rm=TRUE))*100)-100
preds.lower.95.edge <- ((apply(X = preds.edge,MARGIN = 1,FUN = quantile,probs=0.025,
                                 na.rm=TRUE))*100)-100
preds.upper.95.edge <- ((apply(X = preds.edge,MARGIN = 1,FUN = quantile,probs=0.975,
                                 na.rm=TRUE))*100)-100

plot(-9e99,-9e99,xlim=c(0,0.033),
     ylim=c(-100,0),xaxt="n",
     xlab=NA,ylab="P.occur difference (%)",bty="l")
axis(side = 1,at = c(0,0.01,0.02,0.03))
title(xlab="Niche position change")

X.Vec.centre <- c(nd.centre$TEI_delta[2:1001], max(nd.centre$TEI_delta[2:1001]), 
           rev(nd.centre$TEI_delta[2:1001]), min(nd.centre$TEI_delta[2:1001]))
Y.Vec.centre <- c(preds.lower.95.centre[2:1001], tail(preds.upper.95.centre[2:1001], 1), 
           rev(preds.upper.95.centre[2:1001]), (preds.lower.95.centre[2:1001])[1])
Y.Vec2.centre <- c(preds.lower.67.centre[2:1001], tail(preds.upper.67.centre[2:1001], 1), 
            rev(preds.upper.67.centre[2:1001]), (preds.lower.67.centre[2:1001])[1])

X.Vec.edge <- c(nd.edge$TEI_delta[2:1001], max(nd.edge$TEI_delta[2:1001]), 
                  rev(nd.edge$TEI_delta[2:1001]), min(nd.edge$TEI_delta[2:1001]))
Y.Vec.edge <- c(preds.lower.95.edge[2:1001], tail(preds.upper.95.edge[2:1001], 1), 
                  rev(preds.upper.95.edge[2:1001]), (preds.lower.95.edge[2:1001])[1])
Y.Vec2.edge <- c(preds.lower.67.edge[2:1001], tail(preds.upper.67.edge[2:1001], 1), 
                   rev(preds.upper.67.edge[2:1001]), (preds.lower.67.edge[2:1001])[1])

polygon(x = X.Vec.centre,y = Y.Vec2.centre,col="#009F8188",border=NA)
points(x = nd.centre$TEI_delta[2:1001],y = preds.median.centre[2:1001],
       type="l",lwd=2,col="#009F81")

polygon(x = X.Vec.edge,y = Y.Vec2.edge,col="#009F8188",border=NA)
points(x = nd.edge$TEI_delta[2:1001],y = preds.median.edge[2:1001],
       type="l",lwd=2,col="#009F81",lty=2)


abline(h=0,lty=2,col="#00000055")

plot(-9e99,-9e99,xlim=c(0,0.033),
     ylim=c(-100,0),xaxt="n",
     xlab=NA,ylab="P.occur difference (%)",bty="l")
axis(side = 1,at = c(0,0.01,0.02,0.03))
title(xlab="Niche position change")

X.Vec.centre <- c(nd.centre$TEI_delta[1002:2001], max(nd.centre$TEI_delta[1002:2001]), 
                  rev(nd.centre$TEI_delta[1002:2001]), min(nd.centre$TEI_delta[1002:2001]))
Y.Vec.centre <- c(preds.lower.95.centre[1002:2001], tail(preds.upper.95.centre[1002:2001], 1), 
                  rev(preds.upper.95.centre[1002:2001]), (preds.lower.95.centre[1002:2001])[1])
Y.Vec2.centre <- c(preds.lower.67.centre[1002:2001], tail(preds.upper.67.centre[1002:2001], 1), 
                   rev(preds.upper.67.centre[1002:2001]), (preds.lower.67.centre[1002:2001])[1])

X.Vec.edge <- c(nd.edge$TEI_delta[1002:2001], max(nd.edge$TEI_delta[1002:2001]), 
                rev(nd.edge$TEI_delta[1002:2001]), min(nd.edge$TEI_delta[1002:2001]))
Y.Vec.edge <- c(preds.lower.95.edge[1002:2001], tail(preds.upper.95.edge[1002:2001], 1), 
                rev(preds.upper.95.edge[1002:2001]), (preds.lower.95.edge[1002:2001])[1])
Y.Vec2.edge <- c(preds.lower.67.edge[1002:2001], tail(preds.upper.67.edge[1002:2001], 1), 
                 rev(preds.upper.67.edge[1002:2001]), (preds.lower.67.edge[1002:2001])[1])

polygon(x = X.Vec.centre,y = Y.Vec2.centre,col="#A4012288",border=NA)
points(x = nd.centre$TEI_delta[1002:2001],y = preds.median.centre[1002:2001],
       type="l",lwd=2,col="#A40122")

polygon(x = X.Vec.edge,y = Y.Vec2.edge,col="#A4012288",border=NA)
points(x = nd.edge$TEI_delta[1002:2001],y = preds.median.edge[1002:2001],
       type="l",lwd=2,col="#A40122",lty=2)

abline(h=0,lty=2,col="#00000055")

invisible(dev.off())

t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()
