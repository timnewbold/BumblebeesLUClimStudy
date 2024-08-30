#### Setup ####

# Set input directories
inDir <- "0_data/"
dataDir <- "2_PrepareDiversityData/"
modelsDir <- "3_RunModels/"
modelsDirRobustness <- "6_ModelRobustnessChecks/"

outDir <- "0_PrepareBibTex/"

# Create output directory, if it doesn't already exist
if(!dir.exists(outDir)) dir.create(outDir)

# Create log file
sink(paste0(outDir,"log.txt"))

# Start timer
t.start <- Sys.time()

print(t.start)

suppressMessages(suppressWarnings(library(curl)))

# Print session information to log file
sessionInfo()

# Load the compiled species-level data
diversity <- readRDS(paste0(dataDir,"diversity_data.rds"))

allData <- readRDS(paste0(modelsDir,"ModelData.Rds"))

sites <- readRDS(paste0(modelsDir,"AnalysisSites.rds"))

# Set input directory
inDir <- "3_RunModels/"

# Log-transform elevation, which is very skewed
diversity$LogElevation <- log(diversity$Elevation+2)

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
                          'TEI_delta','LogElevation',
                          'Source_ID','SS','SSBS',
                          'Taxon_name_entered','Longitude','Latitude',
                          'Best_guess_binomial','Country',
                          'LogPestToxLow',
                          'NaturalHabitat2k','AgeConv')]

# Remove rows with any NA values
modelData <- na.omit(modelData)

# Read full list of data sources for the PREDICTS database
data.sources <- read.csv(paste0(inDir,"resource.csv"))

# Select data sources included in this study
data.sources <- data.sources[(data.sources$Source_ID %in% modelData$Source_ID),]

# Fix a mistake in one of the DOIs
data.sources$DOI[(data.sources$DOI=="1007/s10980-012-9820-6")] <- "10.1007/s10980-012-9820-6"

h <- new_handle()
handle_setheaders(h, "accept" = "application/x-bibtex")

invisible(sapply(X = data.sources$DOI,FUN = function(doi){
  
  if (doi != ""){
    url <- paste0("https://doi.org/", doi)
    print(paste0("url: ", url))
    bib <- readLines(curl(url,handle = h))
    write(x = bib,file = "data_references.bib",append = TRUE)
  } else {
    
  }

}))

