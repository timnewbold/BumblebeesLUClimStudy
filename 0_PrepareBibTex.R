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

# Add DOI for Hanley et al. (2011), which covers 2 unpublished datasets, and then
# merge the two unpublished datasets by removing the second
data.sources$DOI[(data.sources$BibTeX_reference=="Hanley2011")] <- "10.1111/j.1600-0706.2011.19233.x"
data.sources <- data.sources[-which(data.sources$BibTeX_reference=="Hanley2005"),]

h <- new_handle()
handle_setheaders(h, "accept" = "application/x-bibtex")

labels <- (sapply(X = data.sources$DOI,FUN = function(doi){
  
  if (doi != ""){
    url <- paste0("https://doi.org/", doi)
    print(paste0("url: ", url))
    bib <- readLines(curl(url,handle = h))
    write(x = bib,file = "data_references.bib",append = TRUE)
    label <- paste0("@",strsplit(strsplit(bib,'\\{')[[1]][2],',')[[1]][1])
  } else {
    label <- NA
  }
  
  return(label)

}))

cat("\n")

# Manually create BibTex entries for references without DOIs
bib <- " @phdthesis{Fowler2014, title={An investigation into bee assemblage change along an urban-rural gradient}, author={Fowler, Robert E.}, school={University of Birmingham}, year={2014}}"
write(x = bib,file = "data_references.bib",append = TRUE)
labels <- c(labels,"@Fowler2014")

bib <- " @article{Meyer2007, title={Patch size and landscape effects on pollinators and seed set of the horseshoe vetch, Hippocrepis comosa, in an agricultural landscape of central Europe}, volume={30}, ISSN={0171-0877}, number={2}, journal={Entomologia Generalis}, author={Meyer, Birgit and Gaebele, Volker and Steffan-Dewenter, Ingolf D.}, year={2007}, pages={173-185}}"
write(x = bib,file = "data_references.bib",append = TRUE)
labels <- c(labels,"@Meyer2007")

bib <- " @article{Quaranta2004, title={Wild bees in agroecosystems and semi-natural landscapes. 1997-2000 collection period in Italy}, volume={57}, ISSN={1721-8861}, number={1}, journal={Bulletin of Insectology}, author={Quaranta, Marino and Ambroselli, Sabrina and Barro, Paola and Bella, Salvatore and Carini, Alfredo, and Celli, Giorgio and Cogoi, Piero and Comba, Livio and Comoli, Riccardo and Felicioli, Antonio and Floris, Ignazio and Intoppa, Francesco and Longo, Santi and Maini, Stefano and Manino, Aulo and Mazzeo, Gaetana and Medrzycki, Piotr and Nardi, Emanuela and Niccolini, Lucia and Palmieri, Nicola and Patetta, Augusto and Piatti, Claudia and Piazza, Maria Gioia and Pinzauti, Mauro and Porporato, Marco and Porrini, Claudio and Ricciardelli D'Albore, Giancarlo and Romagnoli, Francesco and Ruiu, Luca and Satta, Alberto and Zandigiacomo, Pietro}, year={2004}, pages={11-61}}"
write(x = bib,file = "data_references.bib",append = TRUE)
labels <- c(labels,"@Quaranta2004")

cat(paste(labels,collapse=", "))

cat("\n")
cat("\n")

#### Finish ####

t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()