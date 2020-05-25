library(reticulate)
# Using a specific python binary
use_python("C:/Users/acsala/AppData/Local/Programs/Python/Python37/python.exe", required = TRUE)

# https://bioconductor.org/packages/release/bioc/html/MOFA.html
library(MOFAdata)
library(MOFA)
#data("CLL_data")

#MOFA with simulated data ####
#https://bioconductor.org/packages/release/bioc/vignettes/MOFA/inst/doc/MOFA_example_simulated.html
set.seed(1234)
data <- makeExampleData()
MOFAobject <- createMOFAobject(data)

MOFAobject

TrainOptions <- getDefaultTrainOptions()
ModelOptions <- getDefaultModelOptions(MOFAobject)
DataOptions <- getDefaultDataOptions()

TrainOptions$DropFactorThreshold <- 0.01


n_inits <- 3
MOFAlist <- lapply(seq_len(n_inits), function(it) {
  
  TrainOptions$seed <- 2018 + it
  
  MOFAobject <- prepareMOFA(
    MOFAobject, 
    DataOptions = DataOptions,
    ModelOptions = ModelOptions,
    TrainOptions = TrainOptions
  )
  
  runMOFA(MOFAobject)
})

compareModels(MOFAlist)

compareFactors(MOFAlist)

MOFAobject <- selectModel(MOFAlist, plotit = FALSE)
MOFAobject

plotVarianceExplained(MOFAobject)

sessionInfo()

#MOFA on 200 patient data ####
# https://rdrr.io/github/bioFAM/MOFA/f/vignettes/MOFA_example_CLL.Rmd
library(MultiAssayExperiment)
library(MOFA)
library(MOFAdata)

#Data loading option 1####
data("CLL_data")
MOFAobject <- createMOFAobject(CLL_data)
MOFAobject

#Data loading option 2####
# Load data
# import list with mRNA, Methylation, Drug Response and Mutation data. 
data("CLL_data") 

# check dimensionalities, samples are columns, features are rows
lapply(CLL_data, dim) 

# Load sample metadata: Sex and Diagnosis
data("CLL_covariates")
head(CLL_covariates)

# Create MultiAssayExperiment object 
mae_CLL <- MultiAssayExperiment(
  experiments = CLL_data, 
  colData = CLL_covariates
)

# Build the MOFA object
MOFAobject <- createMOFAobject(mae_CLL)
MOFAobject


plotDataOverview(MOFAobject)


#fit the model####

DataOptions <- getDefaultDataOptions()
DataOptions 


ModelOptions <- getDefaultModelOptions(MOFAobject)
ModelOptions$numFactors <- 5
ModelOptions


TrainOptions <- getDefaultTrainOptions()

# Automatically drop factors that explain less than 2% of variance in all omics
TrainOptions$DropFactorThreshold <- 0.02

TrainOptions$seed <- 2017

TrainOptions$maxiter <- 50

TrainOptions

MOFAobject <- prepareMOFA(
  MOFAobject, 
  DataOptions = DataOptions,
  ModelOptions = ModelOptions,
  TrainOptions = TrainOptions
)


#MOFAobject <- runMOFA(MOFAobject)
# Loading an existing trained model
filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAdata")

MOFAobject <- loadModel(filepath, MOFAobject)
MOFAobject

#*************************************************************************************************
#STEP 2: Analyse a trained MOFA model####
#*************************************************************************************************

#*************************************************************************************************
# #Part 1: Disentangling the heterogeneity####
#*************************************************************************************************
# Calculation of variance explained by each factor in each view. 
# This is probably the most important plot that MOFA generates, 
# as it summarises the entire heterogeneity of the dataset in a 
# single figure. Here we can see in which view a factor explains 
# variation which can guide further characterisation of the 
# factors by investigating the weights in those views.

# Calculate the variance explained (R2) per factor in each view 
r2 <- calculateVarianceExplained(MOFAobject)
r2$R2Total

# Variance explained by each factor in each view
head(r2$R2PerFactor)

format(r2$R2PerFactor[,3]*100, scientific = F, digits = 2)

sum((r2$R2PerFactor[,3]))

# Plot it
plotVarianceExplained(MOFAobject)


#*************************************************************************************************
#Part 2: Characterisation of individual factors####
#*************************************************************************************************

plotWeightsHeatmap(
  MOFAobject, 
  view = "Mutations", 
  factors = 1:5,
  show_colnames = FALSE
)

plotWeights(
  MOFAobject, 
  view = "Mutations", 
  factor = 2, 
  nfeatures = 10
)

plotWeights(
  MOFAobject, 
  view = "Mutations", 
  factor = 1, 
  nfeatures = 5,
  manual = list(c("BRAF"),c("MED12")),
  color_manual = c("red","blue")
  
)

plotTopWeights(
  MOFAobject, 
  view="Mutations", 
  factor=2, 
  nfeatures = 20
)

plotTopWeights(
  MOFAobject, 
  view = "mRNA", 
  factor = 1
)

plotDataHeatmap(
  MOFAobject, 
  view = "mRNA", 
  factor = 1, 
  features = 20, 
  show_rownames = FALSE
)

#*************************************************************************************************
#Feature set enrichment analysis in the active views####
#*************************************************************************************************

# Load reactome annotations
data("reactomeGS") # binary matrix with feature sets in rows and features in columns
reactomeGS[1:10,1:10]

# perform enrichment analysis
gsea <- runEnrichmentAnalysis(
  MOFAobject,
  view = "mRNA",
  feature.sets = reactomeGS,
  alpha = 0.01
)


plotEnrichmentBars(gsea, alpha=0.01)

interestingFactors <- 4:5

fseaplots <- lapply(interestingFactors, function(factor) {
  plotEnrichment(
    MOFAobject,
    gsea,
    factor = factor,
    alpha = 0.01,
    max.pathways = 10 # The top number of pathways to display
  )
})


cowplot::plot_grid(fseaplots[[1]], fseaplots[[2]],
                   ncol = 1, labels = paste("Factor", interestingFactors))


#Ordination of samples by factors to reveal clusters and gradients in the sample space####
plotFactorScatter(
  MOFAobject,
  factors = 1:2,
  color_by = "IGHV",      # color by the IGHV values that are part of the training data
  shape_by = "trisomy12"  # shape by the trisomy12 values that are part of the training data
)

plotFactorScatters(
  MOFAobject,
  factors = 1:3,
  color_by = "IGHV"
)

plotFactorBeeswarm(
  MOFAobject,
  factors = 1,
  color_by = "IGHV"
)

#*************************************************************************************************
#Customized analysis####
#*************************************************************************************************

MOFAweights <- getWeights(
  MOFAobject, 
  views = "all", 
  factors = "all", 
  as.data.frame = TRUE    # if TRUE, it outputs a long dataframe format. If FALSE, it outputs a wide matrix format
)
head(MOFAweights)


MOFAfactors <- getFactors(
  MOFAobject, 
  factors = c(1,2),
  as.data.frame = FALSE   # if TRUE, it outputs a long dataframe format. If FALSE, it outputs a wide matrix format
)
head(MOFAfactors)

#*************************************************************************************************
#Imputation of missing observations####
#*************************************************************************************************

MOFAobject <- impute(MOFAobject)
imputedDrugs <- getImputedData(MOFAobject, view="Drugs")[[1]]


# training data (incl. missing values)
pheatmap::pheatmap(drugData4Training[1:40,1:20],
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   show_rownames = FALSE, show_colnames = FALSE)

# imputed data
pheatmap::pheatmap(imputedDrugs[1:40,1:20],
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   show_rownames = FALSE, show_colnames = FALSE)




#*************************************************************************************************
#*************************************************************************************************
#*************************************************************************************************
#*************************************************************************************************
#*************************************************************************************************
#*************************************************************************************************
#*************************************************************************************************

#gedoe met python #later solved ####
# library(reticulate)
# #use_python("C:/Users/acsala/AppData/Local/Continuum/anaconda3/envs/r-reticulate/python.exe")
# #use_virtualenv("r-reticulate")
# 
# mofapy <- import("mofapy")
# mofapy
# 
# #py_config()
# #py_install("mofapy", envname = "r-reticulate", method="auto")
# #Sys.which("python")
# 
# 
# devtools::install_github("bioFAM/MOFAdata", build_opts = c("--no-resave-data"))
# devtools::install_github("bioFAM/MOFA", build_opts = c("--no-resave-data"))

library(reticulate)

# Using a specific python binary
#use_python("C:/Users/acsala/AppData/Local/Continuum/anaconda3/python.exe", required = TRUE)
use_python("C:/Users/acsala/AppData/Local/Programs/Python/Python37/python.exe", required = TRUE)
#use_condaenv(conda = "C:/Users/acsala/AppData/Local/Continuum/anaconda3/")

#py_install("mofapy", envname = "r-reticulate", method="auto", conda = "C:/Users/acsala/AppData/Local/Continuum/anaconda3/")

# Using a conda enviroment called "r-reticulate"
#use_condaenv("r-reticulate", required = TRUE)

# install.packages("xml2",
#                  configure.vars = c("INCLUDE_DIR=/trinity/home/cscala/anaconda3/include/libxml2"))
# 
# 
# install.packages("xml2",
#                  configure.vars = c("INCLUDE_DIR=/trinity/home/cscala/anaconda3"))
# 
# library(withr)
# #for linux if C compilation fails
# with_makevars(c(PKG_CFLAGS = "-std=c99"),
#               install.packages("devtools"),
#               assignment = "+=")

