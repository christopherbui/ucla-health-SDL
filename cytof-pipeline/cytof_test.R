# load required packages
library(CATALYST)
library(cowplot)
library(flowCore)
library(ggplot2)
library(SingleCellExperiment)

# load "raw_data" dataset from CATALYST
data("raw_data")

# INSPECT ----------------------------------------
class(raw_data) # flowSet containing 2 flowFrames
sampleNames(raw_data)

# exprs() shows the expression matrix for flowFrame object
# raw_data_1.fcs with 2500 cells & 67 features
ff_1 <- raw_data[[1]]
exprs(ff_1)

# view column names, and shape
colnames(ff_1)
dim(ff_1)
# ------------------------------------------------

# construct SingelCellExperiment object
sce <- prepData(raw_data)
dim(sce)

# view number of events per sample/flowFrame
table(sce$sample_id)

# view non-isotope features
names(int_colData(sce))

