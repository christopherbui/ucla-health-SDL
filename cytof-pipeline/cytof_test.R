# load required packages
library(CATALYST)
library(cowplot)
library(flowCore)
library(ggplot2)
library(SingleCellExperiment)

# load raw_data dataset from CATALYST
data("raw_data")
class(raw_data)
sampleNames(raw_data)
exprs(raw_data[[1]])

# inspect
str(raw_data)

raw_data[,1]

