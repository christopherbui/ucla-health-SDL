source("functions.R")

library(CATALYST)
library(flowCore)
library(ggplot2)

ff_name <- "C:/Users/cdbui/Desktop/PBMC8_30min_patient4_BCR-XL.fcs"
ff_mod_name <- "C:/Users/cdbui/Documents/GitHub/latent-cell-space/cytof-pipeline/data/RawFiles/PBMC8_30min_patient4_BCR-XL_mod.fcs"

ff <- read.FCS(ff_name, transformation = FALSE)
ff_mod <- read.FCS(ff_mod_name, transformation = FALSE)

pData(parameters(ff))
pData(parameters(ff_mod))

ff_col <- colnames(ff)
ff_exp <- exprs(ff)

ff_col

ff_mod_col <- colnames(ff_mod)
ff_mod_exp <- exprs(ff_mod)

ff_mod_col


ff_mod_exp[, "165Di"]
