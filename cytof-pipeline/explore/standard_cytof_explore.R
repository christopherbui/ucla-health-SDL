library(stringr)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(pals)
library(scales)
library(ggplot2)
library(pheatmap)
library(pandoc)   #required by flowAI

library(flowCore)
library(flowAI)
library(flowCut)
library(cytutils)   #clauclting aof/change fcs name


library(flowDensity)
library(CytoNorm)
library(CATALYST)
library(SingleCellExperiment)
library(SummarizedExperiment)

# Adjust path to scripts accordingly
source(paste("U:/cdbui/MPE_Cytof/Rscript/", "Rybakowska_cytof_function.R", sep = ""))
source("C:/Users/cdbui/Documents/GitHub/latent-cell-space/cytof-pipeline/B/MPE_Cytof/Rscript/Rybakowska_cytof_function_LT.R")


################################################################################################
# BASIC EXPLORATION
################################################################################################


ff_path <- "U:/cdbui/MPE_Cytof/CYTOF_data/Cleaned/Myeloid/Myeloid_PID00-MPE_01_batch1_cleaned.fcs"

# flowframe
ff <- read.FCS(ff_path, transformation = FALSE)

# slots of flowframe S4 object
# ff@expr
# ff@parameters$
# ff@description

# channel names & marker info
pData(parameters(ff))

# expression matrix
expr <- exprs(ff)

# channel names of expression matrix
colnames(ff)
# colnames(exprs(ff)) is the same as above

# access slots
ff@exprs
ff@parameters$desc # metadata column from pData(parameters(ff)) i.e. name, desc, range, minRange, maxRange
ff@description$FIL

dim(ff)
nrow(ff)
ncol(ff)

ff@description

################################################################################################
# BAD FCS FILE - OUTLIER DETECTION STEP - ANALYSIS
################################################################################################

ff_bad_path <- "U:/cdbui/MPE_Cytof/CYTOF_data/Cleaned/BAD_FOR_OUTLIER_STEP/Cytokines_BAD/Cytokines_Ref-PBMC_01_batch5_cleaned.fcs"

# load fcs
ff <- read.FCS(ff_bad_path, transformation = FALSE)

# get expression matrix
expr <- exprs(ff)

# initiate vectors to hold bad channel names
channels_all_null <- c()
constant_channels <- c()

channels_with_any_null <- c()

# loop over channels
for (channel in colnames(expr)) {
    values <- expr[, channel]

    # check if channel is all null or constant
    if (all(is.na(values))) {
        channels_all_null <- c(channels_all_null, channel)
    } else if (sd(values, na.rm = TRUE) == 0) {
        constant_channels <- c(constant_channels, channel)
    }

    # check if any nulls in channel
    if (any(is.na(expr[, channel]))) {
        channels_with_any_null <- c(channels_with_any_null, channel)
    }

}

# output results
cat("Channels with all NA values:\n")
print(channels_all_null)

cat("\nChannels with constant values:\n")
print(constant_channels)

cat("Channels with at least one NA:\n")
print(channels_with_any_null)



################################################################################################
# GATING WALKTHROUGH
################################################################################################

hard_cutoff = 4
tinypeak_removal1 = 0.8
tinypeak_removal2 = 0.8
alpha1 = 0.05
alpha2 =0.1


ff_path <- "U:/cdbui/MPE_Cytof/CYTOF_data/Cleaned/Myeloid/Myeloid_PID22-MPE_01_batch6_cleaned.fcs"

# signal cleaned ff
ff <- read.FCS(ff_path, transform = FALSE)

# get file name
file_name <- ff@description$FIL

# transform
ff_t <- flowCore::transform(ff, flowCore::transformList(colnames(ff)[grep("Di", colnames(ff))], CytoNorm::cytofTransform))

# selection mask
selection <- matrix(TRUE,
                    nrow = nrow(ff),
                    ncol = 1,
                    dimnames = list(NULL, c("intact")))

# threshold holder
tr <- list()

for (m in c("Ir193Di", "Ir191Di")) {
    tr[[m]] <- c(flowDensity::deGate(ff_t,
                                      m,
                                      tinypeak.removal = tinypeak_removal1,
                                      upper = FALSE,
                                      use.upper = TRUE,
                                      alpha = alpha1,
                                      verbose = FALSE,
                                      count.lim = 3),
                flowDensity::deGate(ff_t,
                                      m,
                                      tinypeak.removal = tinypeak_removal2,
                                      upper = TRUE,
                                      use.upper = TRUE,
                                      alpha = alpha2,
                                      verbose = FALSE,
                                      count.lim = 3))
}

tr_thres <- tr

for (m in c("Ir193Di", "Ir191Di")) {
    # if deGate() low threshold is bad (low variance in marker intensity)
    if (tr[[m]][1] < 0.1 | is.na(tr[[m]][1])) {
        # if deGate() high threshold is > hard cut off, use the hard_cutoff as selected (lower bound and onward)
        if (tr[[m]][2] > hard_cutoff) {
            tr_thres[[m]][3] <- hard_cutoff
            selection[ff_t@exprs[, m] < tr_thres[[m]][3], "intact"] <- FALSE
        # if deGate() high threshold is < hard cut off,
        # use deGate() high threshold as the selected (lower bound and onward)
        } else {
            tr_thres[[m]][3] <- tr[[m]][2]
            selection[ff_t@exprs[, m] < tr[[m]][2], "intact"] <- FALSE
        }
    # if deGate() low threshold is good,
    # use it as selected (lowwer bound and onward)
    } else {
        tr_thres[[m]][3] <- tr[[m]][1]
        selection[ff_t@exprs[, m] < tr[[m]][1], "intact"] <- FALSE
    }
}

percentage <- (sum(selection)/length(selection))*100
flowDensity::plotDens(ff_t, c("Ir193Di", "Ir191Di"), 
                    main = paste0(basename(file_name)," ( ", format(round(percentage, 2), 
                                                                    nsmall = 2), "% )"))

abline(h = c(tr_thres[["Ir191Di"]]))
abline(v = c(tr_thres[["Ir193Di"]]))
points(ff_t@exprs[!selection[,"intact"], c("Ir193Di", "Ir191Di")], pch = ".")

tr_thres <- as.data.frame(tr_thres)
tr_thres$info <- c("low", "high", "selected")

# tr_thres$file_name <- rep()

fout4thres <- gsub(".fcs", "_intact_thres.txt", file_name)
write.table(tr_thres, fout4thres, sep="\t", quote = FALSE, row.names = FALSE)

ff <- ff[selection[,"intact"], ]