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

source(paste("U:/cdbui/MPE_Cytof/Rscript/", "Rybakowska_cytof_function.R", sep = ""))
source("C:/Users/cdbui/Documents/GitHub/ucla-health-SDL/cytof-pipeline/MPE_Cytof/Rscript/Rybakowska_cytof_function_LT.R")

# progress bar
library(progress)

reloadProgressBar <- function(iterations) {
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = iterations, clear = FALSE, width = 60
  )
  return(pb)
}



# set directories
workFolder <- paste("U:/cdbui/", "MPE_Cytof", sep = "")
setwd(workFolder)

# set panel
selPanel <- "TBNK"

bead_norm_dir <- file.path(workFolder, "CYTOF_data", "BeadNorm", selPanel)
clean_dir <- file.path(workFolder, "CYTOF_data", "Cleaned", selPanel)
gate_dir <- file.path(workFolder, "CYTOF_data", "Gated", selPanel)

# concatenate directories
stage_dirs <- c(bead_norm_dir, clean_dir, gate_dir)
stage_labels <- c("beadNorm", "cleaned", "gated")

# read in sample info to get corename
fin_info <- file.path("Ranalysis", "mpe_cytof_sampleInfo_022625.txt")
allSampleInfo <- read.delim(fin_info, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
allSampleInfo <- allSampleInfo[, -1]

# get corenames for selected panel
tmp <- dplyr::filter(allSampleInfo, panel_id == selPanel)
corenames <- tmp$corename



check_num_events <- function(corenames, dirs, stage_labels) {
    master_df <- data.frame()

    pb <- reloadProgressBar(length(corenames))
    for (corename in corenames) {
        row <- list(corename = corename)

        for (i in seq_along(dirs)) {
            dir <- dirs[i]
            stage <- stage_labels[i]

            # look for file matching name
            fcs_file <- list.files(dir, pattern = paste0(corename, ".*\\.fcs$"), full.names = TRUE)
            print(fcs_file)
            if (length(fcs_file) == 1) {
                ff <- flowCore::read.FCS(fcs_file, transformation = FALSE)
                row[[paste0(stage, "_events")]] <- nrow(ff)
            } else {
                row[[paste0(stage, "_events")]] <- NULL
            }
        }
        master_df <- dplyr::bind_rows(master_df, row)

        pb$tick()
    }
    return(master_df)
}

master_df <- check_num_events(corenames, stage_dirs, stage_labels)

# write to csv
write.csv(master_df, file = paste0("event_counts_", selPanel, ".csv"), row.names = FALSE)

dim(master_df)
