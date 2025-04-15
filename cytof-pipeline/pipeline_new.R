source("functions.R")

library(CATALYST)
library(flowCore)
library(SingleCellExperiment)

library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# main data directory
data_dir <- file.path(getwd(), "data")

# raw & normalized data directory
raw_data_dir <- file.path(data_dir, "RawFiles")
bead_norm_dir <- file.path(data_dir, "BeadNorm")


# -------------------------------------
# NORMALIZE
# -------------------------------------

# get data file names
files <- list.files(file.path(raw_data_dir),
                    pattern = ".fcs$",
                    full.names = TRUE)

for (file in files) {
  # load fcs file
  ff <- read.FCS(file,
                 transformation = FALSE)
  
  # construct SCE
  sce <- prepData(ff)
  
  # normalize
  ff_norm <- normCytof(sce,
                   beads = "dvs",
                   k = 50,
                   assays = c("counts", "exprs"),
                   overwrite = FALSE,
                   remove_beads = TRUE)
  
  # convert sce to .fcs and save 
  fcs <- sce2fcs(ff_norm$data)
  write.FCS(fcs,
            filename = file.path(bead_norm_dir,
                                 gsub(".fcs",
                                      "_norm.fcs",
                                      basename(file), ignore.case = TRUE)))
}



# -------------------------------------
# SIGNAL CLEANING
# -------------------------------------

# set cleaned signal data directory
cleaned_dir <- file.path(data_dir, "Cleaned")
if(!dir.exists(cleaned_dir)) dir.create(cleaned_dir)

# get files to be cleaned
files <- list.files(bead_norm_dir,
                    ignore.case = TRUE,
                    pattern = ".fcs$",
                    full.names = TRUE)

# clean file by file
for (file in files) {
  # read fcs file
  ff <- read.FCS(filename = file,
                 transformation = FALSE)
  
  # clean flow rate
  ff <- clean_flow_rate(flow_frame = ff,
                        out_dir = cleaned_dir,
                        to_plot = TRUE,
                        data_type = "MC")
  
  # clean signal
  ff <- clean_signal(flow_frame = ff,
                     to_plot = "All",
                     out_dir = cleaned_dir,
                     Segment = 500,
                     arcsine_transform = TRUE,
                     data_type = "MC",
                     non_used_bead_ch = NULL)
  
  # write fcs file
  write.FCS(ff,
            file = file.path(cleaned_dir,
                             gsub("_norm", "_cleaned", basename(file))))
}


# -------------------------------------
# IDENTIFY ALIQUOT OUTLIERS
# -------------------------------------

# set outlier output directory
qc_dir <- file.path(data_dir, "QualityControl")
if(!dir.exists(qc_dir)) dir.create(qc_dir)

# get data file names
files <- list.files(cleaned_dir,
                    pattern = ".fcs$",
                    full.names = TRUE)

# get batch_id for each file
# file_batch_id <- stringr::str_match(basename(files),
#                                     "(day[0-9]*).*.fcs")[, 2]


# Save original function
original_rtsne <- get("Rtsne", envir = asNamespace("Rtsne"))

# Unlock the binding so we can overwrite it
unlockBinding("Rtsne", asNamespace("Rtsne"))

# Overwrite the function within the namespace
assign("Rtsne", function(..., check_duplicates = FALSE) {
  original_rtsne(..., check_duplicates = check_duplicates)
}, envir = asNamespace("Rtsne"))




file_quality_check(fcs_files = files, 
                   file_batch_id = NULL, 
                   out_dir = qc_dir,
                   phenotyping_markers = c("Ir","CD", "HLA", "IgD", "Pt"),
                   # phenotyping_markers = c("Ir", "Pt"),
                   arcsine_transform = TRUE, 
                   nClus = 10,
                   sd = 3)


# -------------------------------------
# DEBARCODING
# -------------------------------------

# ------------------------------------------------------------------------------
# AGGREGATION
# ------------------------------------------------------------------------------
































