source("functions.R")

library(CATALYST)
library(flowCore)
library(SingleCellExperiment)

library(ggplot2)
library(ggpubr)
library(flowDensity)
library(RColorBrewer)
library(readxl)



# -------------------------------------
# SETUP
# -------------------------------------

# main data directory
data_dir <- file.path(getwd(), "data")

# raw & normalized data directory
raw_data_dir <- file.path(data_dir, "RawFiles")
if (!dir.exists(raw_data_dir)) dir.create(raw_data_dir)



# -------------------------------------
# NORMALIZE
# -------------------------------------

# set bead normalization directory
bead_norm_dir <- file.path(data_dir, "BeadNorm")
if (!dir.exists(bead_norm_dir)) dir.create(bead_norm_dir)

# get data file names
files <- list.files(file.path(raw_data_dir),
                    pattern = ".fcs$",
                    full.names = TRUE)

panel <- read.csv(file.path(raw_data_dir, "panel.csv"))

# ff <- read.FCS(files[1], transformation = FALSE)
# exprs(ff)[, "Time"]
# colnames(ff)
# sce <- prepData(ff)
# assayNames(sce)
# 
# colnames(ff)

for (file in files) {
  # load fcs file
  ff <- read.FCS(file,
                 transformation = FALSE)
  
  # construct SCE
  sce <- prepData(ff,
                  panel = panel)
  
  # normalize
  ff_norm <- normCytof(sce,
                   beads = "dvs",
                   # beads = c(139, 141, 146, 159, 165, 169, 175),
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
# VISUALIZE NORMALIZED FILES
# -------------------------------------

# NOTE: HAS OUTPUT, BUT NOT FUNCTIONAL


# files before normalization
files_b <- list.files(raw_data_dir,
                      pattern = ".fcs$",
                      ignore.case = TRUE,
                      full.names = TRUE)

# files after normalization
files_a <- list.files(bead_norm_dir,
                      pattern = ".fcs$",
                      ignore.case = TRUE,
                      full.names = TRUE)

plot_marker_quantiles(files_after_norm = files_a, 
                      files_before_norm = files_b, 
                      batch_pattern = "(.*)", 
                      arcsine_transform = TRUE,
                      remove_beads = TRUE,
                      bead_channel = "140", 
                      uncommon_prefix = NULL,
                      markers_to_plot = c("CD", "HLA", "IgD", "IL", "TNF",
                                          "TGF", "GR", "IFNa"),
                      manual_colors = c("darkorchid4", "darkorange", "darkgreen"),
                      out_dir = bead_norm_dir)



# -------------------------------------
# SIGNAL CLEANING
# -------------------------------------

# set cleaned signal data directory
cleaned_dir <- file.path(data_dir, "Cleaned")
if (!dir.exists(cleaned_dir)) dir.create(cleaned_dir)

# get files to be cleaned
files <- list.files(bead_norm_dir,
                    ignore.case = TRUE,
                    pattern = ".fcs$",
                    full.names = TRUE)

# ff <- read.FCS(files[1], transformation = FALSE)
# sce <- prepData(ff)
# colnames(ff)
# files


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
if (!dir.exists(qc_dir)) dir.create(qc_dir)

# get data file names
files <- list.files(cleaned_dir,
                    pattern = ".fcs$",
                    full.names = TRUE)

# get batch_id for each file
# file_batch_id <- stringr::str_match(basename(files),
#                                     "(day[0-9]*).*.fcs")[, 2]


# change Rtsne() function argument: check_duplicates = FALSE
# by default, it thinks fcs files have duplicates.
# I confirmed fcs has no duplicates.

# Save original function
# original_rtsne <- get("Rtsne", envir = asNamespace("Rtsne"))
# 
# # Unlock the binding so we can overwrite it
# unlockBinding("Rtsne", asNamespace("Rtsne"))
# 
# # Overwrite the function within the namespace
# assign("Rtsne", function(..., check_duplicates = FALSE) {
#   original_rtsne(..., check_duplicates = check_duplicates)
# }, envir = asNamespace("Rtsne"))


file_quality_check(fcs_files = files, 
                   file_batch_id = NULL, 
                   out_dir = qc_dir,
                   phenotyping_markers = c("Ir","CD", "HLA", "IgD", "Pt"),
                   # phenotyping_markers = c("Ir", "Pt"),
                   arcsine_transform = TRUE, 
                   nClus = 10,
                   sd = 3)
# Error in Rtsne.default(dimred_data, ...) : 
#   Remove duplicates before running TSNE.


# -------------------------------------
# DEBARCODING (SKIPPED)
# -------------------------------------



# -------------------------------------
# AGGREGATION
# -------------------------------------

# set aggregated output directory
aggregate_dir <- file.path(data_dir, "Aggregated")
if (!dir.exists(aggregate_dir)) dir.create(aggregate_dir)


# get files to be aggregated
files <- list.files(cleaned_dir,
                    pattern = ".fcs$",
                    full.names = TRUE)

# aggregate
agg <- aggregate_files_patch(fcs_files = files,
                             output_file = "Aggregate.fcs",
                             output_dir = aggregate_dir,
                             write_agg_file = TRUE,
                             cTotal_ = 10000)



# -------------------------------------
# GATING
# -------------------------------------

# set gating output directory
gate_dir <- file.path(data_dir, "Gated")
if (!dir.exists(gate_dir)) dir.create(gate_dir)

# get files for gating
files <- list.files(path = aggregate_dir,
                    pattern = ".fcs$",
                    full.names = TRUE)

# gate files & plot gating strategy for each file
n_plots <- 3

png(file.path(gate_dir, paste0("gating.png")),
    width = n_plots * 300, 
    height = length(files) * 300)

layout(matrix(1:(length(files) * n_plots), 
              ncol = n_plots, 
              byrow = TRUE))

for (file in files){
  
  ff <- flowCore::read.FCS(filename = file, 
                           transformation = FALSE)
  
  ff <- gate_intact_cells(flow_frame = ff, 
                          file_name = basename(file))
  
  ff <- gate_singlet_cells(flow_frame = ff,
                           channels = "Event_length",
                           file_name = basename(file))
  
  ff <- gate_live_cells(flow_frame = ff, 
                        viability_channel = "Pt195Di",
                        out_dir = gate_dir)
  
  flowCore::write.FCS(ff, file.path(gate_dir,
                                    gsub(".fcs", "_gated.fcs", basename(file))))
}

dev.off()



# -------------------------------------
# NORMALIZATION USING REFERENCE SAMPLE (SKIPPED, NO BATCH)
# -------------------------------------


# -------------------------------------
# PLOT BATCH EFFECT (SKIPPED)
# -------------------------------------











