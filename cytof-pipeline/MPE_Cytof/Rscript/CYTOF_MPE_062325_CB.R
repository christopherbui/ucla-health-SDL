# ------------------------------------------------------------------------------
# Import packages
# ------------------------------------------------------------------------------

# Run in R v4.3.3

library(stringr)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(pals)
library(scales)
library(ggplot2)
library(pheatmap)
# library(pandoc)   # required by flowAI
library(cowplot)
library(matrixStats)

library(flowCore)
library(flowAI)
library(flowCut)
library(cytutils)   # clauclting aof/change fcs name

library(flowDensity)
library(CytoNorm)
library(CATALYST)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(diffcyt)
library(scater)
library(reshape2)
library(ComplexHeatmap)
library(FlowSOM)

library(progress)   # progress bar

source(paste0("U:/cdbui/MPE_Cytof/Rscript/Christopher_Scripts/", "Rybakowska_cytof_function.R"))
source(paste0("U:/cdbui/MPE_Cytof/Rscript/Christopher_Scripts/", "Rybakowska_cytof_function_LT.R"))
source("C:/Users/cdbui/Documents/GitHub/ucla-health-SDL/cytof-pipeline/MPE_Cytof/Rscript/Cytof_MPE_Utils.R")

# set working directory
workFolder <- paste0("U:/cdbui/", "MPE_Cytof")
setwd(workFolder)


# get sample info of all samples & panels
fin_info <- file.path("Ranalysis","mpe_cytof_sampleInfo_022625.txt")
allSampleInfo <- read.delim(fin_info, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
allSampleInfo <- allSampleInfo[, -1]




# ------------------------------------------------------------------------------
# Check quantile of REFERENCE samples in different batches 
#-------------------------------------------------------------------------------


# selPanel <- c("TBNK")  #*******
# selPanel <- c("Myeloid")  #*******
selPanel <- c("Cytokines")  #*******


# set output directories
panel_output_dir <- file.path(workFolder, "Ranalysis", selPanel)
if(!dir.exists(panel_output_dir)) dir.create(panel_output_dir)


# get sample info
sampleInfo  <- dplyr::filter(allSampleInfo, panel_id == selPanel & patient_id == "Ref")

# get panel info
fin_panel <- paste(selPanel, "_markers_022625.txt", sep = "")
panel_info <- read.delim(file.path("Ranalysis", fin_panel), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# set input folder
bead_norm_dir <- file.path(workFolder, "CYTOF_data", "BeadNorm", selPanel)

# get file names
select_fcs_files <- sampleInfo$file_name    #***must be in the same folder

# define batch id and sample id for each file
batch4extractPatt <- "(?i).*(batch[0-9]*).*.FCS"
batch_pattern <- batch4extractPatt

##file_a <- file.path(bead_norm_dir,select_fcs_files)
##batch_pattern <- str_match(basename(select_fcs_files), batch4extractPatt)[,2]

quantiles_ref <- extract_marker_quantiles_4SDL(
  select_fcs_files = select_fcs_files,
  in_dir = bead_norm_dir,
  batch_pattern = batch4extractPatt,
  arcsine_transform = TRUE,
  markers_to_plot = NULL)

# NOTE: Did not run code below; quantiles_ref already in wide format

# Ftab_quantile_ref <- reshape(quantiles_ref, idvar = c("File","Marker","Batch","Sample"),
#                       timevar = "Quantile", direction = "wide")

Ftab_quantile_ref <- quantiles_ref

# save quantile info
quantiles_ref_file <- paste0(selPanel, "_quantiles_Ref.txt")
write.table(Ftab_quantile_ref, file.path(panel_output_dir, quantiles_ref_file), sep = "\t", quote = FALSE, row.names = FALSE)




# --------------------------------------------------------------
# Check panel expression for reference samples
# 2/28/25
# --------------------------------------------------------------

selPanel <- c("TBNK")  #*******
# selPanel <- c("Myeloid") #*******
# selPanel <- c("Cytokines")

# get sample info
sampleInfo  <- dplyr::filter(allSampleInfo, panel_id == selPanel & patient_id == "Ref")

# get panel info
fin_panel <- paste(selPanel, "_markers_022625.txt", sep = "")
panel_info <- read.delim(file.path("Ranalysis", fin_panel), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# set input folder
bead_norm_dir <- file.path(workFolder, "CYTOF_data", "BeadNorm", selPanel)

# GET ONLY REF FILES; JUST FOR TESTING
all_fcs_files <- list.files(bead_norm_dir,
                            pattern = "_Ref.*\\.fcs$",
                            full.names = TRUE)

# non-transformed data
nfiles <- length(all_fcs_files)
quantile_values <-  c(0.01, 0.25, 0.5, 0.75, 0.99)
Ftab <- NULL


pb <- reloadProgressBar(nfiles)   # progress bar
for (i in c(1:nfiles)){   
  f_raw <- flowCore::read.FCS(all_fcs_files[i],transformation=FALSE)
  orgName <- unlist(keyword(f_raw,"$FIL"))
  
  sel_chanels <- grep("Di",colnames(f_raw),value=TRUE)
  quantiles_res <- t(apply(exprs(f_raw)[,sel_chanels],2,
                           function(x) stats::quantile(x,quantile_values)))
  colnames(quantiles_res) <- paste("quantile_",quantile_values,sep="")
  
  nchannels <- nrow(quantiles_res)
  tmp1 <- data.frame(file = rep(basename(all_fcs_files[i]),nchannels),
                     orgFileName = rep(orgName,nchannels),
                     channel_names = rownames(quantiles_res),
                     quantiles_res)
  Ftab <- bind_rows(Ftab,tmp1)
  pb$tick()
}
rownames(Ftab) <- NULL

Batch <- sapply(Ftab[,1],
                function(x) unlist(strsplit(x,split="_"))[4])
Ftab_out <- Ftab %>% mutate(Batch=Batch)

fout <- paste(selPanel,"_beadNorm_nonTransform_qc.txt",sep="")
write.table(Ftab_out, file.path(panel_output_dir, fout),
            sep="\t", quote=F,row.names=FALSE)



# transformed data
nfiles <- length(all_fcs_files)
quantile_values <-  c(0.01, 0.25, 0.5, 0.75, 0.99)
Ftab <- NULL


pb <- reloadProgressBar(nfiles)   # progress bar
for (i in c(1:nfiles)){   
  f_raw <- flowCore::read.FCS(all_fcs_files[i],transformation=FALSE)
  ff <- flowCore::transform(f_raw, transformList(grep("Di", colnames(f_raw), value = TRUE),
                                                 arcsinhTransform(a = 0, b = 1/5, c = 0)))
  
  orgName <- unlist(keyword(ff,"$FIL"))
  
  sel_chanels <- grep("Di",colnames(ff),value=TRUE)
  quantiles_res <- t(apply(exprs(ff)[,sel_chanels],2,
                           function(x) stats::quantile(x,quantile_values)))
  colnames(quantiles_res) <- paste("quantile_",quantile_values,sep="")
  
  nchannels <- nrow(quantiles_res)
  tmp1 <- data.frame(file = rep(basename(all_fcs_files[i]),nchannels),
                     orgFileName = rep(orgName,nchannels),
                     channel_names = rownames(quantiles_res),
                     quantiles_res)
  Ftab <- bind_rows(Ftab,tmp1)
  pb$tick()
}
rownames(Ftab) <- NULL

Batch <- sapply(Ftab[,1],
                function(x) unlist(strsplit(x,split="_"))[4])
Ftab_out <- Ftab %>% mutate(Batch=Batch)

fout <- paste(selPanel,"_beadNorm_arcsinhTransform_qc.txt",sep="")
write.table(Ftab_out,file.path(panel_output_dir, fout),
            sep="\t",quote=F,row.names=FALSE)




# -----------------------------------------
# Export events number
# -------------------------------------------

selPanel <- c("TBNK")  #*******
# selPanel <- c("Myeloid") #******
# selPanel <- c("Cytokines")


# set input folder
bead_norm_dir <- file.path(workFolder, "CYTOF_data", "BeadNorm", selPanel)

all_fcs_files <- list.files(bead_norm_dir,
                            pattern = ".fcs$",
                            full.names = TRUE)

# normalized, uncleaned
nfiles <- length(all_fcs_files)
Ftab <- matrix(NA, nfiles, 3)

# normalized, cleaned
all_fcs_files_cleaned <- list.files(clean_dir,
                                    pattern = ".fcs$",
                                    full.names = TRUE)

nfiles_cleaned <- length(all_fcs_files_cleaned)
Ftab <- matrix(NA,nfiles,3)


# change read.FCS(file_path) to 'all_fcs_files' OR 'all_fcs_files_cleaned' below
pb <- reloadProgressBar(nfiles_cleaned)
for (i in c(1:nfiles_cleaned)){   
  f_raw <- flowCore::read.FCS(all_fcs_files_cleaned[i],transformation=FALSE)
  orgName <- unlist(keyword(f_raw ,"$FIL"))
  Ftab[i,1] <- basename(all_fcs_files_cleaned[i])
  Ftab[i,2] <- orgName 
  Ftab[i,3] <- nrow(f_raw)
  pb$tick()
}
colnames(Ftab) <- c("file","FIL","event__no")

Ftab_file_name <- paste0(selPanel, "_events_number_cleaned.txt")
write.table(Ftab,
            file.path(panel_output_dir, Ftab_file_name),
            sep="\t",
            quote=F,
            row.names=FALSE)

# ------------------------------------------------------------

library(ggcyto)

fcs_files <- grep("Ref" , all_fcs_files,value = TRUE)
i <- 1
f <- flowCore::read.FCS(fcs_files[i])
keyword(f, "$FIL")

spillover(f)

tmp_channels <- colnames(f)
channel_info <- pData(parameters(f))

sel_idx <- which(grepl("CD45", channel_info$desc) == TRUE)
sel_channel <- channel_info$name[sel_idx]

ggcyto::autoplot(f, "Pt198Di")




# ------------------------------------------------------------------------------
# Signal Cleaning
#-------------------------------------------------------------------------------

list_panels <- c("TBNK")
# list_panels <- c("Myeloid")
# list_panels <- c("Cytokines")

seg_threshold <- 500

for (selPanel in list_panels){
  
  # get sample info
  # GET ONLY REF FILES; JUST FOR TESTING
  sampleInfo  <- dplyr::filter(allSampleInfo,
                               panel_id==selPanel) # & patient_id =="Ref")
  
  # get panel info
  fin_panel <- paste(selPanel,"_markers_022625.txt",sep="")
  panel_info <- read.delim(file.path("Ranalysis",fin_panel),
                           sep="\t",header=T,stringsAsFactors=FALSE)
  
  # set input directory (pathway to the files that are going to be cleaned)
  bead_norm_dir <- file.path(workFolder,"CYTOF_data",
                             "BeadNorm",selPanel)
  
  
  # set and create the directory where cleaned fcs files will be saved
  clean_dir <- file.path(workFolder,"CYTOF_data", 
                         "Cleaned",selPanel)
  if(!dir.exists(clean_dir)) dir.create(clean_dir, recursive = TRUE)
  
  # flowAI: set and create the directory where cleaned fcs files will be saved
  # flowAI_clean_dir <- file.path(getwd(), "CYTOF_data", "Cleaned-flowAI")
  # if (!dir.exists(flowAI_clean_dir)) dir.create(flowAI_clean_dir, recursive = TRUE)
  
  # Define which files will be cleaned
  select_fcs_files <- sampleInfo$file_name    #***must be in the same folder
  files_b <- file.path(bead_norm_dir,select_fcs_files)   #**full_path
  
  # select channel to clean
  # ****can be NULL if want to clean all of them
  dna_ch <- c("Ir191Di","Ir193Di")   #*****
  via_ch <- c("Pt195Di")    #***Cisplatin
  target_ch <- panel_info$fcs_colname
  channels_to_clean <- c(dna_ch,via_ch,target_ch)
  channels_to_exclude <- c("Cd108Di","Cs133Di","La139Di","Ce140Di","Ce142Di",
                           "Gd154Di","Gd157Di","Lu176Di")
  
  # Clean file by file in the loop, saving new file with each loop
  pb <- reloadProgressBar(length(files_b))
  for (file in files_b) {
    
    # read fcs file
    ff <- flowCore::read.FCS(filename = file,
                             transformation = FALSE)
    
    # Note: flowAI doesn't work
    # clean flow rate: use flowAI using Time
    # using single function flow_auto_qc
    # ff <- flowAI::flow_auto_qc(file,
    #            remove_from = "all",   #***FR,FS, and FM
    #           output = 1,            #**only clean signal==1
    #           timeCh = "Time",
    #           ChExcludeFS = channels_to_exclude,
    #           html_report=FALSE,
    #           mini_report=FALSE,
    #           fcs_QC = FALSE,
    #           folder_results = FALSE)
    # ff_2@exprs[,"QCvector"]-> ttmp
    
    ff_3 <- clean_signal(flow_frame = ff,
                         channels_to_clean=channels_to_clean,
                         to_plot = "None",
                         out_dir = clean_dir,
                         Segment = seg_threshold,
                         arcsine_transform = TRUE,
                         data_type = "MC",
                         non_used_bead_ch = "140")
    
    # Write FCS files
    flowCore::write.FCS(ff_3,
                        file = file.path(clean_dir, gsub("_beadNorm","_cleaned", basename(file))))
    
    # Write FCS files from flowAI
    # flowCore::write.FCS(ff,
    #                     file = file.path(flowAI_clean_dir, gsub("_beadNorm", "_flowAIcleaned", basename(file))))
    
    pb$tick()
  }
}




# ------------------------------------------------------------------------------
# Outlier Detection
# 
# NOTE: output is the score from fsom
#-------------------------------------------------------------------------------

list_panels <- c("TBNK")
# list_panels <- c("Myeloid")
# list_panels <- c("Cytokines")

batch4extractPatt <- "(?i).*(batch[0-9]*).*.fcs"
dna_ch_4fsom <- c("191Ir","193Ir")   #*****using desc column
via_ch_4fsom <- c("195Pt")    #***Cisplatin

nCells_thres <- c(500)


for (k in c(1:length(list_panels))){
  
  selPanel <- list_panels[k]
  tmp_nCells_thres <- nCells_thres[1]
  
  # get sample info
  sampleInfo  <- dplyr::filter(allSampleInfo,
                               panel_id==selPanel)  #**** & patient_id =="Ref")
  
  # get panel info --> function and lineage markers
  fin_panel <- paste(selPanel,"_markers_022625.txt",sep="")
  panel_info <- read.delim(file.path("Ranalysis",fin_panel),
                           sep="\t",header=T,stringsAsFactors=FALSE)
  
  lineage_idx <- which(panel_info$marker_class=="lineage")  #****
  lineage_markers <- panel_info$fcs_desc[lineage_idx]  
  function_idx <- which(panel_info$marker_class=="function")  #****
  function_markers <- panel_info$fcs_desc[function_idx] 
  pheno_4fsom <- c(dna_ch_4fsom,via_ch_4fsom,lineage_markers)    #*******
  
  # set input directory where cleaned fcs files will be saved
  clean_dir <- file.path(workFolder,"CYTOF_data", 
                         "Cleaned", selPanel)
  
  # Define out_dir for diagnostic plots
  quality_dir <- file.path(workFolder,"CYTOF_data",
                           "Quality_control", selPanel)
  if(!dir.exists(quality_dir)) dir.create(quality_dir, recursive = TRUE)
  # quality_dir <- file.path("C:/Users/cdbui/Documents", "Quality_control_DEBUG_PlotRedDim", selPanel)
  # if(!dir.exists(quality_dir)) dir.create(quality_dir, recursive = TRUE)
  
  # Define which files will be cleaned
  ###select_fcs_files <- sampleInfo$file_name    #***must be in the same folder
  files_b <- list.files(clean_dir, 
                        pattern = "_cleaned.fcs$", 
                        full.names = TRUE)
  
  # Define batch_id for each file 
  file_batch_id <- stringr::str_match(basename(files_b),batch4extractPatt)[,2] #*****
  
  ###nCells = length(fcs_files)*tmp_nCells_thres
  file_quality_check_4SDL(fcs_files = files_b,
                          file_batch_id = file_batch_id,
                          nCellsPerSample=tmp_nCells_thres,
                          out_dir = quality_dir,
                          phenotyping_markers = pheno_4fsom,
                          arcsine_transform = TRUE,
                          nClus = 10,
                          sd = 3)
}




# ------------------------------------------------------------------------------
# Gating

# Files gating using new strategies:--> *_4SDL functions
# for intact: single low threshold for Ir191 and Ir193 (i.e. >thres)
# for live: only use viability channel(cisplatin - i.e. < thres)
# ****Three gates: intact (Ir193 and Ir191), doublet (Time Event) and live (Cisplatin)
# --> renames to Gated_2gated
#-------------------------------------------------------------------------------

selPanel <- c("TBNK")  #*******
# selPanel <- c("Myeloid") #******
# selPanel <- c("Cytokines")

# singlet deltaT threshold
event_length_cutoff <- 50

# intact thresholds
intact_accepted_thres <- c(1.5, 4)
intact_adjusted_thres <- c(2, 3)

# viability thresholds
live_accepted_thres <- c(3, 6)
live_adjusted_thres <- c(4, 6)


# Set input directory - clean dir
aggregate_dir <- file.path(workFolder, "CYTOF_data", "Cleaned", selPanel)

# Set output directory
gate_dir <- file.path(workFolder, "CYTOF_data", "Gated", selPanel)
if (!dir.exists(gate_dir)) dir.create(gate_dir, recursive = TRUE)

# List files for gating
files <- list.files(path = aggregate_dir,
                    pattern = ".fcs$",
                    full.names = TRUE)

# Gate the files and plot gating strategy for each file
n_plots <- 3 # intact, singlet, live
png(file.path(gate_dir, paste0("gating.png")),
    width = n_plots * 300,
    height = length(files) * 300)

layout(matrix(1:(length(files) * n_plots),
              ncol = n_plots, 
              byrow = TRUE))

message("Started: ", format(Sys.time(), tz = "America/Los_Angeles"))
pb <- reloadProgressBar(length(files))

# master df to hold all gating information
master_df <- NULL

for (file in files) {
  message("\n", (basename(file)))
  
  ff <- flowCore::read.FCS(filename = file, 
                           transformation = FALSE)
  
  res_singlet <- gate_singlet_cells_4SDL(flow_frame = ff,
                                         channels = "Event_length",
                                         file_name = basename(file),
                                         event_length_cutoff = event_length_cutoff)
  
  # function internally considers "Ir193Di", "Ir191Di"
  res_intact <- gate_intact_cells_4SDL(flow_frame = res_singlet$flowFrame,
                                       file_name = basename(file),
                                       accepted_thres = intact_accepted_thres,
                                       adjusted_thres = intact_adjusted_thres)
  
  res_live <- gate_live_cells_4SDL(flow_frame = res_intact$flowFrame,
                                   viability_channel = "Pt195Di", # cisplatin; higher binding ouput with dead cells
                                   out_dir = gate_dir,
                                   accepted_thres = live_accepted_thres,
                                   adjusted_thres = live_adjusted_thres)
  
  
  # gating thresholds
  df_intact <- res_intact$info_df
  df_singlet <- res_singlet$info_df
  df_live <- res_live$info_df
  
  # join gating threshold information
  tmp_merged_df <- merge(df_intact, df_live, by = "file_name", all = TRUE)
  merged_df <- merge(tmp_merged_df, df_singlet, by = "file_name", all = TRUE)
  
  # append gating threshold information to master df
  master_df <- bind_rows(master_df, merged_df)
  
  # output gated fcs file
  flowCore::write.FCS(res_live$flowFrame, file.path(gate_dir, gsub(".fcs", "_gated.fcs", basename(file))))
  
  pb$tick()
}

# save gating threshold information for all fcs files
write.csv(master_df,
          file = file.path(gate_dir, paste0(selPanel, "_Gating_Thresholds.csv")),
          row.names = FALSE)

dev.off()
message("Ended: ", format(Sys.time(), tz = "America/Los_Angeles"))




# ------------------------------------------------------------------------------
# Normalization Using Reference Samples
#
# NOTE: 6/17/2025
#   - SKIPPED normalization due to some reference samples causing bugs during outlier step.
#   - Affected Panels: TBNK, Cytokine
#   - All Myeloid samples passed through outlier step ok.
#-------------------------------------------------------------------------------

# select panel
selPanel <- c("TBNK")
# selPanel <- c("Myeloid")
# selPanel <- c("Cytokines")

# set parameters
sel4ref <- c(1:7)     #***
batch4extractPatt <- "(?i).*(batch[0-9]*).*.fcs"  #****
refPattern <- c("Ref.*_gated.fcs$")

# get channel info
dna_ch <- c("Ir191Di","Ir193Di")   #*****
via_ch <- c("Pt195Di")    #***Cisplatin
target_ch <- panel_info$fcs_colname

channels_to_norm <- c(dna_ch, via_ch, target_ch)

# get panel info
fin_panel <- paste(selPanel, "_markers_022625.txt", sep="")
panel_info <- read.delim(file.path("Ranalysis", fin_panel), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# set input directory
gate_dir <- file.path(workFolder, "CYTOF_data", "Gated", selPanel)

# define reference samples
files_ref <- list.files(gate_dir,
                        pattern = refPattern,
                        full.names = TRUE,
                        recursive = TRUE)

files_ref <- files_ref[sel4ref]   #****only for TBNK

# definte batch labels for each file
labels_ref <- stringr::str_match(basename(files_ref), batch4extractPatt)[, 2]

# define channels to be normalized by reading in fcs to & get channel names
ff <- read.FCS(files_ref[1])
channels <- colnames(ff)[colnames(ff) %in% channels_to_norm]

# set output directory
norm_dir <- file.path(workFolder, "CYTOF_data", "CytoNormed", selPanel)
if (!dir.exists(norm_dir)) (dir.create(norm_dir, recursive = TRUE))

# build normalization model using reference samples & plot quantiles
png(file.path(norm_dir, "005_095_normalization.png"),
    width = length(channels) * 300,
    height = (length(files_ref) * 2 + 1) * 300)

model <- CytoNorm::QuantileNorm.train(files = files_ref,
                                      labels = labels_ref,
                                      channels = channels,
                                      transformList = transformList(channels, CytoNorm::cytofTransform),
                                      nQ = 2,
                                      limit = c(0, 8),
                                      quantileValues = c(0.05, 0.95),
                                      goal = "mean",
                                      plot = TRUE)
dev.off()

# save model
saveRDS(object = model,
        file = file.path(norm_dir, "005_095_model.RDS"))

# set path to files to be normalized
files <- list.files(file.path(gate_dir),
                    pattern = "_gated.fcs$",
                    full.names = TRUE,
                    recursive = TRUE)

# define batch labels for each file
# file's batch label corresponds to reference labels
labels <- stringr::str_match(basename(files), batch4extractPatt)[, 2]

# normalize files
CytoNorm::QuantileNorm.normalize(model = model,
                                 files = files,
                                 labels = labels,
                                 transformList = transformList(channels, CytoNorm::cytofTransform),
                                 transformList.reverse = transformList(channels, CytoNorm::cytofTransform.reverse),
                                 outputDir = norm_dir,
                                 prefix = "")




# ------------------------------------------------------------------------------
# Change $FIL for normalized, gated files
#
# NOTE: 6/17/2025 
#   - Changed $FIL for NON-normalized, gated files to keep workflow consistent
#     while skipping normalization.
#-------------------------------------------------------------------------------

selPanel <- c("TBNK")
# selPanel <- c("Myeloid")
# selPanel <- c("Cytokines")

# input_dir <-  file.path(workFolder, "CYTOF_data", "Gated", selPanel)
input_dir <- file.path(workFolder, "CYTOF_data", "CytoNormed", selPanel)

fcs_files <- list.files(input_dir, 
                        pattern = ".fcs$", 
                        full.names = TRUE)

change_fcs_FIL(fcs_files)




# ------------------------------------------------------------------------------
# Run FSOM & UMAP by CATALYST 
# 
# NOTE:
#   - Using non-normalized, gated data
#-------------------------------------------------------------------------------

# get sample info
fin_info <- file.path("Ranalysis","mpe_cytof_sampleInfo_022625.txt")
allSampleInfo <- read.delim(fin_info, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
allSampleInfo <- allSampleInfo[, -1]

# select panel
selPanel <- c("TBNK")
# selPanel <- c("Myeloid")
# selPanel <- c("Cytokines")

# get panel info
fin_panel <- paste(selPanel, "_markers_022625.txt", sep="")
panel_info <- read.delim(file.path("Ranalysis", fin_panel), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# get channel info for fsom
lineage_idx <- which(panel_info$marker_class == "lineage")
lineage_markers <- panel_info$fcs_desc[lineage_idx]
function_idx <- which(panel_info$marker_class == "function")
function_markers <- panel_info$fcs_desc[function_idx]

# set panel input for CATALYST, requires 3 columns
tmp_panel <- panel_info[, c("fcs_colname", "antigen", "marker_class")]
tmpi <- which(grepl("function", panel_info$marker_class) == TRUE)
tmp_panel$marker_class[tmpi] <- c("state")
tmpi <- which(panel_info$marker_class == "lineage")
tmp_panel$marker_class[tmpi] <- c("type")


# set input directory
norm_dir <- file.path(workFolder, "CYTOF_data", "CytoNormed", selPanel, "updateFileName")

# set parent output directory
analysis_dir <- file.path(workFolder, "CYTOF_data", "Analysis", selPanel)
if (!dir.exists(analysis_dir)) dir.create(analysis_dir, recursive = TRUE)

# set results directory
res_dir <- file.path(analysis_dir, "Results_Full_Lineage")
if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)

# get fcs files
fcs_files <- list.files(norm_dir,
                        pattern = ".fcs$",
                        full.names = TRUE)

# get fcs corename
tmp_core <- sapply(basename(fcs_files), getCoreID, simplify = TRUE)

# match filenames and metadata
# prepare metadata for CATALYST
sampleInfo <- dplyr::filter(allSampleInfo, panel_id == selPanel)
md_info <- left_join(data.frame(file = basename(fcs_files), sample_id = tmp_core),
                     sampleInfo[, -1],
                     by = join_by(sample_id == corename))

md_info$BATCH <- factor(md_info$BATCH)

# import fcs as sce object (arcsinh, transform = TRUE by default)
sce <- CATALYST::prepData(fcs_files,
                          panel = tmp_panel, # panel marker info
                          md = md_info, # sample info
                          panel_cols = list(channel = "fcs_colname", antigen = "antigen", class = "marker_class"),
                          md_cols = list(file = "file", id = "sample_id", factors = c("tissue_type", "patient_id", "parental_sample_id", "BATCH")))

SummarizedExperiment(sce)

colData(sce)$sample_id <- droplevels(colData(sce)$sample_id)

# add "Symbol" marker names for plotting; preferred for publishing
matched_symbols <- panel_info$Symbol[match(rownames(sce), panel_info$antigen)]
rowData(sce)$Symbol <- matched_symbols


# SAVE MAIN SCE OBJECT
# rds_path <- "C:/Users/cdbui/Desktop/SCE_RDS"
# file_name <- "sce_main.rds"
# saveRDS(sce, file = file.path(rds_path, file_name))


# GENERAL SCE OBJECT INFO
#-------------------------------------------------------------------------------
# check distribution of markers
tmp2 <- t(apply(assay(sce, "exprs"), 1, function(x) quantile(x, c(0.5, 0.75, 0.9, 0.95, 1))))
tmp_mean <- rowMeans(assay(sce, "exprs"))
tmpout <- bind_cols(tmp_mean, tmp2)
colnames(tmpout)[1] <- c("mean")
tmpout$marker <- rowData(sce)$Symbol

write.table(tmpout, "tmp.txt", sep="\t", quote=F,row.names=F)

#-----------------------------------

# number of events per file
tmp <- table(colData(sce)$sample_id)    # number of events per file
tmp_2 <- data.frame(tmp)
colnames(tmp_2)[1] <- c("sample_id")

tmp_2 <- left_join(tmp_2, md_info, by = join_by(sample_id))
kruskal.test(tmp_2$Freq~factor(tmp_2$BATCH))
# Kruskal-Wallis rank sum test
# 
# data:  tmp_2$Freq by factor(tmp_2$BATCH)
# Kruskal-Wallis chi-squared = 5.2016, df = 6, p-value = 0.5182
tmp_2$BATCH <- factor(tmp_2$BATCH)

tissue_freq_by_batch_PLOT <- ggplot(tmp_2, aes(x = BATCH, y = Freq)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = tissue_type), width = 0.25) +
  scale_y_continuous(trans = "log10") + theme_bw()

file_name <- paste0(selPanel, "_tissue_freq_by_batch", ".png")
ggsave(filename = file.path(analysis_dir, file_name), plot = tissue_freq_by_batch_PLOT)

write.table(tmp, file.path(analysis_dir, "num_events_per_sample.txt"), sep = "\t", quote = FALSE, col.names = FALSE)

#-----------------------------------

# plot distribution of markers
color_by <- "BATCH"
marker_dist_PLOT <- plotExprs(sce, color_by = color_by)

file_name <- paste0("marker_distribution_by_", color_by, ".png")
ggsave(filename = file.path(analysis_dir, file_name), plot = marker_dist_PLOT, width = 10, height = 7)

#-----------------------------------

# plot mulit-dimension scale (median marker intensities)
color_by <- "BATCH"
label_by <- "sample_id"
pb_mds_PLOT <- pbMDS(sce, color_by = color_by, label_by = label_by)

file_name <- paste0("Ref_pbMDS_by_", color_by, "_", label_by, ".png")
ggsave(filename = file.path(analysis_dir, file_name), plot = pb_mds_PLOT, width = 8, height = 6)

#-----------------------------------

# heatmap of median marker intensities
# in CyTOF applications, a cosine distance shows good performance
# scale = "last": aggregate then scale & trim
marker_heatmap_PLOT <- plotExprHeatmap(sce,
                                       scale = "last",
                                       row_anno = FALSE,
                                       col_anno = FALSE)

file_name <- paste0("expr_heatmap_all_markers", ".png")
png(file.path(analysis_dir, file_name), width = 14, height = 8, units = "in", res = 300)
draw(marker_heatmap_PLOT)
dev.off()


# ------------------------------------------------------------------------------

# Select Markers To Use

# FULL lineage markers
type_markers <- rownames(sce)[rowData(sce)$marker_class == "type"]

# SUBSET lineage markers
# NOTE:
#   - Get associated `antigen` marker_name from rowData() for subsetting sce object 
lin_subset_markers <- readLines(file.path(analysis_dir, "lineage_subset_markers_antigen.txt"))

# select markers

# sel_markers <- type_markers
sel_markers <- lin_subset_markers

# partition sce object for selected markers
# sce_full <- sce[sel_markers, ]
sce_sub <- sce[sel_markers, ]


# SAVE SCE SUBSET
# rds_path <- "C:/Users/cdbui/Desktop/SCE_RDS"
# file_name <- "sce_subset_lineage.rds"
# saveRDS(sce_sub, file = file.path(rds_path, file_name))


# set cluster info directory
cluster_info_dir <- file.path(res_dir, "Cluster_Info")
if (!dir.exists(cluster_info_dir)) dir.create(cluster_info_dir, recursive = TRUE)



# FSOM -------------------------------------------------------------------------
maxk <- 20

# sce_tmp <- sce_full
sce_tmp <- sce_sub

sce_tmp <- cluster(sce_tmp,
                   features = sel_markers,
                   xdim = 10,
                   ydim = 10,
                   maxK = maxk,
                   seed = 1234)
# saveRDS
# rds_path <- "C:/Users/cdbui/Desktop/SCE_RDS"
# file_name <- "sce_subset_lineage_AFTER_FSOM.rds"
# saveRDS(sce_sub, file = file.path(rds_path, file_name))




# matrix showing relationship between SOM (n=100) and consensus cluster (metaX)
tmp_cluster_codes <- cluster_codes(sce_tmp)

# add meta clusters to column data
sce_tmp$meta20 <- cluster_ids(sce_tmp, "meta20")

# sce_full_lin$meta20 <- cluster_ids(sce_full_lin, "meta20")


# export expression of cluster data
# sel_metaK <- c("meta20")   # consensus cluster level or SOM
# sel_aggregate <- c("median")
# 
# raw_clust_expr <- CATALYST:::.agg(sce_sub, by = sel_metaK, assay = "exprs", fun = sel_aggregate)
# raw_clust_expr_df <- as.data.frame(raw_clust_expr)
# raw_clust_expr_df <- tibble::rownames_to_column(raw_clust_expr_df, "marker")
# 
# file_name <- paste0(sel_metaK, "_", sel_aggregate , "_expr_matrix", ".txt")
# write.table(raw_clust_expr_df, file.path(cluster_info_dir, file_name), sep = "\t", quote = FALSE, row.names=FALSE)




# PCA --------------------------------------------------------------------------
pc_prop <- 0.75 #*** proportion of total features to use
no_pcs <- floor(nrow(sce_tmp) * pc_prop)

sce_tmp <- runPCA(sce_tmp, exprs_values = "exprs", ncomponents = floor(nrow(sce_tmp) * pc_prop))
# sce_full_lin <- runPCA(sce_full_lin, exprs_values = "exprs", ncomponents = floor(nrow(sce_full_lin) * pc_prop))


reducedDimNames(sce_tmp)  #"PCA"
dim(reducedDim(sce_tmp,"PCA"))   #ncells x no_pc


# elbow plot
pca_matrix <- reducedDim(sce_tmp, "PCA")

# either for non-normalized (sum<1)
var_explained <- attr(pca_matrix, "varExplained")  #***non-normalize to 1
var_percent <- attr(pca_matrix, "percentVar")   #**in percent
pca_rotation <- attr(pca_matrix, "rotation") #***markers x no_pc

plot(var_percent)

elbow_PLOT <- plot_elbow_plot(var_percent)
file_name <- paste0("PCA_Elbow_Plot", ".png")
ggsave(filename = file.path(res_dir, file_name), plot = elbow_PLOT, height = 8, width = 12)



# UMAP -------------------------------------------------------------------------
pc_to_use <- 10

# parent UMAP directory
umap_dir <- file.path(res_dir, "UMAP")
if (!dir.exists(umap_dir)) dir.create(umap_dir, recursive = TRUE)

# directory specifying PC used for UMAP
umap_pc_dir <- file.path(umap_dir, paste0("PC_", pc_to_use))
if (!dir.exists(umap_pc_dir)) dir.create(umap_pc_dir, recursive = TRUE)

# remove cells = 1e3; use ~10% of avg_cells_per_sample
#avg_cells_per_sample <- mean(table(colData(sce_sub)$sample_id))
n_cells <- 2e4
sce_tmp <- runDR(sce_tmp, "UMAP", cells = n_cells, features = sel_markers, pca = pc_to_use, seed = 1234)

# specify cluster for plotting
sel_metaK <- c("meta20")

umap_PLOT <- plotDR(sce_sub, "UMAP", color_by = cluster_name)
file_name <- paste0("UMAP_", cluster_name, ".png")
ggsave(filename = file.path(umap_pc_dir, file_name), plot = umap_PLOT, width = 10, height = 10, bg = "white")

umap_facet_PLOT <- plotDR(sce_sub, "UMAP", color_by = cluster_name, facet_by = cluster_name)
file_name <- paste0("UMAP_", cluster_name, "_facet_by_cluster", ".png")
ggsave(filename = file.path(umap_pc_dir, file_name), plot = umap_facet_PLOT, width = 12, height = 10, bg = "white")



# Dotplot ----------------------------------------------------------------------
dotplot_dir <- file.path(res_dir, "Dotplot")
if (!dir.exists(dotplot_dir)) dir.create(dotplot_dir, recursive = TRUE)

# no scaling
dotplot_no_scaled_dir <- file.path(dotplot_dir, "Non_Scaled")
if (!dir.exists(dotplot_no_scaled_dir)) dir.create(dotplot_no_scaled_dir, recursive = TRUE)

# scaling
dotplot_scaled_dir <- file.path(dotplot_dir, "Scaled")
if (!dir.exists(dotplot_scaled_dir)) dir.create(dotplot_scaled_dir, recursive = TRUE)


sel_metaK <- c("meta20")
sel_aggregate <- c("median")
scale_option <- TRUE
label_for_dotplot <- paste(ifelse(scale_option, "scaled", "non-scaled"), sel_aggregate, "exprs", sep = " ")


# set output directory; i.e. scaling or no scaling
dp_output_dir <- ifelse(scale_option, dotplot_scaled_dir, dotplot_no_scaled_dir)


tmp_ftab <- dotplotTables(sce_tmp,
                          cluster_name = sel_metaK,
                          assay = "exprs",
                          fun = sel_aggregate,
                          scale = scale_option,
                          q = 0.01)

tmp_fig <- dotplotFig(tmp_ftab$expr_long, lab = label_for_dotplot)


# save dotplot matrices
prefix <- paste0(sel_metaK, "_dotplot_", sel_aggregate, ifelse(scale_option, "_SCALED", "_NON_SCALED"))
write.table(tmp_ftab$expr_wide, file = file.path(dp_output_dir, paste0(prefix, "_wide.txt")),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(tmp_ftab$expr_long, file = file.path(dp_output_dir, paste0(prefix, "_long.txt")),
            sep = "\t", quote = FALSE, col.names = NA)
  

# save dotplot figure
prefix <- paste0(sel_metaK, "_dotplot_", sel_aggregate, ifelse(scale_option, "_SCALED", "_NON_SCALED"))
ggsave(filename = file.path(dp_output_dir, paste0(prefix, ".png")), plot = tmp_fig, width = 10, height = 8)



# bluster silouette ------------------------------------------------------------

# directory of RDS files
rds_path <- file.path(analysis_dir, "RDS")

# load sce object
sce_tmp <- readRDS(file.path(rds_path, "sce_subset_lineage.rds"))


# list of cluster codes
metaK_id  <- colnames(sce_tmp@metadata$cluster_codes)[5:20]

# metacluster id for each cell
tmp_meta_clust <- lapply(as.list(metaK_id), function(id) CATALYST::cluster_ids(sce_tmp, id))

names(tmp_meta_clust) <- metaK_id

mat_expr <- t(as.matrix(assay(sce_tmp,  "exprs")))   #***DO NOT USE assay="exprs" since it gets raw data

sil_expr <- vapply(tmp_meta_clust, function(x) mean(bluster::approxSilhouette(mat_expr, x)$width), 0)

nclusters <- as.numeric(sapply(metaK_id, function(x) substr(x, 5, nchar(x))))

plot(nclusters, sil_expr, xlab="Number of clusters", ylab="Average silhouette width")


# create dataframe required for ggplot
df_sil <- data.frame(
  n_clusters = as.numeric(sub("meta", "", metaK_id)),
  silhouette = sil_expr
)
sil_PLOT <- ggplot(df_sil, aes(x = n_clusters, y = silhouette)) +
  geom_point(color = "#4287f5", size = 2) +
  geom_line(color = "#4287f5", linewidth = 0.9) +
  geom_text(aes(label = round(silhouette, 3)),
            hjust = -0.3, vjust = -0.3, size = 3.5) +
  scale_x_continuous(
    breaks = df_sil$n_clusters,
    labels = df_sil$n_clusters
  ) +
  labs(
    x = "Number of Clusters",
    y = "Average Silhouette Width",
    title = "Silhouette Width by Meta-Cluster Count - Cytokine Subset Lineage"
  ) +
  theme_light(base_size = 14)

ggsave(filename = file.path(res_dir, "Silhouette_Plot.png"), plot = sil_PLOT, width = 12, height = 8)



# ------------------------------------------------------------------------------
# DIFFERENTIAL ANALYSIS
#
#-------------------------------------------------------------------------------

# directory of RDS files
rds_path <- file.path(analysis_dir, "RDS")

# load sce object
# rds_to_use <- "sce_main.rds"
rds_to_use <- "sce_subset_lineage.rds"

sce_tmp <- readRDS(file.path(rds_path, rds_to_use))

sce_full_lineage <- readRDS(file.path(rds_path, "sce_full_lineage.rds"))
sce_sub_lineage <- readRDS(file.path(rds_path, "sce_subset_lineage.rds"))

# add meta6 to colData
# sce_tmp$meta6 <- cluster_ids(sce_tmp, "meta6")  # works if sce_tmp has metadata from clustering

# add cluster ids from sce_full_lineage to sce_main
cluster_name <- c("meta6")  # check silhouette plot & confirm cluster for each panel

colData(sce_full_lineage)[[cluster_name]] <- cluster_ids(sce_full_lineage, cluster_name)

sce_tmp$cluster_id <- sce_full_lineage$cluster_id
colData(sce_tmp)[[cluster_name]] <- colData(sce_full_lineage)[[cluster_name]]

# copy metadata; need for diffcyt
metadata(sce_tmp)$cluster_codes <- metadata(sce_full_lineage)$cluster_codes


sce_tmp_0 <- sce_tmp
# IMPORTANT:
# - Filter out Ref PBMC samples
sce_tmp <- filterSCE(sce_tmp, patient_id != "Ref")

# NOTE:
# - By default, MPE is the reference; Relevel the tissue type
# relevel: set PBMC as reference
old_level <- levels(colData(sce_tmp)$tissue_type)
colData(sce_tmp)$tissue_type <- factor(colData(sce_tmp)$tissue_type, levels = c("PBMC", "MPE"))
new_level <- levels(colData(sce_tmp)$tissue_type)


# diffcyt ----------------------------------------------------------------------

condition <- "tissue_type"
design <- createDesignMatrix(ei(sce_tmp), cols_design = condition)
contrast <- createContrast(c(0, 1))


# differential abundance (DA) of clusters
cluster_name <- c("meta5")
res_DA <- diffcyt(sce_tmp,
                  clustering_to_use = cluster_name,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  design = design,
                  contrast = contrast,
                  verbose = FALSE)

# differential states (DS) within clusters
# "DS" expects functional markers by default; need to specify lineage markers
type_markers <- rowData(sce_tmp)$marker_class == "type"

sel_markers <- rownames(sce_tmp) %in% rownames(sce_tmp) # all markers
# sel_markers <- type_markers


cluster_name <- c("meta5")
res_DS <- diffcyt(sce_tmp,
                  clustering_to_use = cluster_name,
                  analysis_type = "DS",
                  method_DS = "diffcyt-DS-limma",
                  design = design,
                  contrast = contrast,
                  markers_to_test = sel_markers,
                  verbose = FALSE)


rds_path <- paste0("C:/Users/cdbui/Desktop/SCE_RDS_COMP/", selPanel)

# save RDS
file_name <- paste0("res_DA_", cluster_name, "_", condition, ".rds")
saveRDS(res_DA, file = file.path(rds_path, file_name))

file_name <- paste0("res_DS_", cluster_name, "_", condition, ".rds")
saveRDS(res_DS, file = file.path(rds_path, file_name))



# scran - group: clusters ------------------------------------------------------
cluster_name <- c("meta5")
tmp_de_pv <- scran::findMarkers(sce_tmp,
                                groups=cluster_ids(sce_tmp, cluster_name),
                                assay.type = "exprs",
                                test.type = "wilcox",
                                pval.type = "all",   ##c("any", "some", "all"),
                                min.prop = 0.25,   #default
                                direction="up")
tmp_de_all <- bind_rows(
  lapply(names(tmp_de_pv), function(clust) {
    df <- as.data.frame(tmp_de_pv[[clust]])
    df$marker <- rownames(tmp_de_pv[[clust]])
    df$cluster_id <- clust
    df
  })
)
tmp_de_all <- tmp_de_all %>%
  relocate(marker, .before = everything()) %>%
  relocate(cluster_id, .after = marker) %>%
  relocate(AUC.1, .before = AUC.2)

rownames(tmp_de_all) <- NULL



sum_out <- scran::summaryMarkerStats(sce_tmp,
                                     groups = cluster_ids(sce_tmp, cluster_name),
                                     assay.type = "exprs",
                                     average = "median")

sum_out_all <- bind_rows(
  lapply(names(sum_out), function(clust) {
    df <- as.data.frame(sum_out[[clust]])
    df$marker <- rownames(sum_out[[clust]])
    df$cluster_id <- clust
    df
  })
)
sum_out_all <- sum_out_all %>%
  relocate(marker, .before = everything()) %>%
  relocate(cluster_id, .after = marker)

rownames(sum_out_all) <- NULL

# write tables
diff_expr_dir <- file.path(analysis_dir, "Diff_Expr")
if (!dir.exists(diff_expr_dir)) dir.create(diff_expr_dir, recursive = TRUE)

con <- c("Cluster")
diff_expr_condition_dir <- file.path(diff_expr_dir, con)
if (!dir.exists(diff_expr_condition_dir)) dir.create(diff_expr_condition_dir)

file_name <- paste0("scran_de_", cluster_name, "_all_markers.txt")
write.table(tmp_de_all, file = file.path(diff_expr_condition_dir, file_name), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

file_name <- paste0("scran_marker_sum_", cluster_name, "_all_markers.txt")
write.table(sum_out_all, file = file.path(diff_expr_condition_dir, file_name), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)




# scran - group: tissue type ---------------------------------------------------
tmp_de_pv <- scran::findMarkers(sce_tmp,
                                groups=colData(sce_tmp)$tissue_type,
                                assay.type = "exprs",
                                test.type = "wilcox",
                                pval.type = "all",   ##c("any", "some", "all"),
                                min.prop = 0.25,   #default
                                direction="up")
tmp_de_all <- bind_rows(
  lapply(names(tmp_de_pv), function(tissue_type) {
    df <- as.data.frame(tmp_de_pv[[tissue_type]])
    df$marker <- rownames(tmp_de_pv[[tissue_type]])
    df$tissue_type <- tissue_type
    df
  })
)
tmp_de_all <- tmp_de_all %>%
  relocate(marker, .before = everything()) %>%
  relocate(tissue_type, .after = marker) %>%
  relocate(AUC.PBMC, .before = AUC.MPE)

rownames(tmp_de_all) <- NULL




sum_out <- scran::summaryMarkerStats(sce_tmp,
                                     groups = colData(sce_tmp)$tissue_type,
                                     assay.type = "exprs",
                                     average = "median")
sum_out_all <- bind_rows(
  lapply(names(sum_out), function(tissue_type) {
    df <- as.data.frame(sum_out[[tissue_type]])
    df$marker <- rownames(sum_out[[tissue_type]])
    df$tissue_type <- tissue_type
    df
  })
)
sum_out_all <- sum_out_all %>%
  relocate(marker, .before = everything()) %>%
  relocate(tissue_type, .after = marker)

rownames(sum_out_all) <- NULL


# write tables
diff_expr_dir <- file.path(analysis_dir, "Diff_Expr")
if (!dir.exists(diff_expr_dir)) dir.create(diff_expr_dir, recursive = TRUE)

con <- c("Tissue_Type")
diff_expr_condition_dir <- file.path(diff_expr_dir, con)
if (!dir.exists(diff_expr_condition_dir)) dir.create(diff_expr_condition_dir)

file_name <- paste0("scran_de_", "tissue", "_all_markers.txt")
write.table(tmp_de_all, file = file.path(diff_expr_condition_dir, file_name), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

file_name <- paste0("scran_marker_sum_", "tissue", "_all_markers.txt")
write.table(sum_out_all, file = file.path(diff_expr_condition_dir, file_name), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)








# boxplots of cluster abundance by batch/tissue_type ---------------------------

cluster_name <- c("meta5")

# boxplot of cluster abundance within sample, subplot by batch, facet by cluster
abun1_PLOT <- CATALYST::plotAbundances(sce_tmp, k = cluster_name, by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
ggsave(file.path(analysis_dir, "abun_cluster_batch_boxplot.png"), plot = abun1_PLOT, width = 10, height = 6, dpi = 300)

# boxplot of tissue type abundance within sample, subplot by tissue type, facet by cluster
abun2_PLOT <- CATALYST::plotAbundances(sce_tmp, k = cluster_name, by = "cluster_id", group_by = "tissue_type") +
  ggtitle("Abundance by Cluster & Tissue Type")
ggsave(file.path(analysis_dir, "abun_cluster_tissue_type_boxplot.png"), plot = abun2_PLOT, width = 10, height = 6, dpi = 300)

# barplot with x = each patient, subplot by BATCH
CATALYST::plotAbundances(sce_tmp, k = cluster_name, by = "sample_id", group_by = "BATCH")



# use PBMC to check batch effect ------------------------
cluster_name <- c("meta6")

sce_all_pbmc <- filterSCE(sce_tmp_0, tissue_type == "PBMC")

clust_tbl_pbmc <- table(cluster_ids(sce_all_pbmc, cluster_name), colData(sce_all_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
png(file.path(analysis_dir, "pbmc_sample_cluster_heatmap.png"), width = 1200, height = 1000)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")
dev.off()




#  --------------------- determine which cluster to use ------------------------
clust5 <- cluster_ids(sce_tmp_0, "meta5")
clust6 <- cluster_ids(sce_tmp_0, "meta6")
clust8 <- cluster_ids(sce_tmp_0, "meta8")
clust11 <- cluster_ids(sce_tmp_0, "meta11")

sce_tmp_0$meta5 <- cluster_ids(sce_tmp_0, "meta5")
sce_tmp_0$meta8 <- cluster_ids(sce_tmp_0, "meta8")
sce_tmp_0$meta11 <- cluster_ids(sce_tmp_0, "meta11")

c5_dp <- dotplotTables(sce_tmp_0,
                       "meta5",
                       assay = "exprs",
                       fun = "median",
                       scale = TRUE,
                       q = 0.01)

c6_dp <- dotplotTables(sce_tmp_0,
                       "meta6",
                       assay = "exprs",
                       fun = "median",
                       scale = TRUE,
                       q = 0.01)

c8_dp <- dotplotTables(sce_tmp_0,
                       "meta8",
                       assay = "exprs",
                       fun = "median",
                       scale = TRUE,
                       q = 0.01)

c10_dp <- dotplotTables(sce_tmp_0,
                       "meta10",
                       assay = "exprs",
                       fun = "median",
                       scale = TRUE,
                       q = 0.01)

c11_dp <- dotplotTables(sce_tmp_0,
                        "meta11",
                        assay = "exprs",
                        fun = "median",
                        scale = TRUE,
                        q = 0.01)


# intra cluster corr
corr_mat_c5 <- cor(c5_dp$expr_wide)
pheatmap(corr_mat_c5,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Correlation between meta5 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

corr_mat_c6 <- cor(c6_dp$expr_wide)
pheatmap(corr_mat_c6,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Correlation between meta6 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

corr_mat_c8 <- cor(c8_dp$expr_wide)
pheatmap(corr_mat_c8,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Correlation between meta8 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

corr_mat_c10 <- cor(c10_dp$expr_wide)
pheatmap(corr_mat_c10,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Correlation between meta10 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

corr_mat_c11 <- cor(c11_dp$expr_wide)
pheatmap(corr_mat_c11,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Correlation between meta10 clusters",
         fontsize_row = 10,
         fontsize_col = 10)







corr_mat <- cor(c5_dp$expr_wide, c8_dp$expr_wide)
pheatmap(corr_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Correlation between meta5 and meta8 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

corr_mat <- cor(c6_dp$expr_wide, c8_dp$expr_wide)
pheatmap(corr_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Correlation between meta6 and meta8 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

corr_mat <- cor(c6_dp$expr_wide, c10_dp$expr_wide)
pheatmap(corr_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Correlation between meta6 and meta10 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

corr_mat <- cor(c6_dp$expr_wide, c11_dp$expr_wide)
pheatmap(corr_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Correlation between meta6 and meta10 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

