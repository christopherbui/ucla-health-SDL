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

source(paste("U:/cdbui/MPE_Cytof/Rscript/", "Rybakowska_cytof_function.R", sep = ""))
source("C:/Users/cdbui/Documents/GitHub/ucla-health-SDL/cytof-pipeline/MPE_Cytof/Rscript/Rybakowska_cytof_function_LT.R")
source("C:/Users/cdbui/Documents/GitHub/ucla-health-SDL/cytof-pipeline/MPE_Cytof/Rscript/Cytof_MPE_Utils.R")

# set working directory
workFolder <- paste("U:/cdbui/", "MPE_Cytof", sep = "")
setwd(workFolder)


# get sample info of all samples & panels
fin_info <- file.path("Ranalysis","mpe_cytof_sampleInfo_022625.txt")
allSampleInfo <- read.delim(fin_info, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
allSampleInfo <- allSampleInfo[, -1]


# ------------------------------------------------------------------------------
# Check quantile of REFERENCE samples in different batches 
#-------------------------------------------------------------------------------


selPanel <- c("TBNK")  #*******
# selPanel <- c("Myeloid")  #*******
# selPanel <- c("Cytokines")  #*******


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





# define ncols for ggplot2 below
ncols <- 4

sel_marker4plot <- panel_info$fcs_desc[1:14]
tmp4plot <- subset(quantiles_ref ,
    is.element(quantiles_ref$Marker,sel_marker4plot))
p_x <- tmp4plot %>%   
    ggplot2::ggplot(aes(x = Sample,
               y = Value,
               color = Marker)) +
    geom_point() +
     ylab(label = sel_marker4plot)+  
    facet_wrap(~ Marker, ncol = ncols, scales = "free_x") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "bottom")

# plot range (0.25,0.99) and median
sel_marker4plot <- panel_info$fcs_desc[1:9]
tmp4plot <- subset(Ftab_quantile_ref,
    is.element(quantiles_ref$Marker,sel_marker4plot))
p_y <- tmp4plot %>%   
    ggplot2::ggplot(aes(x = Sample,
               y = Value.0.5,
               color = Batch)) +
    geom_pointrange(aes(ymin = Value.0.25, ymax = Value.0.99)) +
    facet_wrap(~ Marker, scales = "free_x") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "bottom")




# --------------------------------------------------------------
# Check panel expression in each reference Cytokine/Myeloid/TBNK cells
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

## get channel information for fsom
# lineage_idx <- which(panel_info$marker_class=="lineage")  #****
# lineage_markers <- panel_info$fcs_desc[lineage_idx]  
# function_idx <- which(panel_info$marker_class=="function")  #****
# function_markers <- panel_info$fcs_desc[function_idx] 
# 
# ## setup panel input for CATALYST, requiring 3 columns
# tmp_panel <- panel_info[,c("fcs_colname", "antigen", "marker_class")]
# 
# tmpi <- which(grepl("function", panel_info$marker_class) == TRUE)
# tmp_panel$marker_class[tmpi] <- c("state")
# 
# tmpi <- which(panel_info$marker_class=="lineage")
# tmp_panel$marker_class[tmpi] <- c("type")

#---------------------

# set input folder
bead_norm_dir <- file.path(workFolder, "CYTOF_data", "BeadNorm", selPanel)

# GET ONLY REF FILES; JUST FOR TESTING
all_fcs_files <- list.files(bead_norm_dir,
                            pattern = "_Ref.*\\.fcs$",
                            full.names = TRUE)


# # set input folder
# bead_norm_dir <- file.path(workFolder,"CYTOF_data",
#                            "Gated",selPanel)
# all_fcs_files <- list.files(bead_norm_dir, 
#                             pattern = ".fcs$", 
#                             full.names = TRUE)
# 
# # set input folder
# bead_norm_dir <- file.path(workFolder,"CYTOF_data",
#                            "Gated_2gates",selPanel)
# all_fcs_files <- list.files(bead_norm_dir, 
#                             pattern = ".fcs$", 
#                             full.names = TRUE)


#---------------------


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
write.table(Ftab_out,file.path(panel_output_dir, fout),
            sep="\t",quote=F,row.names=FALSE)



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

# # ------
# fout <- paste(selPanel,"_cleanedGated_arcsinhTransform_qc.txt",sep="")
# write.table(Ftab_out,file.path("Ranalysis",fout),
#             sep="\t",quote=F,row.names=FALSE)
# # -----
# fout <- paste(selPanel,"_cleanedGated2gates_arcsinhTransform_qc.txt",sep="")
# write.table(Ftab_out,file.path("Ranalysis",fout),
#             sep="\t",quote=F,row.names=FALSE)




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



# setup input folder
# bead_norm_dir <- file.path(workFolder,"CYTOF_data",
#                            "Gated",selPanel)
# all_fcs_files <- list.files(bead_norm_dir, 
#                             pattern = ".fcs$", 
#                             full.names = TRUE)
# 
# bead_norm_dir <- file.path(workFolder,"CYTOF_data",
#                            "Cleaned",selPanel)
# all_fcs_files <- list.files(bead_norm_dir, 
#                             pattern = ".fcs$", 
#                             full.names = TRUE)

# -----------------------


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
# Error in .local(x, ...) : No spillover matrix stored in that sample
# Error in .local(x, ...) : No spillover matrix stored in that sample
# Error in .local(x, ...) : No spillover matrix stored in that sample

tmp_channels <- colnames(f)
channel_info <- pData(parameters(f))

sel_idx <- which(grepl("CD45", channel_info$desc) == TRUE)
sel_channel <- channel_info$name[sel_idx]

ggcyto::autoplot(f, "Pt198Di")



# ------------------------------------------------------------------------------
# Signal Cleaning
#-------------------------------------------------------------------------------

# list_panels <- c("TBNK","Myeloid","Cytokines")

# list_panels <- c("TBNK")
# list_panels <- c("Myeloid")
list_panels <- c("Cytokines")


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
# Files outliers detection
# output is the score from fsom
#-------------------------------------------------------------------------------

# list_panels <- c("TBNK","Myeloid","Cytokines")

# list_panels <- c("TBNK")
list_panels <- c("Myeloid")
# list_panels <- c("Cytokines")

batch4extractPatt <- "(?i).*(batch[0-9]*).*.fcs"
dna_ch_4fsom <- c("191Ir","193Ir")   #*****using desc column
via_ch_4fsom <- c("195Pt")    #***Cisplatin
# nCells_thres <- c(10000,3000,1000)   #*****
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
# Normalization using reference sample
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
#-------------------------------------------------------------------------------
# Normally, we have to change $FIL for normalized gated files.
# But we will use only gated data for changing $FIL & downstream flowSOM & UMAP; Problem with some reference files


# update '$FIL' in fcs files --> output in the subfolder, "updateFileName"
# Set input directory 
# selPanel <- c("TBNK")
selPanel <- c("Myeloid")
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
# Using gated, non-normalized data*******
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

# set output directory
analysis_dir <- file.path(workFolder, "CYTOF_data", "Analysis", selPanel)
if (!dir.exists(analysis_dir)) dir.create(analysis_dir, recursive = TRUE)

# get fcs files
fcs_files <- list.files(norm_dir,
                        pattern = ".fcs$",
                        full.names = TRUE)

# extract sample info from fcs files
getCoreID <- function(x) {
  xlist <- unlist(strsplit(x, split = "_"))[1:4]
  xcore <- paste0(xlist, collapse = "_")
}

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
                          md_cols = list(file = "file",
                                         id = "sample_id",
                                         factors = c("tissue_type", "patient_id", "parental_sample_id", "BATCH")))

SummarizedExperiment(sce)
# class: SummarizedExperiment 
# dim: 41 14001111 
# metadata(0):
#   assays(1): ''
# rownames(41): VISTA CD27 ... CD45 CD11B
# rowData names(0):
#   colnames: NULL
# colData names(0):

colnames(rowData(sce))  # features / markers
dim(colData(sce))       # summary samples data

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


################################################################################
# ONLY REFERENCE ANALYSIS
################################################################################

# inspect data
# inspect Ref-PBMC samples
ref_sam <- grep("Ref", md_info$sample_id, value = TRUE)
sce_ref <- sce[, sce$sample_id %in% ref_sam]  # get data only belonging to Ref samples
colData(sce_ref)$sample_id <- droplevels(colData(sce_ref)$sample_id)

# set output directory for reference file analysis
analysis_ref_dir <- file.path(analysis_dir, "Res_Ref")
if (!dir.exists(analysis_ref_dir)) dir.create(analysis_ref_dir)

# plot distribution of markers
# sce_ref <- sce[, sce$sample_id %in% ref_sam[2]] # selects only batch 2
color_by <- "BATCH"
marker_dist_PLOT <- plotExprs(sce_ref, color_by = color_by)

file_name <- paste0("Ref_marker_distribution_by_", color_by, ".png")
ggsave(filename = file.path(analysis_ref_dir, file_name), plot = marker_dist_PLOT, width = 10, height = 7)



# plot mulit-dimension scale (median marker intensities)
color_by <- "BATCH"
label_by <- "sample_id"
pb_mds_PLOT <- pbMDS(sce_ref, color_by = color_by, label_by = label_by)

file_name <- paste0("Ref_pbMDS_by_", color_by, "_", label_by, ".png")
ggsave(filename = file.path(analysis_ref_dir, file_name), plot = pb_mds_PLOT, width = 8, height = 6)



# get only lineage markers
type_markers <- rownames(sce_ref)[rowData(sce_ref)$marker_class == "type"]


# heatmap of median marker intensities
# in CyTOF applications, a cosine distance shows good performance
# scale = "last": aggregate then scale & trim
marker_heatmap_PLOT <- plotExprHeatmap(sce_ref,
                                       scale = "last",
                                       row_anno = FALSE,
                                       col_anno = FALSE)

file_name <- paste0("Ref_expr_heatmap_all_markers", ".png")
png(file.path(analysis_ref_dir, file_name), width = 14, height = 8, units = "in", res = 300)
draw(marker_heatmap_PLOT)
dev.off()


# identification with FlowSOM: using k = 20 metacluster
# id of som(i.e. 100 clusters) and metacluster are stored in cluster_ids
set.seed(1234)

features = "type"   # type, state, or vector of specific markers
sce_ref <- cluster(sce_ref,
                   features = features,
                   xdim = 10,
                   ydim = 10,
                   maxK = 20,
                   seed = 1234)

# add meta clusters to column data
sce_ref$meta8 <- cluster_ids(sce_ref, "meta8")
sce_ref$meta10 <- cluster_ids(sce_ref, "meta10")
sce_ref$meta15 <- cluster_ids(sce_ref, "meta15")
sce_ref$meta20 <- cluster_ids(sce_ref, "meta20")


marker <- "LINEAGE"  # lineage, functional, all
k <- "meta8"   # meta8, meta10, meta15, meta20
cluster_marker_heatmap <- plotExprHeatmap(sce_ref,
                                          features = "type",
                                          by = "cluster_id",
                                          k = k,
                                          bars = TRUE,
                                          perc = TRUE)

file_name <- paste0(k, "_expr_heatmap_", marker, "_markers", ".png")
png(file.path(analysis_ref_dir, file_name), width = 12, height = 8, units = "in", res = 300)
draw(cluster_marker_heatmap)
dev.off()


names(cluster_codes(sce_ref))    # codes of clusters

# access specific clustering resolution i.e. cellIDS
table(cluster_ids(sce_ref, "som100"))
table(cluster_ids(sce_ref, "meta20"))



# export matrices
k <- "meta20"   # cluster_id (100 clusters), meta8, meta10, meta15, meta20
fun <- "mean"   # mean, median
raw_clust_expr <- CATALYST:::.agg(sce_ref, by = k, assay = "exprs", fun = fun)
raw_clust_expr_df <- as.data.frame(raw_clust_expr)

file_name <- paste0(k, "_", fun, "_expr_matrix", ".txt")
write.table(raw_clust_expr_df, file.path(analysis_ref_dir, file_name), sep = "\t", quote = FALSE, col.names = NA)



# run PCA
pc_prop <- 0.75 # proportion of total features to use
sce_ref <- runPCA(sce_ref, exprs_values = "exprs", ncomponents = floor(length(type_markers) * pc_prop))

# elbow plot
pca_matrix <- reducedDim(sce_ref, "PCA")
pca_variance <- apply(pca_matrix, 2, var)
pca_variance_explained <- pca_variance / sum(pca_variance)

# examine cumulative variance
cumsum(pca_variance_explained)


elbow_df <- data.frame(
  PC = seq_along(pca_variance_explained),
  VarianceExplained = pca_variance_explained
)

elbow_PLOT <- ggplot(elbow_df, aes(x = PC, y = VarianceExplained)) +
  geom_point(color = "#4287f5", size = 2) +
  geom_line(color = "#4287f5", linewidth = 0.9) +
  geom_text(aes(label = round(VarianceExplained, 2)),  # label each point
            hjust = -0.5, vjust = -0.5, size = 3.5) +
  scale_x_continuous(breaks = elbow_df$PC) +           # show all PC numbers on axis
  labs(title = paste0("PCA Elbow Plot", " - ", selPanel, " Reference"),
       x = "Principal Component",
       y = "Proportion of Variance Explained") +
  theme_light(base_size = 14)

file_name <- paste0(selPanel, "_PCA_Elbow_Plot", ".png")
ggsave(filename = file.path(analysis_ref_dir, file_name), plot = elbow_PLOT, height = 8, width = 12)



# UMAP
n_components <- 10   # change as needed

umap_dir <- file.path(analysis_ref_dir, "UMAP", paste0("PC_", n_components))
if (!dir.exists(umap_dir)) dir.create(umap_dir, recursive = TRUE)


sce_ref <- runDR(sce_ref, "UMAP", cells = 1e3, features = type_markers, ncomponents = n_components)

k <- "meta20"  # meta8, meta10, meta15, meta20

umap_PLOT <- plotDR(sce_ref, "UMAP", color_by = k)
file_name <- paste0("UMAP_", k, ".png")
ggsave(filename = file.path(umap_dir, file_name), plot = umap_PLOT, width = 10, height = 10, bg = "white")

umap_facet_PLOT <- plotDR(sce_ref, "UMAP", color_by = k, facet_by = k)
file_name <- paste0("UMAP_", k, "_facet_by_cluster", ".png")
ggsave(filename = file.path(umap_dir, file_name), plot = umap_facet_PLOT, width = 12, height = 10, bg = "white")

umap_facet_sample_PLOT <- plotDR(sce_ref, "UMAP", color_by = k, facet_by = "sample_id")
file_name <- paste0("UMAP_", k, "_facet_by_sample_id", ".png")
ggsave(filename = file.path(umap_dir, file_name), plot = umap_facet_sample_PLOT, width = 14, height = 14, bg = "white")


# dotplot marker expression
dotplot_dir <- file.path(analysis_ref_dir, "Dotplot")
if (!dir.exists(dotplot_dir)) dir.create(dotplot_dir, recursive = TRUE)

k <- "meta8"
fun <- "mean"
dot_PLOT <- dotplot(sce_ref, k = k, fun = fun, output_dir = dotplot_dir)

file_name <- paste0(k, "dotplot_", fun, ".png")
ggsave(filename = file.path(dotplot_dir, file_name), plot = dot_PLOT, width = 11, height = 6)



# analyze cluster proportion by sample
cluster_info_dir <- file.path(analysis_ref_dir, "Cluster_Info")
if (!dir.exists(cluster_info_dir)) dir.create(cluster_info_dir, recursive = TRUE)

k <- "meta8"
df_clusters_per_sample <- get_cluster_prop_within_sample(sce_ref, k)

file_name <- paste0(k, "_cluster_proportion_within_sample.txt")
write.table(df_clusters_per_sample, file.path(cluster_info_dir, file_name), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# track meta cluster grouping
meta_vec <- c("meta8", "meta10", "meta15", "meta20")  # change as needed
df_cluster_map <- get_cluster_mapping(sce_ref, meta_vec = meta_vec)

file_name <- paste0("cluster_mapping.txt")
write.table(df_cluster_map, file.path(cluster_info_dir, file_name), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



################################################################################
# ANALYSIS - ALL SAMPLES
################################################################################

# plot distribution of markers
color_by <- "BATCH"
marker_dist_PLOT <- plotExprs(sce, color_by = color_by)

file_name <- paste0("marker_distribution_by_", color_by, ".png")
ggsave(filename = file.path(analysis_dir, file_name), plot = marker_dist_PLOT, width = 10, height = 7)

# plot mulit-dimension scale (median marker intensities)
color_by <- "BATCH"
label_by <- "sample_id"
pb_mds_PLOT <- pbMDS(sce, color_by = color_by, label_by = label_by)

file_name <- paste0("pbMDS_by_", color_by, "_", label_by, ".png")
ggsave(filename = file.path(analysis_dir, file_name), plot = pb_mds_PLOT, width = 14, height = 14)


# marker heatmap by sample
marker_heatmap_PLOT <- plotExprHeatmap(sce,
                                       scale = "last",
                                       row_anno = FALSE,
                                       col_anno = FALSE)

file_name <- paste0("expr_heatmap_all_markers", ".png")
png(file.path(analysis_dir, file_name), width = 18, height = 12, units = "in", res = 300)
draw(marker_heatmap_PLOT, heatmap_legend_side = "bottom", padding = unit(c(10, 10, 10, 20), "mm"))
dev.off()


################################################################################
# ANALYSIS - ALL SAMPLES ; LINEAGE MARKERS ONLY
################################################################################

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

# set output directory
res_dir <- file.path(analysis_dir, "Results")
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

# get only lineage markers
type_markers <- rownames(sce)[rowData(sce)$marker_class == "type"]

# subset sce for lineage markers; reassign 'sce' to save memory
sce <- sce[type_markers, ]



# plot distribution of markers
color_by <- "BATCH"
marker_dist_PLOT <- plotExprs(sce, color_by = color_by)

marker <- "LINEAGE"   # LINEAGE, FUNCTIONAL, ""

file_name <- paste0("marker_distribution_by_", color_by, "_", marker, "_markers", ".png")
ggsave(filename = file.path(res_dir, file_name), plot = marker_dist_PLOT, width = 10, height = 7)

# plot mulit-dimension scale (median marker intensities)
color_by <- "BATCH"
label_by <- "sample_id"
pb_mds_PLOT <- pbMDS(sce, color_by = color_by, label_by = label_by)

file_name <- paste0("pbMDS_by_", color_by, "_", label_by, "_", marker, "_markers", ".png")
ggsave(filename = file.path(res_dir, file_name), plot = pb_mds_PLOT, width = 14, height = 14)



maxk <- 25  # change as needed



# flowSOM
set.seed(1234)

sce <- cluster(sce,
               features = type_markers,
               xdim = 10,
               ydim = 10,
               maxK = maxk,
               seed = 1234)


# add meta clusters to column data
sce$meta8 <- cluster_ids(sce, "meta8")
sce$meta10 <- cluster_ids(sce, "meta10")
sce$meta15 <- cluster_ids(sce, "meta15")
sce$meta20 <- cluster_ids(sce, "meta20")
sce$meta25 <- cluster_ids(sce, "meta25")


k <- "meta25"   # meta8, meta10, meta15, meta20, meta25
cluster_marker_heatmap <- plotExprHeatmap(sce,
                                          features = type_markers,
                                          by = "cluster_id",
                                          k = k,
                                          bars = TRUE,
                                          perc = TRUE)

marker <- "LINEAGE"  # lineage, functional, all
file_name <- paste0(k, "_expr_heatmap_", marker, "_markers", ".png")
png(file.path(res_dir, file_name), width = 12, height = 8, units = "in", res = 300)
draw(cluster_marker_heatmap)
dev.off()



# export matrices
k <- "meta25"   # cluster_id (100 clusters), meta8, meta10, meta15, meta20
fun <- "mean"   # mean, median
raw_clust_expr <- CATALYST:::.agg(sce, by = k, assay = "exprs", fun = fun)
raw_clust_expr_df <- as.data.frame(raw_clust_expr)

file_name <- paste0(k, "_", fun, "_expr_matrix", ".txt")
write.table(raw_clust_expr_df, file.path(res_dir, file_name), sep = "\t", quote = FALSE, col.names = NA)




# run PCA
pc_prop <- 0.75 # proportion of total features to use
sce <- runPCA(sce, exprs_values = "exprs", ncomponents = floor(length(type_markers) * pc_prop))

# elbow plot
pca_matrix <- reducedDim(sce, "PCA")
pca_variance <- apply(pca_matrix, 2, var)
pca_variance_explained <- pca_variance / sum(pca_variance)

# examine cumulative variance
cumsum(pca_variance_explained)


elbow_df <- data.frame(
  PC = seq_along(pca_variance_explained),
  VarianceExplained = pca_variance_explained
)

elbow_PLOT <- ggplot(elbow_df, aes(x = PC, y = VarianceExplained)) +
              geom_point(color = "#4287f5", size = 2) +
              geom_line(color = "#4287f5", linewidth = 0.9) +
              geom_text(aes(label = round(VarianceExplained, 2)),  # label each point
                        hjust = -0.5, vjust = -0.5, size = 3.5) +
              scale_x_continuous(breaks = elbow_df$PC) +           # show all PC numbers on axis
              labs(title = paste0("PCA Elbow Plot", " - ", selPanel),
                   x = "Principal Component",
                   y = "Proportion of Variance Explained") +
              theme_light(base_size = 14)

file_name <- paste0(selPanel, "_PCA_Elbow_Plot", ".png")
ggsave(filename = file.path(res_dir, file_name), plot = elbow_PLOT, height = 8, width = 12)




# UMAP
n_components <- 4   # change as needed

umap_dir <- file.path(res_dir, "UMAP", paste0("PC_", n_components))
if (!dir.exists(umap_dir)) dir.create(umap_dir, recursive = TRUE)


sce <- runDR(sce, "UMAP", cells = 1e3, features = type_markers, ncomponents = n_components)

k <- "meta20"  # change as needed

umap_PLOT <- plotDR(sce, "UMAP", color_by = k)
file_name <- paste0("UMAP_", k, ".png")
ggsave(filename = file.path(umap_dir, file_name), plot = umap_PLOT, width = 10, height = 10, bg = "white")

umap_facet_PLOT <- plotDR(sce, "UMAP", color_by = k, facet_by = k)
file_name <- paste0("UMAP_", k, "_facet_by_cluster", ".png")
ggsave(filename = file.path(umap_dir, file_name), plot = umap_facet_PLOT, width = 12, height = 10, bg = "white")

# umap_facet_sample_PLOT <- plotDR(sce, "UMAP", color_by = k, facet_by = "sample_id")
# file_name <- paste0("UMAP_", k, "_facet_by_sample_id", ".png")
# ggsave(filename = file.path(umap_dir, file_name), plot = umap_facet_sample_PLOT, width = 20, height = 20, bg = "white")


# dotplot marker expression
dotplot_dir <- file.path(res_dir, "Dotplot")
if (!dir.exists(dotplot_dir)) dir.create(dotplot_dir, recursive = TRUE)

k <- "meta20"   # change as needed
fun <- "mean"
dot_PLOTS <- dotplot(sce, k = k, fun = fun, scale = FALSE, output_dir = dotplot_dir)

file_name <- paste0(k, "_dotplot_", fun, ".png")
ggsave(filename = file.path(dotplot_dir, file_name), plot = dot_PLOTS$dp, width = 10, height = 8)
# file_name <- paste0(k, "_dotplot_", fun, "_values", ".png")
# ggsave(filename = file.path(dotplot_dir, file_name), plot = dot_PLOTS$dp_val, width = 10, height = 8)


# analyze cluster proportion by sample
cluster_info_dir <- file.path(res_dir, "Cluster_Info")
if (!dir.exists(cluster_info_dir)) dir.create(cluster_info_dir, recursive = TRUE)

k <- "meta25"
df_clusters_per_sample <- get_cluster_prop_within_sample(sce, k)

file_name <- paste0(k, "_cluster_proportion_within_sample.txt")
write.table(df_clusters_per_sample, file.path(cluster_info_dir, file_name), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# track meta cluster grouping
meta_vec <- c("meta8", "meta10", "meta15", "meta20", "meta25")  # change as needed
df_cluster_map <- get_cluster_mapping(sce, meta_vec = meta_vec)

file_name <- paste0("cluster_mapping.txt")
write.table(df_cluster_map, file.path(cluster_info_dir, file_name), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


#-------------------------------------------------------------------------------
# ANALYSIS - PIPELINE WITH FSOM IN PC = 6
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

# set output directory
res_dir <- file.path(analysis_dir, "Results")
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


#-------------------------------------------------------------------------------
# FULL Lineage Markers
type_markers <- rownames(sce)[rowData(sce)$marker_class == "type"]

# subset sce for markers; reassign 'sce' to save memory
sce <- sce[type_markers, ]
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# SUBSET Lineage Markers
sce_lin_subset <- CATALYST::prepData(fcs_files,
                          panel = tmp_panel, # panel marker info
                          md = md_info, # sample info
                          panel_cols = list(channel = "fcs_colname", antigen = "antigen", class = "marker_class"),
                          md_cols = list(file = "file", id = "sample_id", factors = c("tissue_type", "patient_id", "parental_sample_id", "BATCH")))

# add "Symbol" marker names for plotting; preferred for publishing
matched_symbols <- panel_info$Symbol[match(rownames(sce_lin_subset), panel_info$antigen)]
rowData(sce_lin_subset)$Symbol <- matched_symbols


subset_markers <- readLines("U:/cdbui/MPE_Cytof/CYTOF_data/Analysis/TBNK/Results_Subset_Lineage/subset_lineage_markers.txt")
sce_lin_subset <- sce_lin_subset[rowData(sce_lin_subset)$Symbol %in% subset_markers, ]
#-------------------------------------------------------------------------------

# run PCA
pc_prop <- 0.75 # proportion of total features to use
sce <- runPCA(sce, exprs_values = "exprs", ncomponents = floor(length(type_markers) * pc_prop))

# elbow plot
pca_matrix <- reducedDim(sce, "PCA")
pca_variance <- apply(pca_matrix, 2, var)
pca_variance_explained <- pca_variance / sum(pca_variance)

# examine cumulative variance
cumsum(pca_variance_explained)


elbow_df <- data.frame(
  PC = seq_along(pca_variance_explained),
  VarianceExplained = pca_variance_explained
)

elbow_PLOT <- ggplot(elbow_df, aes(x = PC, y = VarianceExplained)) +
  geom_point(color = "#4287f5", size = 2) +
  geom_line(color = "#4287f5", linewidth = 0.9) +
  geom_text(aes(label = round(VarianceExplained, 2)),  # label each point
            hjust = -0.5, vjust = -0.5, size = 3.5) +
  scale_x_continuous(breaks = elbow_df$PC) +           # show all PC numbers on axis
  labs(title = paste0("PCA Elbow Plot", " - ", selPanel),
       x = "Principal Component",
       y = "Proportion of Variance Explained") +
  theme_light(base_size = 14)

file_name <- paste0("PCA_Elbow_Plot", ".png")
ggsave(filename = file.path(res_dir, file_name), plot = elbow_PLOT, height = 8, width = 12)


# flowSOM by Catalyst; done on full set of marker dimensions
n_components <- 10   # change as needed; derived from elbow plot

cluster_info_dir <- file.path(res_dir, "Cluster_Info")
if (!dir.exists(cluster_info_dir)) dir.create(cluster_info_dir, recursive = TRUE)

cluster_info_pc_dir <- file.path(cluster_info_dir, paste0("PC_", n_components))
if (!dir.exists(cluster_info_pc_dir)) dir.create(cluster_info_pc_dir, recursive = TRUE)

set.seed(1234)
maxK <- 20

sce <- cluster(sce,
               features = type_markers,
               xdim = 10,
               ydim = 10,
               maxK = maxK,
               seed = 1234)
# add meta clusters to column data
sce$meta20 <- cluster_ids(sce, "meta20")

# flowSOM on PCA
# fsom_dir <- file.path(res_dir, "FSOM", paste0("PC_", n_components))
# if (!dir.exists(fsom_dir)) dir.create(fsom_dir, recursive = TRUE)

# select data with desired PC
pca_data <- reducedDim(sce, "PCA")
pca_data_selected <- as.matrix(pca_data[, 1:n_components])

fsom <- FlowSOM::FlowSOM(
  input = pca_data_selected,
  xdim = 10,
  ydim = 10,
  seed = 1234
)
# initial xdim x ydim SOM nodes
cell_clusters <- fsom$map$mapping[, 1]
colData(sce)[["cluster_id_processed"]] <- as.factor(cell_clusters)

# perform metaclustering
meta_cluster_map <- metaClustering_consensus(fsom$map$codes, k = maxK)

# map processed meta cluster to processed cluster id
meta_cluster_labels <- meta_cluster_map[cell_clusters]

# assign processed metacluster id to cells
meta_cluster_name <- paste0("meta", maxK, "_processed")
colData(sce)[[meta_cluster_name]] <- as.factor(meta_cluster_labels)

# sce column data
df <- as.data.frame(colData(sce))

# unique raw & processed cluster mappings
cluster_cols <- c("cluster_id", "meta20", "cluster_id_processed", "meta20_processed")
cluster_map <- unique(df[, cluster_cols])
cluster_map <- cluster_map[order(cluster_map$cluster_id), ]
rownames(cluster_map) <- NULL
write.table(cluster_map, file = file.path(cluster_info_pc_dir, "meta20_raw_processed_mapping.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# raw cluster mapping
raw_cols <- c("som100", "meta20")
raw_map <- cluster_codes(sce)[, raw_cols]
names(raw_map)[1] <- "cluster_id"
write.table(raw_map, file = file.path(cluster_info_dir, "meta20_raw_mapping.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# processed cluster mapping
proc_cols <- c("cluster_id_processed", "meta20_processed")
proc_map <- unique(df[, proc_cols])
proc_map <- proc_map[order(proc_map$cluster_id), ]
rownames(proc_map) <- NULL
write.table(proc_map, file = file.path(cluster_info_pc_dir, "meta20_processed_mapping.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# view cross cluster table
cross_cluster_table <- table(meta20 = sce$meta20, meta20_processed = sce$meta20_processed)
write.table(as.data.frame.matrix(cross_cluster_table), file = file.path(cluster_info_pc_dir, "meta20_cross_cluster_count.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# row-wise proportion calculation
cross_cluster_table_prop <- prop.table(cross_cluster_table, margin = 1)
write.table(as.data.frame.matrix(cross_cluster_table_prop), file = file.path(cluster_info_pc_dir, "meta20_cross_cluster_prop.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# heatmap of cluster prop table
cross_clust_prop_mat <- as.matrix(cross_cluster_table_prop)
hm <- Heatmap(cross_clust_prop_mat,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = TRUE,
              row_names_side = "left",
              show_column_names = TRUE,
              column_names_side = "bottom",
              column_names_rot = 0,
              row_title = "Raw Clusters",
              column_title = "Processed Clusters",
              column_title_side = "bottom")
png(file.path(cluster_info_pc_dir, "meta20_cross_cluster_prop_heatmap.png"), width = 3000, height = 3000, res = 300)
draw(hm, column_title = "Raw-Processed Cluster Proportions")
dev.off()


cell_dist_raw_count <- table(cluster_id = sce$cluster_id, meta20 = sce$meta20)
cell_dist_raw_prop <- prop.table(cell_dist_raw_count, margin = 1)
write.table(cell_dist_raw_count, file = file.path(cluster_info_dir, "cell_dist_raw_count.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(cell_dist_raw_prop, file = file.path(cluster_info_dir, "cell_dist_raw_prop.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

cell_dist_proc_count <- table(cluster_id_processed = sce$cluster_id_processed, meta20_processed = sce$meta20_processed)
cell_dist_proc_prop <- prop.table(cell_dist_proc_count, margin = 1)
write.table(cell_dist_proc_count, file = file.path(cluster_info_pc_dir, "cell_dist_proc_count.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(cell_dist_proc_prop, file = file.path(cluster_info_pc_dir, "cell_dist_proc_prop.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)


# UMAP
n_components <- 10

umap_dir <- file.path(res_dir, "UMAP")
if (!dir.exists(umap_dir)) dir.create(umap_dir, recursive = TRUE)

umap_pc_dir <- file.path(umap_dir, paste0("PC_", n_components))
if (!dir.exists(umap_pc_dir)) dir.create(umap_pc_dir, recursive = TRUE)

sce <- runDR(sce, "UMAP", cells = 1e3, features = type_markers, ncomponents = n_components, seed = 1234)

cluster_name <- "meta20_processed"

umap_PLOT <- plotDR(sce, "UMAP", color_by = cluster_name)
file_name <- paste0("UMAP_", cluster_name, ".png")
ggsave(filename = file.path(umap_pc_dir, file_name), plot = umap_PLOT, width = 10, height = 10, bg = "white")

umap_facet_PLOT <- plotDR(sce, "UMAP", color_by = cluster_name, facet_by = cluster_name)
file_name <- paste0("UMAP_", cluster_name, "_facet_by_cluster", ".png")
ggsave(filename = file.path(umap_pc_dir, file_name), plot = umap_facet_PLOT, width = 12, height = 10, bg = "white")


# dotplot

# create directories
dotplot_dir <- file.path(res_dir, "Dotplot")
if (!dir.exists(dotplot_dir)) dir.create(dotplot_dir, recursive = TRUE)

dotplot_no_scaled_dir <- file.path(dotplot_dir, "Test/No_Scaled")
if (!dir.exists(dotplot_no_scaled_dir)) dir.create(dotplot_no_scaled_dir, recursive = TRUE)

dotplot_scaled_dir <- file.path(dotplot_dir, "Scaled")
if (!dir.exists(dotplot_scaled_dir)) dir.create(dotplot_scaled_dir, recursive = TRUE)

# create raw & processed for scaled and no_scaled
sel_dir <- dotplot_no_scaled_dir  # change as needed: dotplot_no_scaled_dir, dotplot_scaled_dir

dotplot_raw_dir <- file.path(sel_dir, "Raw")
if (!dir.exists(dotplot_raw_dir)) dir.create(dotplot_raw_dir, recursive = TRUE)

dotplot_proc_dir <- file.path(sel_dir, "Processed", paste0("PC_", n_components))
if (!dir.exists(dotplot_proc_dir)) dir.create(dotplot_proc_dir, recursive = TRUE)


# do dotplot
cluster_name <- "meta20"  # meta20, meta20_processed
fun <- "median"

# raw, no scale
dot_PLOT <- dotplot_OG(sce, cluster_name = cluster_name, fun = fun, scale = FALSE, output_dir = dotplot_raw_dir)
file_name <- paste0(cluster_name, "_dotplot_", fun, ".png")
ggsave(filename = file.path(dotplot_raw_dir, file_name), plot = dot_PLOT, width = 10, height = 8)

# processed, no scale
dot_PLOT <- dotplot_OG(sce, cluster_name = cluster_name, fun = fun, scale = FALSE, output_dir = dotplot_proc_dir)
file_name <- paste0(cluster_name, "_dotplot_", fun, ".png")
ggsave(filename = file.path(dotplot_proc_dir, file_name), plot = dot_PLOT, width = 10, height = 8)




# make original scaling dir
dotplot_scaled_og_dir <- file.path(dotplot_scaled_dir, "Original_Scaling")
if (!dir.exists(dotplot_scaled_og_dir)) dir.create(dotplot_scaled_og_dir, recursive = TRUE)

# make raw & processed dirs
dotplot_scaled_og_raw_dir <- file.path(dotplot_scaled_og_dir, "Raw")
if (!dir.exists(dotplot_scaled_og_raw_dir)) dir.create(dotplot_scaled_og_raw_dir, recursive = TRUE)

dotplot_scaled_og_proc_dir <- file.path(dotplot_scaled_og_dir, "Processed")
if (!dir.exists(dotplot_scaled_og_proc_dir)) dir.create(dotplot_scaled_og_proc_dir, recursive = TRUE)

# raw, scale OG
dot_PLOT <- dotplot_OG(sce, cluster_name = cluster_name, fun = fun, scale = TRUE, output_dir = dotplot_scaled_og_raw_dir)
file_name <- paste0(cluster_name, "_dotplot_", fun, ".png")
ggsave(filename = file.path(dotplot_scaled_og_raw_dir, file_name), plot = dot_PLOT, width = 10, height = 8)

# processed, scale OG
dot_PLOT <- dotplot_OG(sce, cluster_name = cluster_name, fun = fun, scale = TRUE, output_dir = dotplot_scaled_og_proc_dir)
file_name <- paste0(cluster_name, "_dotplot_", fun, ".png")
ggsave(filename = file.path(dotplot_scaled_og_proc_dir, file_name), plot = dot_PLOT, width = 10, height = 8)



# make qq scaling dir
dotplot_scaled_qq_dir <- file.path(dotplot_scaled_dir, "QQ_Scaling")

# make raw & processed dirs
dotplot_scaled_qq_raw_dir <- file.path(dotplot_scaled_qq_dir, "Raw")
if (!dir.exists(dotplot_scaled_qq_raw_dir)) dir.create(dotplot_scaled_qq_raw_dir, recursive = TRUE)

dotplot_scaled_qq_proc_dir <- file.path(dotplot_scaled_qq_dir, "Processed")
if (!dir.exists(dotplot_scaled_qq_proc_dir)) dir.create(dotplot_scaled_qq_proc_dir, recursive = TRUE)

# raw, scale qq
dot_PLOT <- dotplot_4SDL(sce, cluster_name = cluster_name, fun = fun, scale = TRUE, output_dir = dotplot_scaled_qq_raw_dir)
file_name <- paste0(cluster_name, "_dotplot_", fun, ".png")
ggsave(filename = file.path(dotplot_scaled_qq_raw_dir, file_name), plot = dot_PLOT, width = 10, height = 8)

# processed, scale qq
dot_PLOT <- dotplot_4SDL(sce, cluster_name = cluster_name, fun = fun, scale = TRUE, output_dir = dotplot_scaled_qq_proc_dir)
file_name <- paste0(cluster_name, "_dotplot_", fun, ".png")
ggsave(filename = file.path(dotplot_scaled_qq_proc_dir, file_name), plot = dot_PLOT, width = 10, height = 8)








dot_PLOT <- dotplot_4SDL(sce, cluster_name = cluster_name, fun = fun, scale = scale, output_dir = dotplot_dir)
file_name <- paste0(cluster_name, "_dotplot_", fun, ".png")
ggsave(filename = file.path(dotplot_dir, file_name), plot = dot_PLOT, width = 10, height = 8)


# test_dir_scaled <- file.path(dotplot_dir, "TEST/Scaled")
# if (!dir.exists(test_dir_scaled)) dir.create(test_dir_scaled, recursive = TRUE)
# 
# test_dir_no_scaled <- file.path(dotplot_dir, "TEST/No_Scaled")
# if (!dir.exists(test_dir_no_scaled)) dir.create(test_dir_no_scaled, recursive = TRUE)

# k <- "meta20" # meta20, meta20_processed
# fun <- "mean"
# dot_PLOT_TEST <- dotplot_4SDL(sce, cluster_name = k, fun = fun, scale = TRUE, output_dir = test_dir_no_scaled)
# file_name <- paste0(k, "_dotplot_", fun, "_TEST", ".png")
# ggsave(filename = file.path(test_dir_no_scaled, file_name), plot = dot_PLOT, width = 10, height = 8)


###############################
# TEST DOTPLOT
###############################
output_dir <- dotplot_raw_dir

cluster_name <- "meta20"  # meta20, meta20_processed
fun <- "median"
scale <- FALSE
q <- 0.01
assay <- "exprs"

sce$cid_sel <- colData(sce)[[cluster_name]]

es <- assay(sce, assay) # expression matrix
th <- rowMedians(es)  # median threshold

cs <- seq_len(ncol(sce))  # indices vector with length = number of cells
cs <- split(cs, sce$cid_sel)   # split indices into groups based on cid_sel

# return fraction of cells > th per marker
fq <- sapply(cs, function(i) {
  rowMeans(es[, i, drop = FALSE] > th)
})

# compute mean/median expression by cluster
lab <- paste(fun, assay)
ms <- CATALYST:::.agg(sce, by = "cid_sel", assay = assay, fun = fun)

# use "Symbol" marker name
stopifnot(all(rownames(ms) %in% rownames(sce)))   # check
marker_symbols <- rowData(sce)$Symbol[match(rownames(ms), rownames(sce))]
rownames(ms) <- marker_symbols

# save un-scaled ms
write.table(ms, file = file.path(output_dir, paste0(cluster_name, "_dotplot_", "NO_SCALED_expr_wide.txt")), sep = "\t", quote = FALSE, col.names = NA)

# ......

# scaling
lab <- paste("scaled", lab)
ms <- CATALYST:::.scale_exprs(ms, q = q)
# save scaled ms
write.table(ms, file = file.path(output_dir, paste0(cluster_name, "_dotplot_", "SCALED_expr_wide.txt")), sep = "\t", quote = FALSE, col.names = NA)




#-------------------------------------------------------------------------------
# ANALYSIS - DETERMINE WHETHER FSOM CLUSTERING ON FULL MARKER DIM OR PC DIM IS BETTER/CONSISTENT
#-------------------------------------------------------------------------------
no_scaled_output_dir <- file.path(res_dir, "Dotplot/PC_6/TEST/No_Scaled")

# load cross cluster prop matrix
cross_clust_table_prop <- read.delim(file.path(cluster_info_dir, "meta20_cross_cluster_prop.txt"), header = TRUE, check.names = FALSE)
cross_clust_prop_mat <- as.matrix(cross_clust_table_prop)[, -1]

# load raw & processed dotplot matrices
# marker x cluster
raw_dp <- read.delim(file.path(res_dir, "Dotplot/PC_6/meta20_dotplot_mean_expr_matrix_wide.txt"), header = TRUE, check.names = FALSE)
proc_dp <- read.delim(file.path(res_dir, "Dotplot/PC_6/meta20_processed_dotplot_mean_expr_matrix_wide.txt"), header = TRUE, check.names = FALSE)

# only get numeric columns for correlation calculation
raw_dp <- raw_dp[, -1]
proc_dp <- proc_dp[, -1]


# raw_cluster x processed_cluster
dp_corr_mat <- cor(raw_dp, proc_dp)
# NOTE: cluster_info_dir
write.table(dp_corr_mat, file = file.path(no_scaled_output_dir, paste0("meta", maxK, "_cross_cluster_dotplot_corr.txt")), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
hm <- Heatmap(dp_corr_mat,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_names_side = "left",
              row_title = "Raw Clusters",
              column_names_side = "bottom",
              column_names_rot = 0,
              column_title = "Processed Clusters",
              column_title_side = "bottom")
png(file.path(no_scaled_output_dir, "cross_cluster_dotplot_correlation.png"), width = 3000, height = 3000, res = 300)
draw(hm, column_title = "Raw-Processed Cluster Correlation (from Dotplots)")
dev.off()


#-------------------------------------------------------------------------------
# ANALYSIS - DETERMINE WHETHER FSOM CLUSTERING ON FULL MARKER DIM OR PC DIM IS BETTER/CONSISTENT
# SCALED DATA
#-------------------------------------------------------------------------------
scaled_output_dir <- file.path(res_dir, "Dotplot/PC_6/TEST/Scaled")
maxK <- 20

# load cross cluster prop matrix
cross_clust_table_prop <- read.delim(file.path(cluster_info_dir, "meta20_cross_cluster_prop.txt"), header = TRUE, check.names = FALSE)
cross_clust_prop_mat <- as.matrix(cross_clust_table_prop)[, -1]

# load raw & processed dotplot matrices
# marker x cluster
raw_dp <- read.delim(file.path(res_dir, "Dotplot/PC_6/TEST/Scaled/meta20_dotplot_mean_expr_matrix_wide.txt"), header = TRUE, check.names = FALSE)
proc_dp <- read.delim(file.path(res_dir, "Dotplot/PC_6/TEST/Scaled/meta20_processed_dotplot_mean_expr_matrix_wide.txt"), header = TRUE, check.names = FALSE)

# order by markers for consistency
# common_markers <- intersect(raw_dp$Symbol, proc_dp$Symbol)
# raw_dp <- raw_dp[order(raw_dp$Symbol), ]
# proc_dp <- proc_dp[order(proc_dp$Symbol), ]

# only get numeric columns for correlation calculation
raw_dp <- raw_dp[, -1]
proc_dp <- proc_dp[, -1]


# correlation matrix of raw & processed dotplot clusters
# marker x cluster, so we do corr by column

# raw_cluster x processed_cluster
dp_corr_mat <- cor(raw_dp, proc_dp)
write.table(dp_corr_mat, file = file.path(scaled_output_dir, paste0("meta", maxK, "_cross_cluster_dotplot_corr.txt")), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
hm <- Heatmap(dp_corr_mat,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_names_side = "left",
              row_title = "Raw Clusters",
              column_names_side = "bottom",
              column_names_rot = 0,
              column_title = "Processed Clusters",
              column_title_side = "bottom")
png(file.path(scaled_output_dir, "cross_cluster_dotplot_correlation.png"), width = 3000, height = 3000, res = 300)
draw(hm, column_title = "Raw-Processed Cluster Correlation (from Dotplots)")
dev.off()


#-------------------------------------------------------------------------------
# ANALYSIS - RAW & PROCESSED DOTPLOT CORRELATION
#-------------------------------------------------------------------------------
cluster_name <- "meta20"

no_scaled_corr_dir <- file.path(dotplot_no_scaled_dir, "Corr_Info")
if (!dir.exists(no_scaled_corr_dir)) dir.create(no_scaled_corr_dir, recursive = TRUE)

# load cross cluster prop matrix
cross_clust_table_prop <- read.delim(file.path(cluster_info_dir, "PC_10", "meta20_cross_cluster_prop.txt"), header = TRUE, check.names = FALSE)
cross_clust_prop_mat <- as.matrix(cross_clust_table_prop)[, -1]

# load raw & processed dotplot matrices
# marker x cluster
raw_dp <- read.delim(file.path(dotplot_no_scaled_dir, "Raw/meta20_dotplot_NO_SCALED_expr_wide.txt"), header = TRUE, check.names = FALSE)
proc_dp <- read.delim(file.path(dotplot_no_scaled_dir, "Processed/PC_10/meta20_processed_dotplot_NO_SCALED_expr_wide.txt"), header = TRUE, check.names = FALSE)

# order by markers for consistency
# common_markers <- intersect(raw_dp$Symbol, proc_dp$Symbol)
# raw_dp <- raw_dp[order(raw_dp$Symbol), ]
# proc_dp <- proc_dp[order(proc_dp$Symbol), ]

# only get numeric columns for correlation calculation
raw_dp <- raw_dp[, -1]
proc_dp <- proc_dp[, -1]


# correlation matrix of raw & processed dotplot clusters
# marker x cluster, so we do corr by column

# raw_cluster x processed_cluster
dp_corr_mat <- cor(raw_dp, proc_dp)
# NOTE: cluster_info_dir
write.table(dp_corr_mat, file = file.path(no_scaled_corr_dir, paste0("meta", maxK, "_cross_cluster_dotplot_corr.txt")), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
hm <- Heatmap(dp_corr_mat,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_names_side = "left",
              row_title = "Raw Clusters",
              column_names_side = "bottom",
              column_names_rot = 0,
              column_title = "Processed Clusters",
              column_title_side = "bottom")
png(file.path(no_scaled_corr_dir, paste0(cluster_name, "_cross_cluster_dotplot_corr.png")), width = 3000, height = 3000, res = 300)
draw(hm, column_title = "Raw-Processed Cluster Correlation (from Dotplots)")
dev.off()



#-------------------------------------------------------------------------------
# ANALYSIS - SUBSET LINEAGE MARKER
#-------------------------------------------------------------------------------
subset_dir <- file.path(analysis_dir, "Results_Subset_Lineage")
if (!dir.exists(subset_dir)) dir.create(subset_dir, recursive = TRUE)

cluster_info_dir <- file.path(subset_dir, "Cluster_Info")
if (!dir.exists(cluster_info_dir)) dir.create(cluster_info_dir, recursive = TRUE)


maxK <- 20

markers <- rowData(sce_lin_subset)$marker_name[rowData(sce_lin_subset)$Symbol %in% subset_markers]
sce_lin_subset <- cluster(sce_lin_subset,
                          features = markers,
                          xdim = 10,
                          ydim = 10,
                          maxK = maxK,
                          seed = 1234)
# add meta clusters to column data
sce_lin_subset$meta20_subset <- cluster_ids(sce_lin_subset, "meta20")




# raw_subset cluster mapping
raw_cols <- c("som100", "meta20")
raw_map <- cluster_codes(sce)[, raw_cols]
names(raw_map)[1] <- "cluster_id_subset"
write.table(raw_map, file = file.path(cluster_info_dir, "meta20_raw_mapping_subset.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# raw cell distribution
cell_dist_raw_count <- table(cluster_id = sce$cluster_id, meta20 = sce$meta20)
cell_dist_raw_prop <- prop.table(cell_dist_raw_count, margin = 1)
write.table(cell_dist_raw_count, file = file.path(cluster_info_dir, "cell_dist_raw_count_subset.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(cell_dist_raw_prop, file = file.path(cluster_info_dir, "cell_dist_raw_prop_subset.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)


# raw & raw_subset cross distribution
cross_cluster_table <- table(meta20 = sce$meta20, meta20_subset = sce$meta20_subset)
write.table(as.data.frame.matrix(cross_cluster_table), file = file.path(cluster_info_dir, "TEST_meta20_raw_raw_subset_count.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# row-wise proportion calculation
cross_cluster_table_prop <- prop.table(cross_cluster_table, margin = 1)
write.table(as.data.frame.matrix(cross_cluster_table_prop), file = file.path(cluster_info_dir, "TEST_meta20_raw_raw_subset_prop.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# heatmap of cluster prop table
cross_clust_prop_mat <- as.matrix(cross_cluster_table_prop)
hm <- Heatmap(cross_clust_prop_mat,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = TRUE,
              row_names_side = "left",
              show_column_names = TRUE,
              column_names_side = "bottom",
              column_names_rot = 0,
              row_title = "Raw Meta20",
              column_title = "Raw Meta20 Subset",
              column_title_side = "bottom")
png(file.path(cluster_info_dir, "TEST_meta20_raw_raw_subset_prop_heatmap.png"), width = 3000, height = 3000, res = 300)
draw(hm, column_title = "Raw - Raw Subset Cluster Proportions")
dev.off()




# dotplot
# create directories
dotplot_dir <- file.path(subset_dir, "Dotplot")
if (!dir.exists(dotplot_dir)) dir.create(dotplot_dir, recursive = TRUE)

dotplot_no_scaled_dir <- file.path(dotplot_dir, "Test/No_Scaled")
if (!dir.exists(dotplot_no_scaled_dir)) dir.create(dotplot_no_scaled_dir, recursive = TRUE)

# dotplot_scaled_dir <- file.path(dotplot_dir, "Scaled")
# if (!dir.exists(dotplot_scaled_dir)) dir.create(dotplot_scaled_dir, recursive = TRUE)

sel_dir <- dotplot_no_scaled_dir  # change as needed: dotplot_no_scaled_dir, dotplot_scaled_dir

dotplot_raw_dir <- file.path(sel_dir, "Raw")
if (!dir.exists(dotplot_raw_dir)) dir.create(dotplot_raw_dir, recursive = TRUE)

# dotplot_proc_dir <- file.path(sel_dir, "Processed", paste0("PC_", n_components))
# if (!dir.exists(dotplot_proc_dir)) dir.create(dotplot_proc_dir, recursive = TRUE)


# do dotplot
cluster_name <- "meta20_subset"  # meta20_subset
fun <- "median"

# raw, no scale
dot_PLOT <- dotplot_OG(sce_lin_subset, cluster_name = cluster_name, fun = fun, scale = FALSE, output_dir = dotplot_raw_dir)
file_name <- paste0(cluster_name, "_dotplot_", fun, ".png")
ggsave(filename = file.path(dotplot_raw_dir, file_name), plot = dot_PLOT, width = 10, height = 8)

# processed, no scale
# dot_PLOT <- dotplot_OG(sce, cluster_name = cluster_name, fun = fun, scale = FALSE, output_dir = dotplot_proc_dir)
# file_name <- paste0(cluster_name, "_dotplot_", fun, ".png")
# ggsave(filename = file.path(dotplot_proc_dir, file_name), plot = dot_PLOT, width = 10, height = 8)



# load RAW & RAW_SUBSET dotplot matrices
# marker x cluster
raw_dp <- read.delim("U:/cdbui/MPE_Cytof/CYTOF_data/Analysis/TBNK/Results/Dotplot/No_Scaled/Raw/meta20_dotplot_NO_SCALED_expr_wide.txt", header = TRUE, check.names = FALSE)
raw_subset_dp <- read.delim(file.path(subset_dir, "Dotplot/No_Scaled/Raw/meta20_subset_dotplot_NO_SCALED_expr_wide.txt"), header = TRUE, check.names = FALSE)

# select only common markers for raw_dp
common_markers <- raw_subset_dp[, 1]
raw_dp <- raw_dp[raw_dp[, 1] %in% common_markers, ]

# check markers are in same order
check <- identical(raw_dp[, 1], raw_subset_dp[, 1])
if (!check) {
  cat("MISMATCH ROWS")
} else {
  cat("OK")
}

# only get numeric columns for correlation calculation
raw_dp <- raw_dp[, -1]
raw_subset_dp <- raw_subset_dp[, -1]


# correlation matrix of RAW & RAW_SUBSET dotplot clusters
# marker x cluster, so we do corr by column
no_scaled_corr_dir <- file.path(dotplot_no_scaled_dir, "Corr_Info")
if (!dir.exists(no_scaled_corr_dir)) dir.create(no_scaled_corr_dir, recursive = TRUE)

# raw cluster x raw_subset cluster
dp_corr_mat <- cor(raw_dp, raw_subset_dp)
# NOTE: cluster_info_dir
write.table(dp_corr_mat, file = file.path(no_scaled_corr_dir, paste0("meta", maxK, "_raw_raw_subset_dotplot_corr.txt")), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
hm <- Heatmap(dp_corr_mat,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_names_side = "left",
              row_title = "Raw Meta20",
              column_names_side = "bottom",
              column_names_rot = 0,
              column_title = "Raw Meta20 Subset",
              column_title_side = "bottom")
png(file.path(no_scaled_corr_dir, paste0("meta", maxK, "_raw_raw_subset_dotplot_corr.png")), width = 3000, height = 3000, res = 300)
draw(hm, column_title = "Raw - Raw_Subset Cluster Correlation (from Dotplots)")
dev.off()


#-------------------------------------------------------------------------------
# ANALYSIS - CHECK %TILE OF VALUES BY MARKER
#-------------------------------------------------------------------------------

# full lineage
es <- assay(sce, "exprs")
marker <- "KIR2DL3" # CD19, KIR2DL3
marker_index <- which(rowData(sce)$Symbol == marker)

# subset lineage
es <- assay(sce_lin_subset, "exprs")
marker <- "KIR2DL3" # CD19, KIR2DL3
marker_index <- which(rowData(sce_lin_subset)$Symbol == marker)


x <- es[marker_index, ]

marker_stats <- c(
  mean = mean(x),
  `50%` = unname(quantile(x, 0.50)),
  `75%` = unname(quantile(x, 0.75)),
  `90%` = unname(quantile(x, 0.90)),
  `95%` = unname(quantile(x, 0.95)),
  max = max(x)
)



# check aggregate stats per cluster for selected markers

# full lineage
es <- assay(sce, "exprs")
marker <- "KIR2DL3" # CD19, KIR2DL3
marker_index <- which(rowData(sce)$Symbol == marker)

x <- es[marker_index, ]
clust_vec <- colData(sce)$meta20

stats <- t( sapply( split(x, clust_vec), function(v) {
  c(mean = mean(v, na.rm = TRUE),
    `25%` = unname(quantile(v, 0.25, na.rm = TRUE)),
    `50%` = unname(quantile(v, 0.50, na.rm = TRUE)),   # same as median
    `75%` = unname(quantile(v, 0.75, na.rm = TRUE)),
    `90%` = unname(quantile(v, 0.90, na.rm = TRUE)),
    `95%` = unname(quantile(v, 0.95, na.rm = TRUE)),
    max   = max(v, na.rm = TRUE))
}) )

stats <- stats[order(as.numeric(rownames(stats))), ]

write.table(stats, file = file.path(res_dir, paste0("full_lineage_cluster_", marker, "_stats.txt")), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)



# subset lineage
es <- assay(sce_lin_subset, "exprs")
marker <- "KIR2DL3" # CD19, KIR2DL3
marker_index <- which(rowData(sce_lin_subset)$Symbol == marker)

x <- es[marker_index, ]
clust_vec <- colData(sce_lin_subset)$meta20_subset

stats <- t( sapply( split(x, clust_vec), function(v) {
  c(mean = mean(v, na.rm = TRUE),
    `25%` = unname(quantile(v, 0.25, na.rm = TRUE)),
    `50%` = unname(quantile(v, 0.50, na.rm = TRUE)),   # same as median
    `75%` = unname(quantile(v, 0.75, na.rm = TRUE)),
    `90%` = unname(quantile(v, 0.90, na.rm = TRUE)),
    `95%` = unname(quantile(v, 0.95, na.rm = TRUE)),
    max   = max(v, na.rm = TRUE))
}) )

stats <- stats[order(as.numeric(rownames(stats))), ]

write.table(stats, file = file.path(res_dir, paste0("subset_lineage_cluster_", marker, "_stats.txt")), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)


# write.table(stats, file = file.path(res_dir, "dummy.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



# marker %tiles
tbnk_marker_info <- calculate_marker_quantiles(sce)
file_name <- file.path(analysis_dir, "TBNK_Marker_Stats_Full_Lineage.txt")
write.table(tbnk_marker_info, file = file_name, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

myeloid_marker_info <- calculate_marker_quantiles(sce)
file_name <- file.path(analysis_dir, "Myeloid_Marker_Stats_Full_Lineage.txt")
write.table(myeloid_marker_info, file = file_name, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

cytokine_marker_info <- calculate_marker_quantiles(sce)
file_name <- file.path(analysis_dir, "Cytokine_Marker_Stats_Full_Lineage.txt")
write.table(cytokine_marker_info, file = file_name, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
