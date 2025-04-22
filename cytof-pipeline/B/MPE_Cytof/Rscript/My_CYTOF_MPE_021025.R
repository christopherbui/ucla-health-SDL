# **************************************************************************
# Inventory of files
# 2/10/25
# **************************************************************************
# Run in R v4.3.3
workFolder <- paste("C:/Users/cdbui/Box/ChristopherProjects/", "MPE_Cytof/", sep = "")
setwd(workFolder)

parentDir <- paste("C:/Users/cdbui/Box/ChristopherProjects/", "MPE_Cytof/", sep = "")

tmp_listFiles <- list.files(path = parentDir,
                            pattern = ".fcs",
                            recursive = TRUE,
                            full.names = FALSE)

tmp_fileNames <- sapply(tmp_listFiles,
                        function(x) basename(x),
                        simplify = TRUE)

tmp_out <- cbind(tmp_listFiles, tmp_fileNames)
write.table(tmp_out,
            "tmp.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


# Run in R v4.3.3
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

source(paste("C:/Users/cdbui/Box/ChristopherProjects/MPE_Cytof/Rscript/", "Rybakowska_cytof_function.R", sep = ""))
source(paste("C:/Users/cdbui/Box/ChristopherProjects/MPE_Cytof/Rscript/", "Rybakowska_cytof_function_LT.R", sep = ""))

workFolder <- paste("C:/Users/cdbui/Box/ChristopherProjects/MPE_Cytof", sep = "")
setwd(workFolder)

# progress bar
library(progress)

reloadProgressBar <- function(iterations) {
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = iterations, clear = FALSE, width = 60
  )
  return(pb)
}


# Step 1: setting environment
# Get sample info of all samples/panels
fin_info <- file.path("Ranalysis","mpe_cytof_sampleInfo_022625.txt")
allSampleInfo <- read.delim(fin_info, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
allSampleInfo <- allSampleInfo[, -1]

# ------------------------------------------------------------------------------
# Check quantile of reference samples in different batches 
#-------------------------------------------------------------------------------
# -------
# selPanel <- c("TBNK")  #*******
# selPanel <- c("Myeloid")  #*******
selPanel <- c("Cytokines")  #*******


# make output directories
panel_output_dir <- file.path(workFolder, "Ranalysis", selPanel)
if(!dir.exists(panel_output_dir)) dir.create(panel_output_dir)



#--------
# get sample info
sampleInfo  <- dplyr::filter(allSampleInfo,
                 panel_id == selPanel & patient_id == "Ref")

# get panel info
fin_panel <- paste(selPanel, "_markers_022625.txt", sep = "")
panel_info <- read.delim(file.path("Ranalysis", fin_panel),
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# setup input folder
bead_norm_dir <- file.path(workFolder, "CYTOF_data",
                   "BeadNorm", selPanel)
select_fcs_files <- sampleInfo$file_name    #***must be in the same folder

# Define batch id and sample id for each file
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
# check panel expression in each reference Cytokine/Myeloid/TBNK cells
# 2/28/25
# --------------------------------------------------------------
# selPanel <- c("TBNK")  #*******
# selPanel <- c("Myeloid") #*******
selPanel <- c("Cytokines")


sampleInfo  <- dplyr::filter(allSampleInfo,
                             panel_id == selPanel & patient_id == "Ref")

# get panel info
fin_panel <- paste(selPanel, "_markers_022625.txt", sep = "")
panel_info <- read.delim(file.path("Ranalysis", fin_panel),
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)

## get channel information for fsom
lineage_idx <- which(panel_info$marker_class=="lineage")  #****
lineage_markers <- panel_info$fcs_desc[lineage_idx]  
function_idx <- which(panel_info$marker_class=="function")  #****
function_markers <- panel_info$fcs_desc[function_idx] 

## setup panel input for CATALYST, requiring 3 columns
tmp_panel <- panel_info[,c("fcs_colname","antigen","marker_class")]
tmpi <- which(grepl("function",panel_info$marker_class)==TRUE)
tmp_panel$marker_class[tmpi] <- c("state")
tmpi <- which(panel_info$marker_class=="lineage")
tmp_panel$marker_class[tmpi] <- c("type")

#---------------------
# setup input folder
bead_norm_dir <- file.path(workFolder,"CYTOF_data",
                           "BeadNorm",selPanel)

# GET ONLY REF FILES; JUST FOR TESTING
all_fcs_files <- list.files(bead_norm_dir,
                            pattern = "_Ref.*\\.fcs$",
                            full.names = TRUE)


# all_fcs_files <- list.files(bead_norm_dir, 
#                             pattern = ".fcs$", 
#                             full.names = TRUE)
# 
# # setup input folder
# bead_norm_dir <- file.path(workFolder,"CYTOF_data",
#                            "Gated",selPanel)
# all_fcs_files <- list.files(bead_norm_dir, 
#                             pattern = ".fcs$", 
#                             full.names = TRUE)
# 
# # setup input folder
# bead_norm_dir <- file.path(workFolder,"CYTOF_data",
#                            "Gated_2gates",selPanel)
# all_fcs_files <- list.files(bead_norm_dir, 
#                             pattern = ".fcs$", 
#                             full.names = TRUE)

# ---
nfiles <- length(all_fcs_files)
quantile_values <-  c(0.01, 0.25, 0.5, 0.75, 0.99)
Ftab <- NULL

pb <- reloadProgressBar(nfiles)
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

# -------------------
nfiles <- length(all_fcs_files)
quantile_values <-  c(0.01, 0.25, 0.5, 0.75, 0.99)
Ftab <- NULL

pb <- reloadProgressBar(nfiles)
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
#export events number
# -------------------------------------------
# selPanel <- c("TBNK")  #*******
# selPanel <- c("Myeloid") #******
selPanel <- c("Cytokines")

#---------------------
# setup input folder
bead_norm_dir <- file.path(workFolder,"CYTOF_data",
                           "BeadNorm",selPanel)

# GET ONLY REF FILES; JUST FOR TESTING
# all_fcs_files <- list.files(bead_norm_dir,
#                             pattern = "_Ref.*\\.fcs$",
#                             full.names = TRUE)

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
Ftab <- matrix(NA,nfiles,3)

# normalized, cleaned
all_fcs_files_cleaned <- list.files(clean_dir,
                            pattern = ".fcs$",
                            full.names = TRUE)
nfiles_cleaned <- length(all_fcs_files_cleaned)
Ftab <- matrix(NA,nfiles,3)


# change to 'all_fcs_files' OR 'all_fcs_files_cleaned' below
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
write.table(Ftab,file.path(panel_output_dir, Ftab_file_name),
            sep="\t",quote=F,row.names=FALSE)



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
# Signal Cleaning --------------------------------------------------------------
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
# Files outliers detection -----------------------------------------------------
# output is the score from fsom
#-------------------------------------------------------------------------------

# list_panels <- c("TBNK","Myeloid","Cytokines")

list_panels <- c("TBNK")

batch4extractPatt <- "(?i).*(batch[0-9]*).*.fcs"
dna_ch_4fsom <- c("191Ir","193Ir")   #*****using desc column
via_ch_4fsom <- c("195Pt")    #***Cisplatin
nCells_thres <- c(10000,3000,1000)   #*****


for (k in c(1:length(list_panels))){
  selPanel <- list_panels[k]
  tmp_nCells_thres <- nCells_thres[k]
  
  # get sample info
  # GET ONLY REF FILES; JUST FOR TESTING
  sampleInfo  <- dplyr::filter(allSampleInfo,
                               panel_id==selPanel & patient_id =="Ref")
  
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
                         "Cleaned",selPanel)
  
  # Define out_dir for diagnostic plots
  quality_dir <- file.path(workFolder,"CYTOF_data", 
                           "Quality_control",selPanel)
  if(!dir.exists(quality_dir)) dir.create(quality_dir, recursive = TRUE)
  
  # Define which files will be cleaned
  ###select_fcs_files <- sampleInfo$file_name    #***must be in the same folder
  files_b <- list.files(clean_dir, 
                        pattern = "_cleaned.fcs$", 
                        full.names = TRUE)
  files_b <- files_b[-1]
  
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
  
  
  #--------------------------------------
  # BAD OUTPUT
  #
  # Maybe some markers have NULL descriptions? 
  #
  # Yb176Di
  # Ir191Di
  # Pt195Di
  # Pt198Di
  # Bi209Di
  # Error in seq.default(min(x, na.rm = T), max(x, na.rm = T), length.out = n +  :
  #                        'from' must be a finite number
  #                      In addition: Warning messages:
  #                        1: In min(x, na.rm = T) : no non-missing arguments to min; returning Inf
  #                      2: In max(x, na.rm = T) : no non-missing arguments to max; returning -Inf
  #--------------------------------------
  
}



selPanel <- "TBNK"

batch4extractPatt <- "(?i).*(batch[0-9]*).*.fcs"
dna_ch_4fsom <- c("191Ir","193Ir")   #*****using desc column
via_ch_4fsom <- c("195Pt")    #***Cisplatin
nCells_thres <- c(10000,3000,1000)   #*****

# get panel info --> function and lineage markers
fin_panel <- paste(selPanel,"_markers_022625.txt",sep="")
panel_info <- read.delim(file.path("Ranalysis",fin_panel),
                         sep="\t",header=T,stringsAsFactors=FALSE)

lineage_idx <- which(panel_info$marker_class=="lineage")  #****
lineage_markers <- panel_info$fcs_desc[lineage_idx]  
function_idx <- which(panel_info$marker_class=="function")  #****
function_markers <- panel_info$fcs_desc[function_idx] 
pheno_4fsom <- c(dna_ch_4fsom,via_ch_4fsom,lineage_markers)    #*******


fcs_path <- "C:/Users/cdbui/Box/ChristopherProjects/MPE_Cytof/CYTOF_data/Cleaned/TBNK"
fcs_files <- list.files(fcs_path, pattern = ".fcs$", full.names = TRUE)

problem_list <- list()

# check if markers have NULL description
for (f in fcs_files) {
  ff <- read.FCS(f, transformation = FALSE)
  descs <- pData(parameters(ff))$desc
  names(descs) <- pData(parameters(ff))$name
  
  missing_markers <- pheno_4fsom[pheno_4fsom %in% descs[is.na(descs)]]
  print(missing_markers)
  if (length(missing_markers) > 0) {
    problem_list[[basename(f)]] <- missing_markers
  }
}


# check if expression data has nulls
for (f in fcs_files) {
  ff <- read.FCS(f, transformation = FALSE)
  param_info <- pData(parameters(ff))
  desc_map <- param_info$desc
  names(desc_map) <- param_info$name
  
  cat("Checking:", basename(f), "\n")
  missing_markers <- c()
  
  for (marker in pheno_4fsom) {
    channel <- names(desc_map)[desc_map == marker]
    
    # Skip if no channels match the marker
    if (length(channel) == 0) next
    
    # Filter to only existing columns in exprs
    valid_channels <- channel[channel %in% colnames(exprs(ff))]
    
    # Loop through valid channels
    for (ch in valid_channels) {
      expr_data <- exprs(ff)[, ch]
      if (all(is.na(expr_data))) {
        cat(sprintf("  Marker: %-10s | Channel: %-10s | all NA: TRUE | n events: %d\n",
                    marker, ch, length(expr_data)))
        missing_markers <- c(missing_markers, marker)
      }
    }
  }
  
  
  if (length(missing_markers) > 0) {
    problem_list[[basename(f)]] <- missing_markers
  }
}

# Final summary
cat("\nSummary of files with all-NA markers:\n")
print(problem_list)

print(problem_list)





# get panel info --> function and lineage markers
fin_panel <- paste(selPanel,"_markers_022625.txt",sep="")
panel_info <- read.delim(file.path("Ranalysis",fin_panel),
                         sep="\t",header=T,stringsAsFactors=FALSE)

lineage_idx <- which(panel_info$marker_class=="lineage")  #****
lineage_markers <- panel_info$fcs_desc[lineage_idx]  
function_idx <- which(panel_info$marker_class=="function")  #****
function_markers <- panel_info$fcs_desc[function_idx] 
pheno_4fsom <- c(dna_ch_4fsom,via_ch_4fsom,lineage_markers)    #*******






















ff_test <- read.FCS("C:/Users/cdbui/Box/ChristopherProjects/MPE_Cytof/CYTOF_data/Cleaned/TBNK/TBNK_Ref-PBMC_01_batch1_cleaned.fcs", transformation = FALSE)
param <- flowCore::parameters(ff_test)
desc <- flowCore::pData(param)$desc
print(desc)
param <- pData(parameters(ff_test))
param
