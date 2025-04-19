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
selPanel <- c("TBNK")  #*******
# selPanel <- c("Myeloid")  #*******
# selPanel <- c("Cytokines")  #*******

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

write.table(Ftab_quantile_ref, "tmp.txt", sep = "\t", quote = FALSE, row.names = FALSE)


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
# check panel expression in each reference TBNK cells
# 2/28/25
# --------------------------------------------------------------
selPanel <- c("TBNK")  #*******
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
all_fcs_files <- list.files(bead_norm_dir, 
                            pattern = ".fcs$", 
                            full.names = TRUE)

# setup input folder
bead_norm_dir <- file.path(workFolder,"CYTOF_data",
                           "Gated",selPanel)
all_fcs_files <- list.files(bead_norm_dir, 
                            pattern = ".fcs$", 
                            full.names = TRUE)

# setup input folder
bead_norm_dir <- file.path(workFolder,"CYTOF_data",
                           "Gated_2gates",selPanel)
all_fcs_files <- list.files(bead_norm_dir, 
                            pattern = ".fcs$", 
                            full.names = TRUE)

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
write.table(Ftab_out,file.path("Ranalysis",fout),
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
write.table(Ftab_out,file.path("Ranalysis",fout),
            sep="\t",quote=F,row.names=FALSE)

# ------
fout <- paste(selPanel,"_cleanedGated_arcsinhTransform_qc.txt",sep="")
write.table(Ftab_out,file.path("Ranalysis",fout),
            sep="\t",quote=F,row.names=FALSE)
# -----
fout <- paste(selPanel,"_cleanedGated2gates_arcsinhTransform_qc.txt",sep="")
write.table(Ftab_out,file.path("Ranalysis",fout),
            sep="\t",quote=F,row.names=FALSE)


# -----------------------------------------
#export events number
# -------------------------------------------
selPanel <- c("TBNK")  #*******

#---------------------
# setup input folder
bead_norm_dir <- file.path(workFolder,"CYTOF_data",
                           "BeadNorm",selPanel)
all_fcs_files <- list.files(bead_norm_dir, 
                            pattern = ".fcs$", 
                            full.names = TRUE)

# setup input folder
bead_norm_dir <- file.path(workFolder,"CYTOF_data",
                           "Gated",selPanel)
all_fcs_files <- list.files(bead_norm_dir, 
                            pattern = ".fcs$", 
                            full.names = TRUE)

bead_norm_dir <- file.path(workFolder,"CYTOF_data",
                           "Cleaned",selPanel)
all_fcs_files <- list.files(bead_norm_dir, 
                            pattern = ".fcs$", 
                            full.names = TRUE)

# -----------------------
nfiles <- length(all_fcs_files)
Ftab <- matrix(NA,nfiles,3)

pb <- reloadProgressBar(nfiles)
for (i in c(1:nfiles)){   
  f_raw <- flowCore::read.FCS(all_fcs_files[i],transformation=FALSE)
  orgName <- unlist(keyword(f_raw ,"$FIL"))
  Ftab[i,1] <- basename(all_fcs_files[i])
  Ftab[i,2] <- orgName 
  Ftab[i,3] <- nrow(f_raw)
  pb$tick()
}
colnames(Ftab) <- c("file","FIL","event__no")
write.table(Ftab,file.path("Ranalysis","tmp.txt"),
            sep="\t",quote=F,row.names=FALSE)

# ------------------------------------------------------------

library(ggcyto)

fcs_files <- grep("Ref" , all_fcs_files,value = TRUE)
i <- 1
f <- flowCore::read.FCS(fcs_files[i])
keyword(f, "$FIL")

spillover(f)
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
list_panels <- c("TBNK")


# filter for only REF here


for (selPanel in list_panels){
  
  # get sample info
  sampleInfo  <- dplyr::filter(allSampleInfo,
                               panel_id==selPanel)  #**** & patient_id =="Ref")
  
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
  flowAI_clean_dir <- file.path(getwd(), "CYTOF_data", "Cleaned-flowAI")
  if (!dir.exists(flowAI_clean_dir)) dir.create(flowAI_clean_dir, recursive = TRUE)
  
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
  for (file in files_b) {
    
    # read fcs file
    # ff <- flowCore::read.FCS(filename = file, 
    #                          transformation = FALSE)
    
    
    # clean flow rate: use flowAI using Time
    # using single function flow_auto_qc
    ff <- flowAI::flow_auto_qc(file,
               remove_from = "all",   #***FR,FS, and FM
              output = 1,            #**only clean signal==1
              timeCh = "Time",
              ChExcludeFS = channels_to_exclude,
              html_report=FALSE,
              mini_report=FALSE,
              fcs_QC = FALSE,
              folder_results = FALSE)
    # ff_2@exprs[,"QCvector"]-> ttmp
    
    # ff_3 <- clean_signal(flow_frame = ff,
    #                      channels_to_clean=channels_to_clean,
    #                      to_plot = "None",
    #                      out_dir = clean_dir,
    #                      Segment = 1000,
    #                      arcsine_transform = TRUE,
    #                      data_type = "MC",
    #                      non_used_bead_ch = "140")
    
    # Write FCS files
    # flowCore::write.FCS(ff_3,
    #                     file = file.path(clean_dir, gsub("_beadNorm","_cleaned", basename(file)))) 
    
    # Write FCS files from flowAI
    flowCore::write.FCS(ff,
                        file = file.path(flowAI_clean_dir, gsub("_beadNorm", "_flowAIcleaned", basename(file))))
    
    pb$tick()
  }
}






