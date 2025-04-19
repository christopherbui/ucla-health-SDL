# **************************************************************************
# Inventory of files
# 2/10/25
# **************************************************************************
# Run in R v4.3.3
workFolder <- paste("D:/","MPE_Cytof/",sep="")
setwd(workFolder)


parentDir <- paste("D:/","MPE_Cytof/",sep="")
tmp_listFiles <- list.files(path = parentDir, pattern=".fcs", 
           recursive = TRUE, full.names = FALSE)
tmp_fileNames <- sapply(tmp_listFiles,function(x) basename(x),simplify=TRUE)

tmp_out <- cbind(tmp_listFiles,tmp_fileNames)
write.table(tmp_out,"tmp.txt",sep="\t",quote=F,row.names=F,col.names=F)


# ***************************************************************************
# Check channel
# 2/24/25
# ***************************************************************************
# Run in R v4.3.3
workFolder <- paste("D:/","MPE_Cytof/Ranalysis/",sep="")
setwd(workFolder)

library(flowCore)

parentInputFiolder <- c("D:/MPE_Cytof/") #****

fin <- c("mpecytof_original_fcsFiles_021025.txt") #****
file_list <- read.delim(fin,sep="\t",header=T,
                 stringsAsFactors=FALSE)

selPanel <- c("CYTOKINE")     #****
selPanel <- c("TBNK")     #****
selPanel <- c("MYELOID")     #****

tmpidx <- which(grepl(selPanel,file_list[,2])==TRUE)
tmpSel_fcs <- file_list[tmpidx,1]
nfiles <- length(tmpidx)

tmp_fin <- paste(parentInputFiolder,tmpSel_fcs[1],sep="")
tmp1 <- read.FCS(tmp_fin, transformation=FALSE)

panel_info <- pData(parameters(tmp1))
Ftab_matched <- matrix(NA,nfiles,3)
Ftab_max <- matrix(0,nrow(panel_info ),nfiles)
Ftab_max[,1] <- panel_info$maxRange
for (i in c(2:nfiles)){
   tmp_fin <- paste(parentInputFiolder,tmpSel_fcs[i],sep="")
   x <- read.FCS(tmp_fin, transformation=FALSE)
   tmp_channel_info <- pData(parameters(x))
   Ftab_max[,i] <- tmp_channel_info$maxRange 
   Ftab_matched[i,1] <- all(panel_info$name==tmp_channel_info$name)
   Ftab_matched[i,2] <- all(panel_info$desc==tmp_channel_info$desc,na.rm=T)
   Ftab_matched[i,3] <- cor(panel_info$maxRange,tmp_channel_info$maxRange,
                             method="spearman")
}
colnames(Ftab_matched) <- c("MatchedName","MatchedDesc","MaxBasedCorr")
rownames(Ftab_matched) <- file_list[tmpidx,2]

rownames(Ftab_max) <- panel_info$name
colnames(Ftab_max) <- file_list[tmpidx,2]


write.table(Ftab_matched,"tmp.txt",sep="\t",quote=F)
write.table(Ftab_max,"tmp2.txt",sep="\t",quote=F)

# ************************************************************************
# Change fcs filename (beadNormlaization)
# 2/26/25
# ************************************************************************
# Run in R v4.3.3
workFolder <- paste("D:/",
       "MPE_Cytof/Ranalysis",sep="")
setwd(workFolder)

inputList <- c("mpe_cytof_nameDictionary_022625.txt")
fileList <- read.delim(inputList,sep="\t",header=T,
              stringsAsFactors=FALSE)

nfiles <- nrow(fileList)
parentInputFolder <- paste("D:/", "MPE_Cytof",sep="")

for (i in c(1:nfiles)){
   old_name <- file.path(parentInputFolder,fileList$originalName[i])
   new_name <- file.path(fileList$new_folder[i],fileList$new_filename[i])

   file.rename(old_name, new_name)
   
}

# **************************************************************************
# Using Rybakowska et al pipeline
# github: https://github.com/prybakowska/CyTOF_analysis_Pipeline1
# Skip
# 1. bead normalization
# Note 
# 1. for channel ID
#     - find_mass_ch(ff_tmp,value = TRUE) --> names of channels(e.g. Cd106Di)
#     - FlowSOM::GetMarkers(.,value=TRUE) --> channel desc(e.g. "106Cd_CD27")
# 2. Use channel 140 as non-bead, non-used channel --> high background signal
# 3. flowAI: https://www.bioconductor.org/packages/devel/bioc/vignettes/flowAI/inst/doc/flowAI.html
#            single function flow_auto_qc whith input from flowFrame
#            (single sample per time) --> loop if from flowSet
#            set html_report=FALSE, mini_report=FALSE to avoid pandoc issue
#    --> skip this step due to memory issue
# 4. flowCut (clean signal) --> myeloid and cytokine data of reference PBMC
#    batch 3,4, and 6 is very small!!!!    
# 5. file_quality_check by fsom 
#    - default cell number for each batch is nsamples*10K --> reduced to 2000
#       if cell number after flow Cut is small (2000 cells)
#   - nClust >=10, error if reduced   
# 6. Normalization using reference sample
#    - select references (i.e. healthy/coltrol) of all batches
#    - using channels is not either Pd/Rh/140 (mpe only use 140
#    - build model applying CytoNorm::QuantileNorm.train on all reference files
#    - normalize using model built from reference files
# 2/26/25
# ***************************************************************************
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

source(paste("D:/PharmG/Rsamples/","Rybakowska_cytof_function.R",sep=""))
source(paste("D:/PharmG/Rsamples/","Rybakowska_cytof_function_LT.R",sep=""))

workFolder <- paste("D:/","MPE_Cytof",sep="")
setwd(workFolder)

# Step 1: setting environment
# Get sample info of all samples/panels
fin_info <- file.path("Ranalysis","mpe_cytof_sampleInfo_022625.txt")
allSampleInfo <- read.delim(fin_info,sep="\t",header=T,stringsAsFactors=FALSE)
allSampleInfo <- allSampleInfo[,-1]

# ------------------------------------------------------------------------------
# Check quantile of reference samples in different batches 
#-------------------------------------------------------------------------------
# -------
selPanel <- c("TBNK")  #*******
selPanel <- c("Myeloid")  #*******
selPanel <- c("Cytokines")  #*******

#--------
# get sample info
sampleInfo  <- dplyr::filter(allSampleInfo,
                 panel_id==selPanel & patient_id =="Ref")
 
# get panel info
fin_panel <- paste(selPanel,"_markers_022625.txt",sep="")
panel_info <- read.delim(file.path("Ranalysis",fin_panel),
      sep="\t",header=T,stringsAsFactors=FALSE)


# setup input folder
bead_norm_dir <- file.path(workFolder,"CYTOF_data",
                   "BeadNorm",selPanel)
select_fcs_files <- sampleInfo$file_name    #***must be in the same folder

# Define batch id and sample id for each file
batch4extractPatt <- "(?i).*(batch[0-9]*).*.FCS"
batch_pattern <- batch4extractPatt

##file_a <- file.path(bead_norm_dir,select_fcs_files)
##batch_pattern <- str_match(basename(select_fcs_files), batch4extractPatt)[,2]

quantiles_ref <- extract_marker_quantiles_4SDL(
      select_fcs_files = select_fcs_files,
      in_dir = bead_norm_dir,
      batch_pattern= batch4extractPatt,
      arcsine_transform = TRUE, 
      markers_to_plot = NULL)


Ftab_quantile_ref <- reshape(quantiles_ref, idvar = c("File","Marker","Batch","Sample"),
                      timevar = "Quantile", direction = "wide")

write.table(Ftab_quantile_ref,"tmp.txt",sep="\t",quote=F,row.names=F)


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
                 panel_id==selPanel & patient_id =="Ref")

# get panel info
fin_panel <- paste(selPanel,"_markers_022625.txt",sep="")
panel_info <- read.delim(file.path("Ranalysis",fin_panel),
      sep="\t",header=T,stringsAsFactors=FALSE)
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
for (i in c(1:nfiles)){   
   f_raw <- flowCore::read.FCS(all_fcs_files[i],transformation=FALSE)
   orgName <- unlist(keyword(f_raw ,"$FIL"))
   Ftab[i,1] <- basename(all_fcs_files[i])
   Ftab[i,2] <- orgName 
   Ftab[i,3] <- nrow(f_raw)
}
colnames(Ftab) <- c("file","FIL","event__no")
write.table(Ftab,file.path("Ranalysis","tmp.txt"),
   sep="\t",quote=F,row.names=FALSE)

# ------------------------------------------------------------

library(ggcyto)

##fcs_files <- grep("Ref",all_fcs_files,value=TRUE)
i <- 1
f <- flowCore::read.FCS(fcs_files[i])
keyword(f,"$FIL")

spillover(f)

tmp_channels <- colnames(f)
channels_info <- pData(parameters(f))

sel_idx <- which(grepl("CD45",channels_info$desc)==TRUE)
sel_channel <- channels_info$name[sel_idx]

 autoplot(f, "Pt198Di")


# ------------------------------------------------------------------------------
# Signal Cleaning --------------------------------------------------------------
#-------------------------------------------------------------------------------
list_panels <- c("TBNK","Myeloid","Cytokines")

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
   if(!dir.exists(clean_dir)) dir.create(clean_dir)

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
      ff <- flowCore::read.FCS(filename = file, 
                           transformation = FALSE)
      
  
     # clean flow rate: use flowAI using Time
     # using single function flow_auto_qc
     ##ff <- flowAI::flow_auto_qc(ff,
     ##            remove_from = "all",   #***FR,FS, and FM
     ##           output = 1,            #**only clean signal==1
     ##           timeCh = "Time",
     ##           ChExcludeFS = channels_to_exclude,
     ##           html_report=FALSE,
     ##           mini_report=FALSE,
     ##           fcs_QC = FALSE, 
     ##           folder_results = FALSE)
     ###ff_2@exprs[,"QCvector"]-> ttmp

     ff_3 <- clean_signal(flow_frame = ff,
                     channels_to_clean=channels_to_clean,
                     to_plot = "None",
                     out_dir = clean_dir,
                     Segment = 1000,
                     arcsine_transform = TRUE,
                     data_type = "MC",
                     non_used_bead_ch = "140")

     # Write FCS files
     flowCore::write.FCS(ff_3,
            file = file.path(clean_dir, gsub("_beadNorm","_cleaned", basename(file)))) 

   }
}

# ------------------------------------------------------------------------------
# Files outliers detection -----------------------------------------------------
# output is the score from fsom
#-------------------------------------------------------------------------------
list_panels <- c("TBNK","Myeloid","Cytokines")
batch4extractPatt <- "(?i).*(batch[0-9]*).*.fcs"
dna_ch_4fsom <- c("191Ir","193Ir")   #*****using desc column
via_ch_4fsom <- c("195Pt")    #***Cisplatin
nCells_thres <- c(10000,3000,1000)   #*****

for (k in c(1:length(list_panels))){
   selPanel <- list_panels[k]
   tmp_nCells_thres <- nCells_thres[k]

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
                   "Cleaned",selPanel)

   # Define out_dir for diagnostic plots
   quality_dir <- file.path(workFolder,"CYTOF_data", 
                  "Quality_control",selPanel)

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

# ------------------------
k <- 3   #errors at batch 4 and later
   selPanel <- list_panels[k]
   tmp_nCells_thres <- nCells_thres[k]

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
                   "Cleaned",selPanel)

   # Define out_dir for diagnostic plots
   quality_dir <- file.path(workFolder,"CYTOF_data", 
                  "Quality_control",selPanel)

   # Define which files will be cleaned
   ###select_fcs_files <- sampleInfo$file_name    #***must be in the same folder
   files_b <- list.files(clean_dir, 
                    pattern = "_cleaned.fcs$", 
                    full.names = TRUE)

   # Define batch_id for each file 
   file_batch_id <- stringr::str_match(basename(files_b),batch4extractPatt)[,2] #*****
        
   fcs_files = files_b
   file_batch_id = file_batch_id 
   nCellsPerSample = tmp_nCells_thres
   out_dir = quality_dir
   phenotyping_markers = pheno_4fsom
   arcsine_transform = TRUE
   nClus = 10
   sd = 3

  
  if (!is.null(file_batch_id)) {
    scores <- list()
    for (batch in unique(file_batch_id)){
      print(batch)
      
      tmp_files <- fcs_files[file_batch_id == batch]
      nCells <- nCellsPerSample*length(tmp_files)
      fsom <- fsom_aof_4SDL(fcs_files = tmp_files, 
                       phenotyping_markers = phenotyping_markers, 
                       nCells = nCells,
                       out_dir = out_dir, 
                       arcsine_transform = arcsine_transform,
                       nClus = nClus,
                       batch = batch)
    
      scores[[batch]] <- aof_scoring_4SDL(fcs_files = tmp_files, 
                                     phenotyping_markers = phenotyping_markers,
                                     fsom = fsom, out_dir = out_dir, batch = batch)
    }
    
  } else {
    tmp_files <- fcs_files
    fsom <- fsom_aof_4SDL(fcs_files = tmp_files, phenotyping_markers = phenotyping_markers, 
                     out_dir = out_dir, arcsine_transform = arcsine_transform, 
                     nClus = nClus,
                     batch = NULL)
    
    scores <- aof_scoring_4SDL(fcs_files = tmp_files, 
                          phenotyping_markers = phenotyping_markers,
                          fsom = fsom, out_dir = out_dir, batch = NULL)
  }
  
  final_score <- file_outlier_detecion(scores = scores, out_dir = out_dir, 
                                       sd = sd)

# ------------------------------------------------------------------------------
# Files gating -----------------------------------------------------------------
# Three gates: intact, doublet and live 
# --> renames to Gated_3gated
#-------------------------------------------------------------------------------
selPanel <- c("TBNK")  #*******

# Set input directory - clean_dir
aggregate_dir <- file.path(workFolder,"CYTOF_data", 
                   "Cleaned",selPanel)

# Set output directory 
gate_dir <- file.path(workFolder,"CYTOF_data", 
                  "Gated",selPanel)
if (!dir.exists(gate_dir)) { 
  dir.create(gate_dir)
}

# List files for gating 
files <- list.files(path = aggregate_dir, 
                    pattern = ".fcs$", 
                    full.names = TRUE)

# Gate the files and plot the gating strategy for each file 
n_plots <- 6  
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


# ----------------------------------------------------------------------------------------
# update '$FIL' in fcs files --> output in the subfolder, "updateFileName"
# Set input directory 
selPanel <- c("TBNK")  #*******
input_dir <-  file.path(workFolder,"CYTOF_data", 
                  "Gated",selPanel)

fcs_files <- list.files(input_dir, 
                    pattern = ".fcs$", 
                    full.names = TRUE)

change_fcs_FIL(fcs_files)

# ------------------------------------------------------------------------------
# Files gating -----------------------------------------------------------------
# ****Two gates: doublet and live 
# --> renames to Gated_2gated
#-------------------------------------------------------------------------------
selPanel <- c("TBNK")  #*******

# Set input directory - clean_dir
aggregate_dir <- file.path(workFolder,"CYTOF_data", 
                   "Cleaned",selPanel)

# Set output directory 
gate_dir <- file.path(workFolder,"CYTOF_data", 
                  "Gated",selPanel)
if (!dir.exists(gate_dir)) { 
  dir.create(gate_dir)
}

# List files for gating 
files <- list.files(path = aggregate_dir, 
                    pattern = ".fcs$", 
                    full.names = TRUE)

# Gate the files and plot the gating strategy for each file 
n_plots <- 6  
png(file.path(gate_dir, paste0("gating.png")),
    width = n_plots * 300, 
    height = length(files) * 300)
layout(matrix(1:(length(files) * n_plots), 
              ncol = n_plots, 
              byrow = TRUE))

for (file in files){
  
  ff <- flowCore::read.FCS(filename = file, 
                           transformation = FALSE)
  
  ##ff <- gate_intact_cells(flow_frame = ff, 
  ##                        file_name = basename(file))
  
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

# ----------------------------------------------------------------------------------------
# update '$FIL' in fcs files --> output in the subfolder, "updateFileName"
# Set input directory 
selPanel <- c("TBNK")  #*******
input_dir <-  file.path(workFolder,"CYTOF_data", 
                  "Gated",selPanel)

fcs_files <- list.files(input_dir, 
                    pattern = ".fcs$", 
                    full.names = TRUE)

change_fcs_FIL(fcs_files)

# ------------------------------------------------------------------------------
# Files gating using new strategies:--> *_4SDL functions
# for intact: single low threshold for Ir191 and Ir193 (i.e. >thres)
# for live: only use viability channel(cisplatin - i.e. < thres)
# ****Three gates: intact (Ir193 and Ir191), doublet (Time Event) and live (Cisplatin)
# --> renames to Gated_2gated
#-------------------------------------------------------------------------------
selPanel <- c("TBNK")  #*******
hard_cutoff = 4   #**hard_cutoff for Ir channels in gating intact cells

# Set input directory - clean_dir
aggregate_dir <- file.path(workFolder,"CYTOF_data", 
                   "Cleaned",selPanel)

# Set output directory 
gate_dir <- file.path(workFolder,"CYTOF_data", 
                  "Gated",selPanel)
if (!dir.exists(gate_dir)) { 
  dir.create(gate_dir)
}

# List files for gating 
files <- list.files(path = aggregate_dir, 
                    pattern = ".fcs$", 
                    full.names = TRUE)

# Gate the files and plot the gating strategy for each file 
n_plots <- 3  #intact,singlet,live 
png(file.path(gate_dir, paste0("gating.png")),
    width = n_plots * 300, 
    height = length(files) * 300)
layout(matrix(1:(length(files) * n_plots), 
              ncol = n_plots, 
              byrow = TRUE))

for (file in files){
  
  ff <- flowCore::read.FCS(filename = file, 
                           transformation = FALSE)
  
  ff <- gate_intact_cells_4SDL(flow_frame = ff, 
                          file_name = basename(file),
                          hard_cutoff = hard_cutoff)
  
  ff <- gate_singlet_cells(flow_frame = ff,
                           channels = "Event_length",
                           file_name = basename(file))
  
  ff <- gate_live_cells_4SDL(flow_frame = ff, 
                        viability_channel = "Pt195Di",
                        out_dir = gate_dir)
  
  flowCore::write.FCS(ff, file.path(gate_dir,
                                    gsub(".fcs", "_gated.fcs", basename(file))))
}

dev.off()


# ------------------------------------------------------------------------------
# Normalization using reference sample -----------------------------------------
#-------------------------------------------------------------------------------
selPanel <- c("TBNK")  #*******
sel4ref <- c(1:7)     #***
batch4extractPatt <- "(?i).*(batch[0-9]*).*.fcs"  #****
refPattern <- c("Ref.*_gated.fcs$")
dna_ch <- c("Ir191Di","Ir193Di")   #*****
via_ch <- c("Pt195Di")    #***Cisplatin

# get panel info
fin_panel <- paste(selPanel,"_markers_022625.txt",sep="")
panel_info <- read.delim(file.path("Ranalysis",fin_panel),
     sep="\t",header=T,stringsAsFactors=FALSE)
 
target_ch <- panel_info$fcs_colname
channels_to_norm <- c(dna_ch,via_ch,target_ch)

# Set input directory - clean_dir
gate_dir <- file.path(workFolder,"CYTOF_data", 
                  "Gated",selPanel)

# Define reference samples
files_ref <- list.files(gate_dir, 
                        pattern = refPattern , 
                        full.names = TRUE, 
                        recursive = T)
# ---
files_ref <- files_ref[sel4ref]   #****only for TBNK


# Define batch labels for each files
labels_ref <- stringr::str_match(basename(files_ref), batch4extractPatt)[,2]
                                 
# Define channels to be normalized
ff <- read.FCS(files_ref[1])
##channels <- grep("Pd|Rh|140",grep("Di", colnames(ff),value = T),value = T, invert = T) 
channels <- colnames(ff)[colnames(ff) %in% channels_to_norm]

# Define out_dir for normalized files
norm_dir <- file.path(workFolder,"CYTOF_data", 
                  "CytoNormed",selPanel)
if(!dir.exists(norm_dir))(dir.create(norm_dir))

# Build the normalization model using reference samples and plot quantiles 
png(file.path(norm_dir, "005_095_normalization.png"),
    width = length(channels) * 300,
    height = (length(files_ref) * 2 + 1) * 300)
model <- CytoNorm::QuantileNorm.train(files = files_ref,
                labels = labels_ref, 
                channels = channels, 
                transformList = transformList(channels, CytoNorm::cytofTransform), 
                nQ = 2, 
                limit = c(0,8),
                quantileValues = c(0.05, 0.95), 
                goal = "mean", 
                plot = TRUE)
dev.off()

# save the model
saveRDS(object = model, 
        file = file.path(norm_dir, "005_095_model.RDS"))

# Define path to the files for normalization
files <- list.files(file.path(gate_dir), 
                    pattern = "_gated.fcs$", 
                    full.names = TRUE, recursive = T)

# Define batch labels for each files, note that they need to corresponds to 
# reference labels 
labels <- stringr::str_match(basename(files),batch4extractPatt)[,2]

# Normalize files 
CytoNorm::QuantileNorm.normalize(model = model, 
                  files = files, 
                  labels = labels, 
                  transformList = transformList(channels,CytoNorm::cytofTransform),
                  transformList.reverse = transformList(channels,CytoNorm::cytofTransform.reverse), 
                  outputDir = norm_dir)

# ----------------------------------------------------------------------------------------
# update '$FIL' in fcs files --> output in the subfolder, "updateFileName"
# Set input directory 
norm_dir <-  file.path(workFolder,"CYTOF_data", 
                  "CytoNormed",selPanel)

fcs_files <- list.files(norm_dir, 
                    pattern = ".fcs$", 
                    full.names = TRUE)

change_fcs_FIL(fcs_files)

# *****************************************************************************************
# Run UMAP by CATALYST 
# https://www.bioconductor.org/packages/release/bioc/vignettes/CATALYST/inst/doc/differential.html
# Note: '$FIL" in the input fcs must be matched with $file coumn in the meta-data
#        check $FIL --> change the field value using change_fcs_FIL
#           f <- flowCore::read.FCS(fcs_files[i])
#           keyword(f,"$FIL")
# 2/27/25
#-------------------------------------------------------------------------------
# Run in R v4.3.3
library(stringr)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(pals)
library(scales)
library(ggplot2)
library(pheatmap)
#library(pandoc)   #required by flowAI
library(cowplot)

library(flowCore)
library(flowAI)
library(flowCut)
library(cytutils)   #clauclting aof/change fcs name


library(flowDensity)
library(CytoNorm)
library(CATALYST)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(diffcyt)
library(scater)

source(paste("D:/PharmG/Rsamples/","Rybakowska_cytof_function.R",sep=""))
source(paste("D:/PharmG/Rsamples/","Rybakowska_cytof_function_LT.R",sep=""))

workFolder <- paste("D:/","MPE_Cytof",sep="")
setwd(workFolder)

# Step 1: import meta data
# Get sample info of all samples/panels
fin_info <- file.path("Ranalysis","mpe_cytof_sampleInfo_022625.txt")
allSampleInfo <- read.delim(fin_info,sep="\t",header=T,stringsAsFactors=FALSE)
allSampleInfo <- allSampleInfo[,-1]


# Step 2: import data for the selected panel
# Step 2.1: setting parameters 
selPanel <- c("TBNK")  #*******
sel4ref <- c(1:7)     #***
batch4extractPatt <- ".*(batch[0-9]*).*.fcs"  #****
refPattern <- c("Ref.*_gated.fcs$")
dna_ch <- c("Ir191Di","Ir193Di")   #*****
via_ch <- c("Pt195Di")    #***Cisplatin

# Step 2.2. get panel info
fin_panel <- paste(selPanel,"_markers_022625.txt",sep="")
panel_info <- read.delim(file.path("Ranalysis",fin_panel),
     sep="\t",header=T,stringsAsFactors=FALSE)
 
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

# Step 2.3.get fcs data and their metadata
## Set input directory of fcs
norm_dir <-  file.path(workFolder,"CYTOF_data", 
                  "CytoNormed",selPanel,"updateFileName")   #******

## Set output directory 
analysis_dir <- file.path(workFolder,"CYTOF_data", 
                   "Analysis",selPanel)  

if (!dir.exists(analysis_dir)) { 
  dir.create(analysis_dir)
}

fcs_files <- list.files(norm_dir, 
                    pattern = ".fcs$", 
                    full.names = TRUE)
batch_pattern <- c("batch[0-9]*")

## extract sample information from fcs files
getCoreID <- function(x){
   xlist <- unlist(strsplit(x,split="_"))[2:5]
   xcore <- paste0(xlist,collapse="_")
}

tmp_core <- sapply(basename(fcs_files),getCoreID,simplify=TRUE)

## match filenames and their metadata, then prepare metadata for CATALYST, 
## requiring 2 columns
sampleInfo  <- dplyr::filter(allSampleInfo,
                 panel_id==selPanel)
md_info <- left_join(data.frame(file=basename(fcs_files),sample_id=tmp_core),
                     sampleInfo[,-1],
                     by = join_by(sample_id == corename))
md_info$BATCH <- factor(md_info$BATCH)

## import fcs as sce object: (arcsinh-) transform=TRUE by default
sce <- CATALYST::prepData(fcs_files,
           panel = tmp_panel,
           md = md_info,
           panel_cols = list(channel = "fcs_colname", antigen = "antigen", class = "marker_class"),
           md_cols = list(file = "file", id = "sample_id", 
               factors = c("tissue_type","patient_id","parental_sample_id","BATCH")))


SummarizedExperiment(sce)
colnames(rowData(sce))   #excess features/markers
dim(colData(sce))        #summary samples data

tmp <- table(colData(sce)$sample_id)    #numbers of events/cells per file
tmp_2 <- data.frame(tmp)
colnames(tmp_2)[1] <- c("sample_id")

tmp_2 <- left_join(tmp_2,md_info,by=join_by(sample_id))
kruskal.test(tmp_2$Freq~factor(tmp_2$BATCH))

tmp_2$BATCH <- factor(tmp_2$BATCH)
ggplot(tmp_2, aes(x=BATCH,y=Freq))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(color=tissue_type),width = 0.25)+
  scale_y_continuous(trans='log10')+ theme_bw()

write.table(tmp,"tmp.txt",sep="\t",quote=F,col.names=F)

# Step 3. Inspect the data
# Step 3.1. inspect Ref-PBMC samples
ref_sam <- grep("Ref",md_info$sample_id,value=T)
sce_ref <- sce[,sce$sample_id %in% ref_sam]
colData(sce_ref)$sample_id <- droplevels(colData(sce_ref)$sample_id)

# plot distribution of markers
p_ref <- plotExprs(sce, color_by = "BATCH")

sce_ref <- sce[,sce$sample_id %in% ref_sam[8]]
plotExprs(sce_ref, color_by = "BATCH")

# plot multi-dimension scale (median marker intensities
pbMDS(sce_ref, color_by = "BATCH", label_by = "sample_id")

# heatmap of median marker intensities
# in CyTOF applications, a cosine distance shows good performance
# scale= "last": aggregate then scale & trim
plotExprHeatmap(sce_ref, scale = "last",
    row_anno = FALSE,
    col_anno = FALSE)

# identification with FlowSOM: using k = 20 metacluster
# id of som(i.e. 100 clusters) and metacluster are stored in cluster_ids
set.seed(1234)
sce_ref <- cluster(sce_ref, features = "type",
    xdim = 10, ydim = 10, maxK = 20, seed = 1234)
plotExprHeatmap(sce_ref, features = "type", 
    by = "cluster_id", k = "meta20", 
    bars = TRUE, perc = TRUE)

names(cluster_codes(sce_ref))   #codes of clusters

##access specific clustering resolution, i.e. cellIDs
table(cluster_ids(sce_ref, "som100"))
table(cluster_ids(sce_ref, "meta20"))

# run t-SNE/UMAP on at most 500/1000 cells per sample
set.seed(1234)
sce_ref <- runDR(sce_ref, "UMAP", cells = 1e3, features = "type")
plotDR(sce_ref, "UMAP", color_by = "CD4")

# compare FSOM and UMAP clusters
p1 <- plotDR(sce_ref, "UMAP", color_by = "meta20")
lgd <- get_legend(p2)
p2 <- p2 + theme(legend.position = "none")


p <- plotExprs(sce, color_by = "tissue_type")


p <- plotDR(sce_ref, "UMAP",
  color_by = "meta15", 
  facet_by = "sample_id")

# Step 4: Select samples from batch 2+3+5
sel_batches <- c(2,3,5)
sel_sam <- md_info$sample_id[md_info$BATCH %in% sel_batches]

sce_sel <- sce[,sce$sample_id %in% sel_sam]
colData(sce_sel)$sample_id <- droplevels(colData(sce_sel)$sample_id)


# run FSOM and t-SNE/UMAP on at most 500/1000 cells per sample
set.seed(1234)
tmp_lineage=c("CD4","CD8A","CD45","CD45RA","CD45RO","CD19","CD3","CD56","CD94")
table(rowData(sce_sel)$marker_name %in% tmp_lineage)
length(tmp_lineage)


sce_sel <- cluster(sce_sel, features = "type",
    xdim = 10, ydim = 10, maxK = 20, seed = 1234)

sce_sel <- runDR(sce_sel, "UMAP", cells = 1e3, features = "type")

# plot FSOM cluster
plotExprHeatmap(sce_sel, features = "type", 
    by = "cluster_id", k = "meta20", 
    bars = TRUE, perc = TRUE)

sce_sel_2 <- cluster(sce_sel, features = tmp_lineage,
    xdim = 10, ydim = 10, maxK = 20, seed = 1234)
sce_sel_2 <- runDR(sce_sel_2, "UMAP", cells = 1e3, features = tmp_lineage)

# using "meta8"
plotExprHeatmap(sce_sel_2, features = tmp_lineage, 
    by = "cluster_id", k = "meta8", 
    bars = TRUE, perc = TRUE)

plotDR(sce_sel_2, "UMAP", color_by = "meta8")

u <- filterSCE(sce_sel_2, k = "meta8",
    cluster_id %in% c(1, 7, 6))
plot_grid(
    plotDR(sce_sel_2, color_by = "meta8"),
    plotDR(u, color_by = "meta8"))
       
plotDR(sce_sel_2, "UMAP", color_by = "meta8",
    facet_by="tissue_type")

plotDR(sce_sel_2, "UMAP", color_by = "FOXP3",
    facet_by="tissue_type")

table(cluster_ids(sce_sel_2, "meta8"))

plotExprHeatmap(sce_sel, features = "state", 
    scale = "last", q = 0, bars = FALSE)

plotPbExprs(sce_sel_2, k = "meta8", facet_by = "cluster_id", 
   color_by = "tissue_type",ncol = 4)

plotPbExprs(sce_sel_2, k = "meta8", facet_by = "cluster_id", 
   features="type",
   color_by = "tissue_type",ncol = 4)

u <- filterSCE(sce_sel_2, k = "meta8",    
    cluster_id %in% c(2,3,5))
plotPbExprs(u, k = "meta8", facet_by = "cluster_id",
   features = c("GRANZYME_B","PERFORIN"), 
   color_by = "tissue_type",ncol = 1)



u <- filterSCE(sce_sel_2, k = "meta8",    
    cluster_id %in% c(2,4,5))
plotPbExprs(u, k = "meta8", facet_by = "cluster_id",
   features = c("CD39","CD28"), 
   color_by = "tissue_type",ncol = 1)


u <- filterSCE(sce_sel_2, k = "meta8",    
    cluster_id %in% c(2))
# Step 5 diffcyt
design <- createDesignMatrix(ei(sce_sel_2), cols_design = "tissue_type")
contrast <- createContrast(c(0, 1))

# test for
# - differential abundance (DA) of clusters
# - differential states (DS) within clusters
res_DA <- diffcyt(sce_sel_2, clustering_to_use = "meta8",
    analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
    design = design, contrast = contrast, verbose = FALSE)
res_DS <- diffcyt(sce_sel_2, clustering_to_use = "meta8",
    analysis_type = "DS", method_DS = "diffcyt-DS-limma",
    design = design, contrast = contrast, verbose = FALSE)

# extract result tables
tbl_DA <- rowData(res_DA$res)
tbl_DS <- rowData(res_DS$res)






# *****************************************************************************************
# Run UMAP by CATALYST 
# ******Use gated data*******
#-------------------------------------------------------------------------------
# Run in R v4.3.3
library(stringr)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(pals)
library(scales)
library(ggplot2)
library(pheatmap)
#library(pandoc)   #required by flowAI
library(cowplot)

library(flowCore)
library(flowAI)
library(flowCut)
library(cytutils)   #clauclting aof/change fcs name


library(flowDensity)
library(CytoNorm)
library(CATALYST)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(diffcyt)
library(scater)

source(paste("D:/PharmG/Rsamples/","Rybakowska_cytof_function.R",sep=""))
source(paste("D:/PharmG/Rsamples/","Rybakowska_cytof_function_LT.R",sep=""))

workFolder <- paste("D:/","MPE_Cytof",sep="")
setwd(workFolder)

# Step 1: import meta data
# Get sample info of all samples/panels
fin_info <- file.path("Ranalysis","mpe_cytof_sampleInfo_022625.txt")
allSampleInfo <- read.delim(fin_info,sep="\t",header=T,stringsAsFactors=FALSE)
allSampleInfo <- allSampleInfo[,-1]


# Step 2: import data for the selected panel
# Step 2.1: setting parameters 
selPanel <- c("TBNK")  #*******
sel4ref <- c(1:7)     #***
batch4extractPatt <- ".*(batch[0-9]*).*.fcs"  #****
refPattern <- c("Ref.*_gated.fcs$")
dna_ch <- c("Ir191Di","Ir193Di")   #*****
via_ch <- c("Pt195Di")    #***Cisplatin

# Step 2.2. get panel info
fin_panel <- paste(selPanel,"_markers_022625.txt",sep="")
panel_info <- read.delim(file.path("Ranalysis",fin_panel),
     sep="\t",header=T,stringsAsFactors=FALSE)
 
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

# Step 2.3.get fcs data and their metadata
## Set input directory of fcs
## ***************
norm_dir <-  file.path(workFolder,"CYTOF_data", 
                  "Gated",selPanel,"updateFileName")   #******

## extract sample information from fcs files
getCoreID <- function(x){
   xlist <- unlist(strsplit(x,split="_"))[1:4]    #*****
   xcore <- paste0(xlist,collapse="_")
}
# *************************************************

## Set output directory 
analysis_dir <- file.path(workFolder,"CYTOF_data", 
                   "Analysis",selPanel)  

if (!dir.exists(analysis_dir)) { 
  dir.create(analysis_dir)
}

fcs_files <- list.files(norm_dir, 
                    pattern = ".fcs$", 
                    full.names = TRUE)
batch_pattern <- c("batch[0-9]*")

tmp_core <- sapply(basename(fcs_files),getCoreID,simplify=TRUE)

## match filenames and their metadata, then prepare metadata for CATALYST, 
## requiring 2 columns
sampleInfo  <- dplyr::filter(allSampleInfo,
                 panel_id==selPanel)
md_info <- left_join(data.frame(file=basename(fcs_files),sample_id=tmp_core),
                     sampleInfo[,-1],
                     by = join_by(sample_id == corename))
md_info$BATCH <- factor(md_info$BATCH)

## import fcs as sce object: (arcsinh-) transform=TRUE by default
sce <- CATALYST::prepData(fcs_files,
           panel = tmp_panel,
           md = md_info,
           panel_cols = list(channel = "fcs_colname", antigen = "antigen", class = "marker_class"),
           md_cols = list(file = "file", id = "sample_id", 
               factors = c("tissue_type","patient_id","parental_sample_id","BATCH")))


SummarizedExperiment(sce)
colnames(rowData(sce))   #excess features/markers
dim(colData(sce))        #summary samples data

tmp <- table(colData(sce)$sample_id)    #numbers of events/cells per file
tmp_2 <- data.frame(tmp)
colnames(tmp_2)[1] <- c("sample_id")

tmp_2 <- left_join(tmp_2,md_info,by=join_by(sample_id))
kruskal.test(tmp_2$Freq~factor(tmp_2$BATCH))

tmp_2$BATCH <- factor(tmp_2$BATCH)
ggplot(tmp_2, aes(x=BATCH,y=Freq))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(color=tissue_type),width = 0.25)+
  scale_y_continuous(trans='log10')+ theme_bw()

write.table(tmp,"tmp.txt",sep="\t",quote=F,col.names=F)

# Step 3. Inspect the data
# Step 3.1. inspect Ref-PBMC samples
ref_sam <- grep("Ref",md_info$sample_id,value=T)
sce_ref <- sce[,sce$sample_id %in% ref_sam]
colData(sce_ref)$sample_id <- droplevels(colData(sce_ref)$sample_id)

# plot distribution of markers
p_ref <- plotExprs(sce_ref, color_by = "BATCH")

sce_ref <- sce[,sce$sample_id %in% ref_sam[2]]
plotExprs(sce_ref, color_by = "BATCH")

# plot multi-dimension scale (median marker intensities
pbMDS(sce_ref, color_by = "BATCH", label_by = "sample_id")

# heatmap of median marker intensities
# in CyTOF applications, a cosine distance shows good performance
# scale= "last": aggregate then scale & trim
plotExprHeatmap(sce_ref, scale = "last",
    row_anno = FALSE,
    col_anno = FALSE)

# identification with FlowSOM: using k = 20 metacluster
# id of som(i.e. 100 clusters) and metacluster are stored in cluster_ids
set.seed(1234)
sce_ref <- cluster(sce_ref, features = "type",
    xdim = 10, ydim = 10, maxK = 20, seed = 1234)
plotExprHeatmap(sce_ref, features = "type", 
    by = "cluster_id", k = "meta20", 
    bars = TRUE, perc = TRUE)

names(cluster_codes(sce_ref))   #codes of clusters

##access specific clustering resolution, i.e. cellIDs
table(cluster_ids(sce_ref, "som100"))
table(cluster_ids(sce_ref, "meta20"))

# run t-SNE/UMAP on at most 500/1000 cells per sample
set.seed(1234)
sce_ref <- runDR(sce_ref, "UMAP", cells = 1e3, features = "type")
plotDR(sce_ref, "UMAP", color_by = "CD4")

# compare FSOM and UMAP clusters
p1 <- plotDR(sce_ref, "UMAP", color_by = "meta20")
lgd <- get_legend(p2)
p2 <- p2 + theme(legend.position = "none")


p <- plotExprs(sce, color_by = "tissue_type")


p <- plotDR(sce_ref, "UMAP",
  color_by = "meta15", 
  facet_by = "sample_id")


# *****************************************************************************************
# Run UMAP by CATALYST 
# ******Use gated data*******
#-------------------------------------------------------------------------------
# Run in R v4.3.3
library(stringr)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(pals)
library(scales)
library(ggplot2)
library(pheatmap)
#library(pandoc)   #required by flowAI
library(cowplot)

library(flowCore)
library(flowAI)
library(flowCut)
library(cytutils)   #clauclting aof/change fcs name


library(flowDensity)
library(CytoNorm)
library(CATALYST)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(diffcyt)
library(scater)

source(paste("D:/PharmG/Rsamples/","Rybakowska_cytof_function.R",sep=""))
source(paste("D:/PharmG/Rsamples/","Rybakowska_cytof_function_LT.R",sep=""))

workFolder <- paste("D:/","MPE_Cytof",sep="")
setwd(workFolder)

# Step 1: import meta data
# Get sample info of all samples/panels
fin_info <- file.path("Ranalysis","mpe_cytof_sampleInfo_022625.txt")
allSampleInfo <- read.delim(fin_info,sep="\t",header=T,stringsAsFactors=FALSE)
allSampleInfo <- allSampleInfo[,-1]


# Step 2: import data for the selected panel
# Step 2.1: setting parameters 
selPanel <- c("TBNK")  #*******
sel4ref <- c(1:7)     #***
batch4extractPatt <- ".*(batch[0-9]*).*.fcs"  #****
refPattern <- c("Ref.*_gated.fcs$")
dna_ch <- c("Ir191Di","Ir193Di")   #*****
via_ch <- c("Pt195Di")    #***Cisplatin

# Step 2.2. get panel info
fin_panel <- paste(selPanel,"_markers_022625.txt",sep="")
panel_info <- read.delim(file.path("Ranalysis",fin_panel),
     sep="\t",header=T,stringsAsFactors=FALSE)
 
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

# Step 2.3.get fcs data and their metadata
## Set input directory of fcs
## ***************
norm_dir <-  file.path(workFolder,"CYTOF_data", 
                  "Gated",selPanel,"updateFileName")   #******

## extract sample information from fcs files
getCoreID <- function(x){
   xlist <- unlist(strsplit(x,split="_"))[1:4]    #*****
   xcore <- paste0(xlist,collapse="_")
}
# *************************************************

## Set output directory 
analysis_dir <- file.path(workFolder,"CYTOF_data", 
                   "Analysis",selPanel)  

if (!dir.exists(analysis_dir)) { 
  dir.create(analysis_dir)
}

fcs_files <- list.files(norm_dir, 
                    pattern = ".fcs$", 
                    full.names = TRUE)
batch_pattern <- c("batch[0-9]*")

tmp_core <- sapply(basename(fcs_files),getCoreID,simplify=TRUE)

## match filenames and their metadata, then prepare metadata for CATALYST, 
## requiring 2 columns
sampleInfo  <- dplyr::filter(allSampleInfo,
                 panel_id==selPanel)
md_info <- left_join(data.frame(file=basename(fcs_files),sample_id=tmp_core),
                     sampleInfo[,-1],
                     by = join_by(sample_id == corename))
md_info$BATCH <- factor(md_info$BATCH)

## import fcs as sce object: (arcsinh-) transform=TRUE by default
sce <- CATALYST::prepData(fcs_files,
           panel = tmp_panel,
           md = md_info,
           panel_cols = list(channel = "fcs_colname", antigen = "antigen", class = "marker_class"),
           md_cols = list(file = "file", id = "sample_id", 
               factors = c("tissue_type","patient_id","parental_sample_id","BATCH")))


SummarizedExperiment(sce)
colnames(rowData(sce))   #excess features/markers
dim(colData(sce))        #summary samples data

tmp <- table(colData(sce)$sample_id)    #numbers of events/cells per file
tmp_2 <- data.frame(tmp)
colnames(tmp_2)[1] <- c("sample_id")

tmp_2 <- left_join(tmp_2,md_info,by=join_by(sample_id))
kruskal.test(tmp_2$Freq~factor(tmp_2$BATCH))

tmp_2$BATCH <- factor(tmp_2$BATCH)
ggplot(tmp_2, aes(x=BATCH,y=Freq))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(color=tissue_type),width = 0.25)+
  scale_y_continuous(trans='log10')+ theme_bw()

write.table(tmp,"tmp.txt",sep="\t",quote=F,col.names=F)

# Step 3. Inspect the data
# Step 3.1. inspect Ref-PBMC samples
ref_sam <- grep("Ref",md_info$sample_id,value=T)
sce_ref <- sce[,sce$sample_id %in% ref_sam]
colData(sce_ref)$sample_id <- droplevels(colData(sce_ref)$sample_id)

# plot distribution of markers
p_ref <- plotExprs(sce_ref, color_by = "BATCH")

sce_ref <- sce[,sce$sample_id %in% ref_sam[8]]
plotExprs(sce_ref, color_by = "BATCH")

# plot multi-dimension scale (median marker intensities
pbMDS(sce_ref, color_by = "BATCH", label_by = "sample_id")

# heatmap of median marker intensities
# in CyTOF applications, a cosine distance shows good performance
# scale= "last": aggregate then scale & trim
plotExprHeatmap(sce_ref, scale = "last",
    row_anno = FALSE,
    col_anno = FALSE)

# identification with FlowSOM: using k = 20 metacluster
# id of som(i.e. 100 clusters) and metacluster are stored in cluster_ids
set.seed(1234)
sce_ref <- cluster(sce_ref, features = "type",
    xdim = 10, ydim = 10, maxK = 20, seed = 1234)
p1 <- plotExprHeatmap(sce_ref, features = "type", 
    by = "cluster_id", k = "meta20", 
    bars = TRUE, perc = TRUE)+theme(legend.position = "none")

names(cluster_codes(sce_ref))   #codes of clusters

##access specific clustering resolution, i.e. cellIDs
table(cluster_ids(sce_ref, "som100"))
table(cluster_ids(sce_ref, "meta20"))

# run t-SNE/UMAP on at most 500/1000 cells per sample
set.seed(1234)
sce_ref <- runDR(sce_ref, "UMAP", cells = 1e3, features = "type")
plotDR(sce_ref, "UMAP", color_by = "CD4")

# compare FSOM and UMAP clusters
plotDR(sce_ref, "UMAP", color_by = "meta20")

# ----------------
# Step 3.2. batch 2,3,4,5

p <- plotExprs(sce, color_by = "tissue_type")

plotDR(sce_ref, "UMAP", color_by = "meta20")

p <- plotDR(sce_ref, "UMAP",
  color_by = "meta15", 
  facet_by = "sample_id")













### import flowSet using FlowSOM::AggregateFlowFrames--> requires number of cells
##ff_agg <- FlowSOM::AggregateFlowFrames(fileNames = fcs_files,
##                                         cTotal = length(fcs_files) * cells_total,
##                                       verbose = TRUE,
##                                       writeMeta = TRUE,
##                                       writeOutput = TRUE,
##                                       outputFile = file.path(out_dir, paste0("aggregated_for_UMAP_analysis.fcs")))

## using CATALYST
tmp_panel <- panel_info[,c("fcs_colname","antigen","marker_class")]
tmpi <- which(grepl("function",panel_info$marker_class)==TRUE)
tmp_panel$marker_class[tmpi] <- c("state")
tmpi <- which(panel_info$marker_class=="lineage")
tmp_panel$marker_class[tmpi] <- c("type")

frames<-lapply(fcs_files,read.FCS) 
fs<-as(frames,"flowSet")
tmp <- fsApply(fs,keyword, "$FILENAME")


sce <- CATALYST::prepData(fcs_files,md = md,md_cols = list(file = "file", id = "sample_id", 
                   factors = c("tissue_type","patient_id")))
sce <- CATALYST::prepData(fcs_files,
               panel = tmp_panel,
               md = md,
               panel_cols = list(channel = "fcs_colname", antigen = "antigen", class = "marker_class"),
               md_cols = list(file = "file", id = "sample_id", 
                   factors = c("tissue_type","patient_id","BATCH")))



# Build UMAP on aggregated files
UMAP_res <- UMAP(fcs_files = files, 
                 clustering_markers = c("CD", "HLA", "IgD"),
                 functional_markers = c("IL", "TNF", "TGF", "Gr", "IFNa", "MIP", "MCP1"),
                 out_dir = analysis_dir,
                 batch_pattern = "day[0-9]*", 
                 arcsine_transform = TRUE, 
                 cells_total = 5000)



## import all data files in the fcs.dir --> list of flowFrame objects
frames <- lapply(dir(fcs.dir, full.names=TRUE), read.FCS)   ##import data
names(frames) <- sapply(frames, keyword, "SAMPLE ID")   ##extract "SAMPLE ID" parameter using keyword approach then mame it.

fs <- as(frames, "flowSet")   ##set as flowset object
sampleNames(fs)               ## check sample ID

phenoData(fs)$Filename <- fsApply(fs,keyword, "$FIL")    #create phenoData using $FIL (filename) keyword
pData(phenoData(fs))          ##access phenoData




# ------------------------------------------------------------------------------
# Plot batch effect ------------------------------------------------------------
# Not complete
#-------------------------------------------------------------------------------
selPanel <- c("TBNK")  #*******
sel4ref <- c(1:7)     #***
batch4extractPatt <- "(?i).*(batch[0-9]*).*.fcs"  #****
refPattern <- c("Ref.*_gated.fcs$")
dna_ch <- c("Ir191Di","Ir193Di")   #*****
via_ch <- c("Pt195Di")    #***Cisplatin

# get panel info
fin_panel <- paste(selPanel,"_markers_022625.txt",sep="")
panel_info <- read.delim(file.path("Ranalysis",fin_panel),
     sep="\t",header=T,stringsAsFactors=FALSE)
 
# get channel information for fsom
lineage_idx <- which(panel_info$marker_class=="lineage")  #****
lineage_markers <- panel_info$fcs_desc[lineage_idx]  
function_idx <- which(panel_info$marker_class=="function")  #****
function_markers <- panel_info$fcs_desc[function_idx] 
pheno_4fsom <- lineage_markers   #******c(dna_ch_4fsom,via_ch_4fsom,lineage_markers)    #*******

# Define batch pattern
batch_pattern <- "batch[0-9]*"   #****

# Define files before normalization and order them according to the batch 
# Set input directory for files before CytoNorm normalization 
gate_dir <- file.path(workFolder,"CYTOF_data", 
                  "Gated",selPanel)
files_before_norm <- list.files(gate_dir, 
                                pattern = ".fcs", 
                                full.names = T)

batch <- str_match(files_before_norm, batch_pattern)[,1]
files_before_norm <- files_before_norm[order(factor(batch))]

# Define files after normalization and order them according to the batch 
# Set input directory for files after CytoNorm normalization
norm_dir <- file.path(workFolder,"CYTOF_data", 
                  "CytoNormed",selPanel)
files_after_norm <- list.files(norm_dir, 
               pattern = ".fcs", 
               full.names = T)
batch <- stringr::str_match(files_after_norm, batch_pattern)[,1]
files_after_norm <- files_after_norm[order(factor(batch))]

# Plot batch effect 
plot_batch(files_before_norm = files_before_norm, 
           files_after_norm = files_after_norm,
           out_dir = norm_dir, 
           batch_pattern = batch_pattern, 
           clustering_markers = pheno_4fsom,
           manual_colors = c("darkorchid4", "darkorange", "chartreuse4"))

batch_pattern <- "batch[0-9]*"
plot_marker_quantiles(files_after_norm = files_after_norm, 
           files_before_norm = files_before_norm, 
           batch_pattern = batch_pattern, 
           arcsine_transform = TRUE,
           markers_to_plot = c(lineage_markers,function_markers),
           manual_colors = c("darkorchid4", "darkorange", "darkgreen"),
           out_dir = norm_dir)

# Use FlowSOM to extract cell frequency and MSI

# Create a list with files before and after normalization
all_files <- list("before" = files_before_norm,
                 "after" = files_after_norm)

# perform FlowSOM clustering and extract cell frequency and msi per cluster and metacluster
mx <- extract_pctgs_msi_per_flowsom(file_list = all_files, 
            nCells = 5000, 
            phenotyping_markers = pheno_4fsom,
            functional_markers = function_markers
            xdim = 10, 
            ydim = 10,
            n_metaclusters = 35, 
            out_dir = norm_dir, 
            arcsine_transform = TRUE)

# set title for each plot
title_gg <- c("Cluster frequencies" = "cl_pctgs",
              "Metacluster frequencies" = "mcl_pctgs",
              "Cluster MSIs" = "cl_msi",
              "Metacluster MSIs" = "mcl_msi")

# create the list to store the plots
plots <- list()
for (i in seq_along(title_gg)){
  df_plot <- prepare_data_for_plotting(frequency_msi_list = mx, 
                   matrix_type = title_gg[1],
                   n_neighbours = 11)
  
  # extract annotation for plotting
  df_plot$day <- str_match(rownames(df_plot), batch_pattern)[,1]
  df_plot$reference <- ifelse(grepl("Ref", rownames(df_plot)),"ref", "")
  df_plot$sample <- ifelse(grepl("p1", rownames(df_plot)),"p1", 
                           ifelse(grepl("p2", rownames(df_plot)), "p2", "ref"))
  df_plot$stimulation <- str_match(rownames(df_plot), "RSQ|LPS|IMQ|CPG|UNS")
  df_plot$normalization <- factor(ifelse(grepl("Norm", rownames(df_plot)),
                                         "Normalized", "Raw"), 
                                  levels = c("Raw", "Normalized"))
 
  # plot
  gg <- ggplot(df_plot, aes(x = dim1, y = dim2))+
    geom_point(data=df_plot, aes_string(x="dim1", y="dim2", fill = "day", 
                                        shape = "sample", color = "day"), 
               size = 3)+
    facet_wrap(~normalization)+
    ggtitle(names(title_gg)[which(title_gg%in% name)])+
    scale_shape_manual(values = c(22, 21, 24))+
    scale_fill_manual(values = c("darkorchid4", "darkorange", "chartreuse4"))+
    scale_color_manual(values = c("darkorchid4", "darkorange", "chartreuse4"))+
    theme(panel.background = element_rect(fill = "white", colour = "black",
                                          size = 1, linetype = "solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "right",
          legend.key=element_blank(),
          title = element_text(size = 10),
          strip.text = element_blank(), 
          strip.background = element_rect(fill = "white", colour = "white"))
  
  gg <- ggplotGrob(gg)
  plots[[name]] <- gg
}
gg_a <- ggarrange(plotlist = plots,
                  ncol = 2,
                  nrow = 2)

ggplot2::ggsave(filename ="batch_effect_frequency_MSI.png",
                device = "png",
                path = norm_dir,
                plot = gg_a,
                units = "cm",
                width = 19,
                height = 10, dpi = 300)


# *********************************************
# Ref

Ftab_quantile_ref <- reshape(quantiles, idvar = c("File","Marker","Batch","Sample"),
                      timevar = "Quantile", direction = "wide")

sel_marker4plot <- panel_info$fcs_desc[1:9]
tmp4plot <- subset(Ftab_quantile_ref,
    is.element(Ftab_quantile_ref$Marker,sel_marker4plot)
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





plot_marker_quantiles_SDL(select_fcs_files=file_a, 
                      batch_pattern = batch_pattern, 
                      arcsine_transform = TRUE,
                      remove_beads = TRUE,
                      bead_channel = "140", 
                      uncommon_prefix = "_beadNorm.fcs|.FCS", 
                      markers_to_plot = c("CD", "HLA", "IgD", "IL", "TNF",
                                          "TGF", "GR", "IFNa"),
                      manual_colors = c("darkorchid4", "darkorange", "darkgreen"),
                      out_dir = bead_norm_dir)



 quantiles <- cbind(quantiles, "Batch" = stringr::str_match(
    basename(as.character(quantiles$File)), batch4extractPatt)[,1])


ttmp <- stringr::str_match(basename(as.character(quantiles$File)), batch4extractPatt)[,1]

grep("191",names(norm_markers))

sel_marker4plot <- norm_markers[54]
p_x <- quantiles %>% 
    dplyr::filter(Normalization == "YES" & Marker==sel_marker4plot) %>%   
    ggplot2::ggplot(aes(x = Sample,
               y = Value,
               color = Batch)) +
    geom_point(aes(alpha = ifelse(Quantile == "0.5", 2, 0))) +
    geom_line(aes(alpha = ifelse(Quantile != "0.5", 2, 0),
                  group = interaction(Batch, Quantile),
                  size = ifelse(Quantile %in% c("0.05", "0.95"), 1,
                                ifelse(Quantile == "0.5", 0, 2))),
              alpha = 0.5) +
    ylab(label = sel_marker4plot)+  
    scale_size_identity() +
    scale_alpha_identity() +
    facet_wrap(~ Batch, ncol = ncols, scales = "free_x") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "bottom")
  
sel_marker4plot <- norm_markers[1:12]
p_x <- quantiles %>% 
    dplyr::filter(Normalization == "YES" & Marker %in% sel_marker4plot) %>%   
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
p_x

# ---------------------------------------------------------------------------
   if(!dir.exists(out_dir)) dir.create(out_dir)
  
  if (!is.null(file_batch_id)) {
    scores <- list()
    for (batch in unique(file_batch_id)){
      print(batch)
      
      files <- fcs_files[file_batch_id == batch]
      nCells = length(fcs_files)*nCells_thres
      fsom <- fsom_aof(fcs_files = files, 
                       phenotyping_markers = phenotyping_markers, 
                        nCells = nCells,
                       out_dir = out_dir, 
                       arcsine_transform = arcsine_transform,
                       nClus = nClus,
                       batch = batch)
    
      scores[[batch]] <- aof_scoring(fcs_files = files, 
                                     phenotyping_markers = phenotyping_markers,
                                     fsom = fsom, out_dir = out_dir, batch = batch)
    }
    
  } else {
    files <- fcs_files
    fsom <- fsom_aof(fcs_files = files, phenotyping_markers = phenotyping_markers, 
                     out_dir = out_dir, arcsine_transform = arcsine_transform, 
                     nClus = nClus,
                     batch = NULL)
    
    scores <- aof_scoring(fcs_files = files, 
                          phenotyping_markers = phenotyping_markers,
                          fsom = fsom, out_dir = out_dir, batch = NULL)
  }


# quality check

