library(dplyr)


# ******************************************************************************
# Check object files have same columns & in same order within and across batches
# *****************************************************************************
inputParentFolder1 <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/",
                            "20241121 CIRM Panel 2 Set 1 analysis/",
                            "20241121 CIRM Panel 2 Set 1 analysis ObjectData", sep="")
fout1 <- c("check_columns_panel2_batch1.txt")

inputParentFolder2 <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/",
                            "20241121 CIRM Panel 2 Set 2 analysis/",
                            "20241121 CIRM Panel 2 Set 2 analysis ObjectData", sep="")
fout2 <- c("check_columns_panel2_batch2.txt")


output_dir <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20250819_CIRM_Panel2_Analysis")


objects_batch1 <- list.files(inputParentFolder1, pattern="object_results.csv", full.names=TRUE)
objects_batch2 <- list.files(inputParentFolder2, pattern="object_results.csv", full.names=TRUE)

# for comparing columns across batches
all_objects <- c(objects_batch1, objects_batch2)
fout_all <- c("check_columns_across_batches.txt")

# check if columns match within batch
tmpout <- output_dir
# tmp_fout <- fout1 #******
# tmp_fout <- fout2 #******
tmp_fout  <- fout_all #******

# tmp_objects <- objects_batch1 #******
# tmp_objects <- objects_batch2 #******
tmp_objects <- all_objects

column_names <- lapply(tmp_objects, function(x) names(read.csv(x, nrows=1, stringsAsFactors=TRUE, header=TRUE, check.names=FALSE)))

reference <- column_names[[1]]

match_status <- sapply(seq_along(column_names), function(i) {
  if (i == 1) {
    return("ref")
  } else if (identical(column_names[[i]], reference)) {
    return("yes")
  } else {
    return("no")
  }
})

match_table <- data.frame(
  originalFilePath = tmp_objects,
  match = match_status
)

write.table(match_table, file=file.path(tmpout, tmp_fout), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, fileEncoding="UTF-8")



# *****************************************************
# Inspect & quality check data
# *****************************************************

# ----- PANEL 2, SET 1 ---------------------------------------------------------
inputParentFolder <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/",
                           "20241121 CIRM Panel 2 Set 1 analysis/",
                           "20241121 CIRM Panel 2 Set 1 analysis ObjectData", sep="")
name_summary_file <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/",
                           "20241121 CIRM Panel 2 Set 1 analysis/",
                           "20241121 CIRM Panel 2 Set 1 analysis Total_Summary_Results.csv", sep="")
output_dir <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20250819_CIRM_Panel2_Analysis")
fout <- c("qcFile_panel2_batch1_08192025.txt")

# ----- PANEL 2, SET 2 ---------------------------------------------------------
inputParentFolder <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/",
                           "20241121 CIRM Panel 2 Set 2 analysis/",
                           "20241121 CIRM Panel 2 Set 2 analysis ObjectData", sep="")
name_summary_file <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/",
                           "20241121 CIRM Panel 2 Set 2 analysis/",
                           "20241121 CIRM Panel 2 Set 2 analysis.csv", sep="")
output_dir <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20250819_CIRM_Panel2_Analysis")
fout <- c("qcFile_panel2_batch2_08192025.txt")

# ------------------------------------------------------------------------------


# get file names
tmpin <- inputParentFolder
tmpf <- list.files(path=tmpin, pattern="object_results.csv", full.names=TRUE)
tmps <- list.files(path=tmpin, pattern="object_results.csv", full.names=FALSE)
tmpout <- output_dir

# get summary file
summary_file <- name_summary_file
data_summary <- read.csv(summary_file, stringsAsFactors=FALSE, header=TRUE)

# extract meta data from file names
tmps_cond1 <- t(sapply(tmps, function(x) unlist(strsplit(x, split="\\."))[1:2]))

# new ROI
tmps_slideID <- matrix(NA, nrow(tmps_cond1), 5)
for (i in c(1:nrow(tmps_cond1))) {
  tmp_str <- unlist(strsplit(tmps_cond1[i, 1], split=" "))
  if (length(tmp_str) > 2) {
    # PANEL 2, BATCH 1
    # IN OBJECT FILE NAME, THERE IS SPACE: 'hpanel 2'
    # tmps_slideID[i, 1] <- paste0(tmp_str[[2]], tmp_str[[3]]) # panel
    # tmps_slideID[i, 2] <- tmp_str[[4]] #******               # patientID
    # tmps_slideID[i, 3] <- tmp_str[[5]] #******               # time
    # tmps_slideID[i, 4] <- tmp_str[[6]] #******               # quadrant
    # 
    # new_roi <- paste(tmp_str[[4]], tmp_str[[5]], tmp_str[[6]], sep="_") #******
    # tmps_slideID[i, 5] <- new_roi
    
    # PANEL 2, BATCH 2
    # IN OBJECT FILE NAME, THERE IS NO SPACE: 'hpanel2'
    tmps_slideID[i, 1] <- tmp_str[2] #******  # panel
    tmps_slideID[i, 2] <- tmp_str[3] #******  # patientID
    tmps_slideID[i, 3] <- tmp_str[4] #******  # time
    tmps_slideID[i, 4] <- tmp_str[5] #******  # quadrant
    new_roi <- paste(tmp_str[[3]], tmp_str[[4]], tmp_str[[5]], sep="_") #******
    tmps_slideID[i, 5] <- new_roi

  } else {
    tmps_slideID[i, ] <- tmp_str
  }
}
rownames(tmps_slideID) <- NULL

tmps_jobID <- sapply(tmps_cond1[, 2], function(x) substr(x, 8, 11))
names(tmps_jobID) <- NULL

nfiles <- length(tmpf)
listFiles <- NULL
roi_column_interest <- c("Analysis.Region")

for (j in c(1:nfiles)) {
  tmp_data_in <- read.csv(tmpf[j], stringsAsFactors=FALSE, header=TRUE)
  image_colnames <- colnames(tmp_data_in)
  if (is.element(roi_column_interest, image_colnames)) { # if roi column name exists in object's columns
    roi_info <- table(tmp_data_in[, roi_column_interest]) # table of values in object's roi column
    n_roi <- length(roi_info) # number of unique roi in object
    is_sum_avai <- is.element(names(roi_info), data_summary$Analysis.Region) # if any of object's unique roi values appear in summary's roi column
    tmp4sum <- ifelse(is_sum_avai==TRUE, name_summary_file, "NA")
    tmp_info <- cbind(rep(tmpf[j], n_roi),
                      rep(tmps_slideID[j, 1], n_roi),
                      rep(tmps_slideID[j, 2], n_roi),
                      rep(tmps_slideID[j, 3], n_roi),
                      rep(tmps_slideID[j, 4], n_roi),
                      rep(tmps_jobID[j], n_roi),
                      names(roi_info),
                      rep(tmps_slideID[j, 5], n_roi),
                      as.numeric(roi_info),
                      is_sum_avai,
                      tmp4sum)
    listFiles <- rbind(listFiles, tmp_info)
  } else {
    stop("roi_column_interest not found in image")
  }
}

colnames(listFiles) <- c("originalFileName", "Panel", "PatientID", "Time", "Quadrant",
                         "jobID", "originalROI", "newROI", "cellNumber", "IsSummaryAvai", "SumFileName")

write.table(listFiles, file=file.path(tmpout, fout), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, fileEncoding="UTF-8")

# combine
df1 <- read.table("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20250819_CIRM_Panel2_Analysis/qcFile_panel2_batch1_08192025.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
df2 <- read.table("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20250819_CIRM_Panel2_Analysis/qcFile_panel2_batch2_08192025.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
df1$Batch <- 1
df2$Batch <- 2
combined_qc <- rbind(df1, df2)
write.table(combined_qc, file=file.path(output_dir, "qcFile_panel2_all_batch_081925.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, fileEncoding="UTF-8")



# ******************************************************************************
# Regenerate data
# ******************************************************************************

# ----- PANEL 2, SET 1 ---------------------------------------------------------
inputParentFolder <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/",
                           "20241121 CIRM Panel 2 Set 1 analysis/",
                           "20241121 CIRM Panel 2 Set 1 analysis ObjectData", sep="")
output_dir <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20250819_CIRM_Panel2_Analysis")
name_summary_file <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/",
                           "20241121 CIRM Panel 2 Set 1 analysis/",
                           "20241121 CIRM Panel 2 Set 1 analysis Total_Summary_Results.csv", sep="")
qcFilePath <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20250819_CIRM_Panel2_Analysis/",
                    "qcFile_panel2_batch1_08192025.txt", sep="")


# ----- PANEL 2, SET 2 ---------------------------------------------------------
inputParentFolder <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/",
                           "20241121 CIRM Panel 2 Set 2 analysis/",
                           "20241121 CIRM Panel 2 Set 2 analysis ObjectData", sep="")
output_dir <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20250819_CIRM_Panel2_Analysis")
name_summary_file <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/",
                           "20241121 CIRM Panel 2 Set 2 analysis/",
                           "20241121 CIRM Panel 2 Set 2 analysis.csv", sep="")
qcFilePath <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20250819_CIRM_Panel2_Analysis/",
                    "qcFile_panel2_batch2_08192025.txt", sep="")

# ------------------------------------------------------------------------------

fin_header_object <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20250819_CIRM_Panel2_Analysis/",
                           "object_header_panel2_081925.csv", sep="")
fin_header_summary <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20250819_CIRM_Panel2_Analysis/",
                            "summary_header_panel2_081925.csv", sep="")
fin_conversion <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20250819_CIRM_Panel2_Analysis/",
                        "MIF_panel2hu_all_batch_fileConversion_081925.txt", sep="")

# fout_conv <- c("MIF_panel2hu_batch1_oldnew_QuantifiedFile.txt")
fout_conv <- c("MIF_panel2hu_batch2_oldnew_QuantifiedFile.txt")


outputParentFolder <- paste("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20250819_CIRM_Panel2_Analysis/Regenerated_Data")
suffix_obj <- c("_hpanel2_object_results.csv") #******
suffix_summary <- c("_hpanel2_Summary_Results.csv") #******

# ------------------------------------------------------------------------------

# Create Conversion File
batch_num <- 2 #******

qcFile <- read.table(qcFilePath, sep="\t", stringsAsFactors=FALSE, header=TRUE)
qcFile <- qcFile[qcFile$IsSummaryAvai == TRUE, ]

conversion_file <- qcFile %>%
  select(originalFileName, SumFileName, originalROI, newROI, Panel, PatientID, Time, Quadrant, jobID, cellNumber)
conversion_file$SumFileName <- basename(conversion_file$SumFileName)
conversion_file$Batch <- batch_num

tmpout <- output_dir

# write.table(conversion_file, file=file.path(tmpout, "MIF_panel2hu_batch1_fileConversion_081925.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, fileEncoding="UTF-8")
write.table(conversion_file, file=file.path(tmpout, "MIF_panel2hu_batch2_fileConversion_081925.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, fileEncoding="UTF-8")


# combine
df1 <- read.table(file.path(output_dir, "MIF_panel2hu_batch1_fileConversion_081925.txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
df2 <- read.table(file.path(output_dir, "MIF_panel2hu_batch2_fileConversion_081925.txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
combined_conversion <- rbind(df1, df2)
write.table(combined_conversion, file=file.path(output_dir, "MIF_panel2hu_all_batch_fileConversion_081925.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, fileEncoding="UTF-8")

# ------------------------------------------------------------------------------

# Step 1: check existing files
sampleInfo <- read.delim(fin_conversion, sep="\t", header=TRUE, stringsAsFactors=FALSE)

sampleInfo <- sampleInfo %>%
  filter(Batch == 2)

# Step 1.1: check quantified data
nfiles <- nrow(sampleInfo)
Ftab <- data.frame(fileName=sampleInfo[, 1], is_exist=rep(NA, nfiles))
for (i in c(1:nfiles)) {
  Ftab$is_exist[i] <- file.exists(sampleInfo[i, 1])
}
table(Ftab$is_exist)

# Step 1.2: check existence of summary data
tmp_sumfile <- paste0(dirname(dirname(sampleInfo$originalFileName)), "/", sampleInfo$SumFileName, sep="")
sampleInfo$OriginalSummaryFile_fullPath <- tmp_sumfile
list4sum <- sampleInfo$OriginalSummaryFile_fullPath[!duplicated(sampleInfo$OriginalSummaryFile_fullPath)]
nsum <- length(list4sum)
Fsum <- data.frame(fileName=list4sum, is_exist=rep(NA,nsum))
for (i in c(1:nsum)){
  Fsum$is_exist[i] <- file.exists(list4sum[i])
}
table(Fsum$is_exist)


# Step 1.3 insepct summary data for ORIGINAL ROIs and job IDs
Fsum_roi <- NULL
for (i in c(1:nsum)) {
  tmp_sum <- read.csv(list4sum[i], stringsAsFactors=FALSE, header=TRUE)
  tmp <- cbind(rep(list4sum[i], nrow(tmp_sum)), tmp_sum[, c("Analysis.Region", "Job.Id")])
  Fsum_roi <- rbind(Fsum_roi, tmp)
}
colnames(Fsum_roi) <- c("fileName", "Analysis.Region", "Job.Id")


# Step 1.4: check overlap between summary files & object results
nrow(Fsum_roi)
nrow(sampleInfo)
table(duplicated(Fsum_roi$Analysis.Region))

table(is.element(Fsum_roi$Analysis.Region, sampleInfo$originalROI))
table(is.element(sampleInfo$originalROI, Fsum_roi$Analysis.Region))

table(is.element(Fsum_roi$Job.Id, sampleInfo$jobID))
table(is.element(sampleInfo$jobID, Fsum_roi$Job.Id))
Fsum_roi[!is.element(Fsum_roi$Job.Id, sampleInfo$jobID), ]



# Step 2: Upload header files (summary & object) and prepare for level 3 regeneration
# Step 2.1: get headers
obj_header <- read.csv(fin_header_object, header=TRUE, stringsAsFactors=FALSE)
sum_header <- read.csv(fin_header_summary, header=TRUE, stringsAsFactors=FALSE)


# Step 3: Rewrite the object_result files and change names
nfiles <- nrow(sampleInfo)
Ftab_new <- data.frame(fileName_Old=sampleInfo$originalFileName, fileName_New=rep(NA, nfiles), summary_new=rep(NA, nfiles))

for (i in c(1:nfiles)) {
  message(paste0("Doing ", i, "/", nfiles, "..."))
  
  fin_tmp <- sampleInfo$originalFileName[i]
  fin_sum <- sampleInfo$OriginalSummaryFile_fullPath[i]
  
  sampleID_new <- sampleInfo$newROI[i]
  sampleID_old <- sampleInfo$originalROI[i]
  sampleID_jobID <- sampleInfo$jobID[i]
  sampleID_local <- paste(sampleInfo$newROI[i], ".tif", sep="")
  
  fout_tmp <- paste(outputParentFolder, "/", sampleID_new, "_job", sampleID_jobID, suffix_obj, sep="")
  fout_sum <- paste(outputParentFolder, "/", sampleID_new, "_job", sampleID_jobID, suffix_summary, sep="")
  Ftab_new$fileName_New[i] <- fout_tmp
  Ftab_new$summary_new[i] <- fout_sum
  
  # get object data
  data_in <- read.csv(fin_tmp, stringsAsFactors=FALSE, header=TRUE)
  tmp_selR <- which(data_in$Analysis.Region==sampleID_old)
  
  data_out <- data_in[tmp_selR, ]
  data_out$Analysis.Region <- sampleID_new
  data_out[, 1] <- sampleID_local
  
  # get summary data
  tmp_sum <- read.csv(fin_sum, stringsAsFactors=FALSE, header=TRUE)
  isMatchID <- tmp_sum$Analysis.Region==sampleID_old & tmp_sum$Job.Id==sampleID_jobID
  data_sum <- tmp_sum[which(isMatchID == TRUE), ]
  data_sum[, 1] <- sampleID_local # image location
  data_sum[, 2] <- sampleID_local # image tag
  data_sum$Analysis.Region <- sampleID_new
  
  
  # select
  tmp_obj_headers <- obj_header[,1]
  keptCol_obj <- which(obj_header$isKept==1) #******
  # data_out$Classifier.Label <- c("Tissue") # original labeled as "Tumor"
  
  tmp_sum_headers <- sum_header[,1]
  keptCol_sum <- which(sum_header$isKept==1) #******
  
  colnames(data_out) <- tmp_obj_headers
  new_data_out <- data_out[, keptCol_obj]
  
  colnames(data_sum) <- tmp_sum_headers
  new_data_sum <- data_sum[, keptCol_sum]
  
  # 8/19/25
  write.csv(new_data_out, fout_tmp, quote=FALSE, row.names=FALSE, fileEncoding="UTF-8")
  write.csv(new_data_sum, fout_sum, quote=FALSE, row.names=FALSE, fileEncoding="UTF-8")
}

write.table(Ftab_new, file.path(outputParentFolder, fout_conv), sep="\t", quote=FALSE, row.names=FALSE)


# df1 <- read.table("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20250819_CIRM_Panel2_Analysis/Regenerated_Data/MIF_panel2hu_batch1_oldnew_QuantifiedFile.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
# df2 <- read.table("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20250819_CIRM_Panel2_Analysis/Regenerated_Data/MIF_panel2hu_batch2_oldnew_QuantifiedFile.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
# df1$Batch <- 1
# df2$Batch <- 2
# combined_Ftab <- rbind(df1, df2)
# 
# write.table(combined_Ftab, file=file.path(outputParentFolder, "MIF_panel2hu_all_batch_oldnew_QuantifiedFile_081925.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, fileEncoding="UTF-8")


# ******************************************************************************
# Check marker combinations
# ******************************************************************************
library(SingleCellExperiment) 
library(reshape2)
library(dplyr)
library(apcluster)
library(scater)
library(dbscan)   #for identify)cell_communities
library(pheatmap)

source(paste0("D:/CHRISTOPHER_BUI/CIRM_MIF/Rscript/SPIAT_modified_LT_122324.R"))

output_dir <- c("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20250819_CIRM_Panel2_Analysis")

# Step 0: Setting workplace and input folders
# Step 0.1: check if input files and those in sampleInfo match
inputParentFolder <- c("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20250819_CIRM_Panel2_Analysis/Regenerated_Data/")

fin_clin <- c("D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20250819_CIRM_Panel2_Analysis/MIF_panel2hu_sampleInfo_081925.txt")
clin.data <- read.delim(fin_clin, header=TRUE, stringsAsFactors=FALSE)
clin.data <- clin.data[order(clin.data$fileName_New), ] # need to sort, since tmpf sorts

tmp_sc_files <- NULL

tmpin <- inputParentFolder
tmpf <- list.files(tmpin, pattern="object_results.csv", full.names=TRUE)
tmp_sc_files <- c(tmp_sc_files, tmpf)

sum(clin.data$fileName_New == tmp_sc_files)


# Step 0.2: setting input, using files listed in fin_clin
# Step 0.2.1: setting channel info and phenotype
channelInfo_fpath <- paste(inputParentFolder,
                           "Panel2_Human_ChannelInformation.csv", sep="")

# Step 1: Import single cell
# Step 1.1: import and check all positive combination
allPhenotypes_mks <- NULL
nsamples <- nrow(clin.data)
for (i in c(1:nsamples)) {
  message(paste0("Doing ", i, "/", nsamples, "..."))
  image_fpath <- clin.data$fileName_New[i]
  sce <- format_halo_to_sce_4DL(image_fpath,
                                channelInfo_fpath,
                                haloPhenotype_fpath=NULL)
  
  # remove duplicated rows in some files
  tmpdup <- duplicated(colnames(sce))
  if (sum(tmpdup) > 0) {
    sce <- sce[, !tmpdup]
  }
  list_phenotype <- table(colData(sce)$Phenotype)
  mat_phenotype <- as.data.frame(list_phenotype)
  colnames(mat_phenotype) <- c("phenotype", clin.data$newROI[i])
  rownames(mat_phenotype) <- NULL
  
  if (i==1) {
    allPhenotypes_mks <- mat_phenotype
  } else {
    tmp_old <- allPhenotypes_mks
    allPhenotypes_mks <- full_join(tmp_old, mat_phenotype,
                                   by=c("phenotype"="phenotype"))
  }
}
dim(allPhenotypes_mks)

tmp2 <- t(t(allPhenotypes_mks[,-1])/colSums(allPhenotypes_mks[,-1],na.rm=TRUE))
rownames(tmp2) <- allPhenotypes_mks[,1]
tmp2 <- rbind(tmp2,colSums(allPhenotypes_mks[,-1],na.rm=TRUE))

tmp_out <- cbind(c(as.character(allPhenotypes_mks[,1]),"totalCells"), tmp2)
colnames(tmp_out) <- c("markers", colnames(tmp2))

write.table(tmp_out, file.path(output_dir, "panel2_possibleCombination_082025.txt"), sep="\t", quote=FALSE, row.names=FALSE)








