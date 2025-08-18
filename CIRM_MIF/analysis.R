library(dplyr)


# *****************************************************
# Inspect & quality check data
# *****************************************************

# ----- PANEL 2, SET 1 --------------------------------
inputParentFolder <- r'(D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20241121 CIRM Panel 2 Set 1 analysis/20241121 CIRM Panel 2 Set 1 analysis ObjectData)'
name_summary_file <- r'(D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20241121 CIRM Panel 2 Set 1 analysis/20241121 CIRM Panel 2 Set 1 analysis Total_Summary_Results.csv)'
output_dir <- r'(D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20241121 CIRM Panel 2 Set 1 analysis/Output)'
fout <- c("listHaloOut_panel2_set1_08182025.txt")

# ----- PANEL 2, SET 2 --------------------------------
inputParentFolder <- r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 2 analysis\20241121 CIRM Panel 2 Set 2 analysis ObjectData)'
name_summary_file <- r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 2 analysis\20241121 CIRM Panel 2 Set 2 analysis.csv)'
output_dir <- r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 2 analysis\Output)'
fout <- c("listHaloOut_panel2_set2_08182025.txt")

# -----------------------------------------------------

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
tmps_slideID <- matrix(NA, nrow(tmps_cond1), 3)
for (i in c(1:nrow(tmps_cond1))) {
  tmp_str <- unlist(strsplit(tmps_cond1[i, 1], split=" "))
  
  patID <- tmp_str[[4]]
  day <- tmp_str[[5]]
  quadrant <- tmp_str[[6]]
  
  new_roi <- paste(patID, day, quadrant, sep="-")
  tmps_slideID[i, ] <- c(patID, day, new_roi)
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
                      rep(tmps_jobID[j], n_roi),
                      names(roi_info),
                      rep(tmps_slideID[j, 3], n_roi),
                      as.numeric(roi_info),
                      is_sum_avai,
                      tmp4sum)
    listFiles <- rbind(listFiles, tmp_info)
  } else {
    stop("roi_column_interest not found in image")
  }
}

colnames(listFiles) <- c("originalFileName", "PathologicalID", "BlockID",
                         "jobID", "originalROI", "newROI", "cellNumber", "IsSummaryAvai", "SumFileName")

write.table(listFiles, file=file.path(tmpout, fout), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, fileEncoding="UTF-8")

# *****************************************************
# Regenerate data
# *****************************************************

# ----- PANEL 2, SET 1 --------------------------------
inputParentFolder <- r'(D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20241121 CIRM Panel 2 Set 1 analysis/20241121 CIRM Panel 2 Set 1 analysis ObjectData)'
name_summary_file <- r'(D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20241121 CIRM Panel 2 Set 1 analysis/20241121 CIRM Panel 2 Set 1 analysis Total_Summary_Results.csv)'
fin_conversion <- r'(D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20241121 CIRM Panel 2 Set 1 analysis/Output/listHaloOut_panel2_set1_08182025.txt)'

fin_header_object <- r'(D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20241121 CIRM Panel 2 Set 1 analysis/object_header_panel2_set1_081825.csv)'
fin_header_summary <- r'(D:/CHRISTOPHER_BUI/CIRM_MIF/Data/summary_header_panel2_081825.csv)'

outputParentFolder <- r'(D:/CHRISTOPHER_BUI/CIRM_MIF/Data/20241121 CIRM Panel 2 Set 1 analysis/Regenderated Data)'
suffix_obj <- c("_panel2_v1_object_results.csv") #******
suffix_summary <- c("_panel2_v1_Summary_Results.csv") #******

# Step 1: check existing files
sampleInfo <- read.delim(fin_conversion, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Step 1.1: check quantified data
nfiles <- nrow(sampleInfo)
Ftab <- data.frame(fileName=sampleInfo[, 1], is_exist=rep(NA, nfiles))
for (i in c(1:nfiles)) {
  Ftab$is_exist[i] <- file.exists(sampleInfo[i, 1])
}
table(Ftab$is_exist)

# Step 1.2: check existence of summary data
tmp_sum <- read.csv(name_summary_file, stringsAsFactors=FALSE, header=TRUE)
dim(tmp_sum)

tmp_sumfile <- sampleInfo$SumFileName
sampleInfo$OriginalSummaryFile_fullPath <- tmp_sumfile
list4sum <- sampleInfo$OriginalSummaryFile_fullPath[!duplicated(sampleInfo$OriginalSummaryFile_fullPath)]
nsum <- length(list4sum)
Fsum <- data.frame(fileName=list4sum, is_exist=rep(NA,nsum))
for (i in c(1:nsum)){
  Fsum$is_exist[i] <- file.exists(list4sum[i])
}

# Step 1.3: check overlap between summary file & object results
table(is.element(tmp_sum$Analysis.Region, sampleInfo$originalROI))
table(is.element(sampleInfo$originalROI, tmp_sum$Analysis.Region))




# Step 2: rewrite object_result files & change names
# Step 2.1: get headers
obj_header <- read.csv(fin_header_object, header=TRUE, stringsAsFactors=FALSE)
sum_header <- read.csv(fin_header_summary, header=TRUE, stringsAsFactors=FALSE)

# Step 2.2: reload summary file 
data_summary <- read.csv(name_summary_file, stringsAsFactors=FALSE, header=TRUE)
ncol(data_summary)
nrow(sum_header)

nfiles <- nrow(sampleInfo)
Ftab_new <- data.frame(fileName_Old=sampleInfo$originalFileName, fileName_New=rep(NA, nfiles), summary_new=rep(NA, nfiles))

for (i in c(1:nfiles)) {
  fin_tmp <- sampleInfo$originalFileName[i]
  fin_sum <- sampleInfo$SumFileName[i]
  
  sampleID_old <- sampleInfo$originalROI[i]
  sampleID_new <- sampleInfo$newROI[i]
  sampleID_jobID <- sampleInfo$jobID[i]
  
  fout_tmp <- paste(outputParentFolder, "/", sampleID_new, suffix_obj, sep="")
  fout_sum <- paste(outputParentFolder, "/", sampleID_new, suffix_summary, sep="")
  Ftab_new$fileName_New[i] <- fout_tmp
  Ftab_new$summary_new[i] <- fout_sum
  
  # get object data
  data_in <- read.csv(fin_tmp, stringsAsFactors=FALSE, header=TRUE)
  tmp_selR <- which(data_in$Analysis.Region==sampleID_old)
  
  data_out <- data_in[tmp_selR, ]
  data_out$Analysis.Region <- sampleID_new
  
  # get summary data
  tmp_sum <- read.csv(fin_sum, stringsAsFactors=FALSE, header=TRUE)
  isMatchID <- tmp_sum$Analysis.Region==sampleID_old & tmp_sum$Job.Id==sampleID_jobID
  
}













