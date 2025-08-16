library(dplyr)


####### PANEL 2, SET 1 ################################
inputParentFolder <- r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 1 analysis\20241121 CIRM Panel 2 Set 1 analysis ObjectData)'
name_summary_file <- r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 1 analysis\20241121 CIRM Panel 2 Set 1 analysis Total_Summary_Results.csv)'
output_dir <- r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 1 analysis\Output)'
fout <- c("listHaloOut_panel2_set1_08152025.txt")

####### PANEL 2, SET 2 ################################
inputParentFolder <- r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 2 analysis\20241121 CIRM Panel 2 Set 2 analysis ObjectData)'
name_summary_file <- r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 2 analysis\20241121 CIRM Panel 2 Set 2 analysis.csv)'
output_dir <- r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 2 analysis\Output)'
fout <- c("listHaloOut_panel2_set2_08152025.txt")

#######################################################

# get file names
tmpin <- inputParentFolder
tmpf <- list.files(path=tmpin, pattern="object_results.csv", full.names=TRUE)
tmps <- list.files(path=tmpin, pattern="object_results.csv", full.names=FALSE)
tmpout <- output_dir

# get summary file
summary_file <- name_summary_file
data_summary <- read.csv(summary_file, stringsAsFactors = FALSE, header = TRUE)

# extract meta data from file names
tmps_cond1 <- t(sapply(tmps, function(x) unlist(strsplit(x, split="\\."))[1:2]))


# ------------------------------------------------------------------------------
# new ROI
tmps_roi <- matrix(NA, nrow(tmps_cond1), 2)
for (i in c(1:nrow(tmps_cond1))) {
  tmp_str <- unlist(strsplit(tmps_cond1[i, 1], split=" "))
  
  patID <- tmp_str[4]
  day <- tmp_str[5]
  quadrant <- tmp_str[6]
  
  new_roi <- paste(patID, day, quadrant, sep="-")
  tmps_roi[i] <- new_roi
}
rownames(tmps_roi) <- NULL

tmps_jobID <- sapply(tmps_cond1[, 2], function(x) substr(x, 5, 11))
names(tmps_jobID) <- NULL

nfiles <- length(tmpf)
listFiles <- NULL
roi_column_interest <- c("Analysis.Region")

for (j in c(1:nfiles)) {
  tmp_data_in <- read.csv(tmpf[j], stringsAsFactors=FALSE, header=TRUE)
  image_colnames <- colnames(tmp_data_in)
  if (is.element(roi_column_interest, image_colnames)) {
    roi_info <- table(tmp_data_in[, roi_column_interest])
    n_roi <- length(roi_info)
    is_sum_avai <- is.element(names(roi_info), data_summary$Analysis.Region)
    tmp4sum <- ifelse(is_sum_avai==TRUE, name_summary_file, "NA")
    tmp_info <- cbind(rep(tmpf[j], n_roi))
  }
}

# ------------------------------------------------------------------------------

# extract slide ID: different format from first two batches and others
##tmps_slideID <- t(sapply(tmps_cond1[,1],function(x)
##                   unlist(strsplit(x,split=" "))[3:4]))
tmps_slideID <- matrix(NA,nrow(tmps_cond1),2)
for (i in c(1:nrow(tmps_cond1))){
  tmp_str <- unlist(strsplit(tmps_cond1[i,1],split=" ")) 
  if (length(tmp_str)>2){
    tmps_slideID[i,] <- tmp_str[3:4]   #*****
  }else{
    tmps_slideID[i,] <- tmp_str
  }
}
rownames(tmps_slideID) <- NULL

tmps_jobID <- sapply(tmps_cond1[,2],function(x)
  substr(x,5,11))

names(tmps_jobID) <- NULL


nfiles <- length(tmpf)
listFiles <- NULL
roi_column_interest <- c("Analysis.Region")

for (j in c(1:nfiles)){
  tmp_data_in <- read.csv(tmpf[j],stringsAsFactors=FALSE,header=T)
  image_colnames <- colnames(tmp_data_in)
  if (is.element(roi_column_interest,image_colnames)){
    roi_info <- table(tmp_data_in[,roi_column_interest])
    n_roi <- length(roi_info)
    # check existence of ROI in the summary file
    is_sum_avai <- is.element(names(roi_info),data_summary$Analysis.Region)
    tmp4sum <- ifelse(is_sum_avai==TRUE,summary_file,"NA")
    tmp_info <- cbind(rep(tmpf[j],n_roi),
                      rep(tmps_slideID[j,1],n_roi),
                      rep(tmps_slideID[j,2],n_roi),
                      rep(tmps_jobID[j],n_roi),
                      names(roi_info),as.numeric(roi_info),
                      is_sum_avai,tmp4sum)
    listFiles <-  rbind(listFiles,tmp_info)
    
  } else {
    stop("roi_column_interest not found in image")
  }
}

colnames(listFiles) <- c("originalFileName","PathologicalID","BlockID",
                         "jobID","originalROI","cellNumber","IsSummaryAvai","SumFileName")

write.table(listFiles, file.path(output_dir, fout), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)




tmp_sum <- data_summary
# Step 1.3: Check overlap between summary files and object_results
# !!! Need Sample Info !!!
table(is.element(tmp_sum$Analysis.Region,sampleInfo$Original_HTAN_ROI_ID))
table(is.element(sampleInfo$Original_HTAN_ROI_ID,tmp_sum$Analysis.Region))
