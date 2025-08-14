library(dplyr)

################################################################
# PANEL 2, SET 1
################################################################
# output directory
output_dir <- r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\Output)'

# folder with input files
inputParentFolder <- r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 1 analysis\20241121 CIRM Panel 2 Set 1 analysis ObjectData)'

# get file names
tmpin <- inputParentFolder
tmpf <- list.files(path=tmpin, pattern="object_results.csv", full.names=TRUE)
tmps <- list.files(path=tmpin, pattern="object_results.csv", full.names=FALSE)

# get summary file
summary_file <- paste0(r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 1 analysis\20241121 CIRM Panel 2 Set 1 analysis Total_Summary_Results.csv)')
data_summary <- read.csv(summary_file, stringsAsFactors = FALSE, header = TRUE)

# output file
fout <- c("listHaloOut_panel2_batch1_08142025.txt")

# extract meta data from file names
tmps_cond1 <- t(sapply(tmps, function(x) unlist(strsplit(x, split="\\."))[1:2]))

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



# CHECK IF COLUMN NAMES MATCH WITHIN BATCH ------------------------------------------

colnames_list <- list()
file_names <- basename(tmpf)

all_cols_check <- list()

sel_cols <- readLines(r'(C:\Users\cdbui\Desktop\panel2_columns.txt)')
sel_cols_check <- list()

reference_columns <- colnames(read.csv(tmpf[1], stringsAsFactor=FALSE, check.names=FALSE, nrows=1))
for (i in seq_along(tmpf)) {
  obj_res <- tmpf[i]
  column_names <- colnames(read.csv(obj_res, stringsAsFactors=FALSE, check.names=FALSE, nrows=1))
  
  all_cols_check[[file_names[i]]] <- identical(column_names, reference_columns)
  
  colnames_list[[file_names[i]]] <- column_names
  
  # check if object result's columns are in list of selected markers
  sel_cols_check[[obj_res]] <- sel_cols[!sel_cols %in% column_names]
}

for (i in names(all_cols_check)) {
  if (!all_cols_check[[i]]) {
    print(i)
  }
}
print(paste0(sum(unlist(all_cols_check)), "/", length(all_cols_check), " Files Matched"))

for (i in names(sel_cols_check)) {
  if (length(sel_cols_check[[i]]) > 0) {
    print(sel_cols_check[[i]])
  }
}



################################################################
# PANEL 2, SET 2
################################################################
# output directory
output_dir <- r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\Output)'

# folder with input files
inputParentFolder <- r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 2 analysis\20241121 CIRM Panel 2 Set 2 analysis ObjectData)'

# get file names
tmpin <- inputParentFolder
tmpf <- list.files(path=tmpin, pattern="object_results.csv", full.names=TRUE)
tmps <- list.files(path=tmpin, pattern="object_results.csv", full.names=FALSE)

# get summary file
summary_file <- paste0(r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 2 analysis\20241121 CIRM Panel 2 Set 2 analysis.csv)')
data_summary <- read.csv(summary_file, stringsAsFactors = FALSE, header = TRUE)

# output file
fout <- c("listHaloOut_panel2_batch2_08142025.txt")

# extract meta data from file names
tmps_cond1 <- t(sapply(tmps, function(x) unlist(strsplit(x, split="\\."))[1:2]))

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



# CHECK IF COLUMN NAMES MATCH WITHIN BATCH ------------------------------------------

colnames_list <- list()
file_names <- basename(tmpf)

all_cols_check <- list()

sel_cols <- readLines(r'(C:\Users\cdbui\Desktop\panel2_columns.txt)')
sel_cols_check <- list()

reference_columns <- colnames(read.csv(tmpf[1], stringsAsFactor=FALSE, check.names=FALSE, nrows=1))
for (i in seq_along(tmpf)) {
  obj_res <- tmpf[i]
  column_names <- colnames(read.csv(obj_res, stringsAsFactors=FALSE, check.names=FALSE, nrows=1))
  
  all_cols_check[[file_names[i]]] <- identical(column_names, reference_columns)
  
  colnames_list[[file_names[i]]] <- column_names
  
  # check if object result's columns are in list of selected markers
  sel_cols_check[[obj_res]] <- sel_cols[!sel_cols %in% column_names]
}

for (i in names(all_cols_check)) {
  if (!all_cols_check[[i]]) {
    print(i)
  }
}
print(paste0(sum(unlist(all_cols_check)), "/", length(all_cols_check), " Files Matched"))

for (i in names(sel_cols_check)) {
  if (length(sel_cols_check[[i]]) > 0) {
    print(sel_cols_check[[i]])
  }
}

################################################################
# CHECK IF COLUMN NAMES MATCH BETWEEN BATCHES
################################################################
inputFolderP2B1 <- r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 1 analysis\20241121 CIRM Panel 2 Set 1 analysis ObjectData)'
inputFolderP2B2 <- r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 2 analysis\20241121 CIRM Panel 2 Set 2 analysis ObjectData)'

tmpf1 <- list.files(inputFolderP2B1, pattern="object_results.csv", full.names=TRUE)
tmpf2 <- list.files(inputFolderP2B2, pattern="object_results.csv", full.names=TRUE)

column_names1 <- colnames(read.csv(tmpf1[1], stringsAsFactors=FALSE, check.names=FALSE))
column_names2 <- colnames(read.csv(tmpf2[1], stringsAsFactors=FALSE, check.names=FALSE))

batch_cols_check <- identical(column_names1, column_names2)
print(batch_cols_check)




################################################################
# LEVEL 3 REGENERATION
################################################################