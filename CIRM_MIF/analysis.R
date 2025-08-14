library(dplyr)

################################################################
# PANEL 2, SET 1
################################################################
# folder with input files
inputParentFolder <- r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 1 analysis\20241121 CIRM Panel 2 Set 1 analysis ObjectData)'

# get file names
tmpin <- inputParentFolder
tmpf <- list.files(path=tmpin, pattern="object_results.csv", full.names=TRUE)
tmps <- list.files(path=tmpin, pattern="object_results.csv", full.names=FALSE)

# get summary file
summary_file <- paste0(r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 1 analysis\20241121 CIRM Panel 2 Set 1 analysis Total_Summary_Results.csv)')
data_summary <- read.csv(summary_file, stringsAsFactors = FALSE, header = TRUE)

# extract meta data from file names
tmps_cond1 <- t(sapply(tmps, function(x) unlist(strsplit(x, split="\\."))[1:2]))


# CHECK IF EACH FILE NAME MATCHES DATA IN SUMMARY FILE ----------------

# image tag
img_tag <- paste0(tmps_cond1[, 1], ".tif")
match_imgtag <- img_tag %in% data_summary$Image.Tag

# job id
job_id <- sapply(tmps_cond1[, 2], function(x) {
  sub(".*_job(\\d+)CIRM.*", "\\1", x)
})
match_jobid <- job_id %in% data_summary$Job.Id

# algorithm name
algo_name <- sapply(tmps_cond1[, 2], function(x) {
  sub(".*(CIRM.*)_object_results.*", "\\1", x)
})
match_algo <- algo_name %in% data_summary$Algorithm.Name

print(paste0(sum(match_imgtag), "/", length(match_imgtag), " Files Matched"))
print(paste0(sum(match_jobid), "/", length(match_jobid), " Files Matched"))
print(paste0(sum(match_algo), "/", length(match_algo), " Files Matched"))


# CHECK IF COLUMN NAMES MATCH ------------------------------------------

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

# TEST RUN TEMPLATE SCRIPT ------------------------------------
tmpin <- inputParentFolder
tmpf <- list.files(path=tmpin, pattern="object_results", full.names=TRUE)
tmps <- list.files(path=tmpin, pattern="object_results", full.names=FALSE)

fin_summary_file




################################################################
# PANEL 2, SET 2
################################################################
# folder with input files
inputParentFolder <- r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 2 analysis\20241121 CIRM Panel 2 Set 2 analysis ObjectData)'

# get file names
tmpin <- inputParentFolder
tmpf <- list.files(path=tmpin, pattern="object_results.csv", full.names=TRUE)
tmps <- list.files(path=tmpin, pattern="object_results.csv", full.names=FALSE)

# get summary file
summary_file <- paste0(r'(D:\CHRISTOPHER_BUI\CIRM_MIF\Data\20241121 CIRM Panel 2 Set 2 analysis\20241121 CIRM Panel 2 Set 2 analysis.csv)')
data_summary <- read.csv(summary_file, stringsAsFactors = FALSE, header = TRUE)

# extract meta data from file names
tmps_cond1 <- t(sapply(tmps, function(x) unlist(strsplit(x, split="\\."))[1:2]))


# CHECK IF EACH FILE NAME MATCHES DATA IN SUMMARY FILE ----------------

# image tag
img_tag <- paste0(tmps_cond1[, 1], ".tif")
match_imgtag <- img_tag %in% data_summary$Image.Tag

# job id
job_id <- sapply(tmps_cond1[, 2], function(x) {
  sub(".*_job(\\d+)CIRM.*", "\\1", x)
})
match_jobid <- job_id %in% data_summary$Job.Id

# algorithm name
algo_name <- sapply(tmps_cond1[, 2], function(x) {
  sub(".*(CIRM.*)_object_results.*", "\\1", x)
})
match_algo <- algo_name %in% data_summary$Algorithm.Name

print(paste0(sum(match_imgtag), "/", length(match_imgtag), " Files Matched"))
print(paste0(sum(match_jobid), "/", length(match_jobid), " Files Matched"))
print(paste0(sum(match_algo), "/", length(match_algo), " Files Matched"))


# CHECK IF COLUMN NAMES MATCH ------------------------------------------

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


