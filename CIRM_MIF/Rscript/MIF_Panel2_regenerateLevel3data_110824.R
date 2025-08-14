# ***********************************************************
# Convert table (from table(...) to data frame w/o stack:

tmpAdj_short <- sign(table(Ftab_all2$SampleID[tmpidx],Ftab_all2$SYMBOLS[tmpidx]))
tmpAdj_short2 <- data.frame(rbind(tmpAdj_short))   #*****

# **********************************************************************************
# Regenerate meta data for Panel 2 data--> include PCGA2 ID
# First, list directories since each batch is a folder --> multiple folders
# Using output from version 1.2
# 11/8/2024
# **********************************************************************************
# Run on R version 4.3.3
library(dplyr)

workFolder = c("D:/MIF/Human_SU2C/Ranalysis/")   #*****
setwd(workFolder)

# ----********non-GGO cohort ******************************
inputParentFolder <- paste("E:/",
                        "MIF_UCLA/HALO_out/Panel2_SU2C/",sep="")

fout <- c("tmp.txt")  #****

# -------
tmp_folder <- list.dirs(path = inputParentFolder, full.names = TRUE, recursive = TRUE)
write.table(tmp_folder,fout,sep="\t",quote=F,
             row.names=F,col.names=T)

tmp_summ <- list.files(path = inputParentFolder, pattern = "Total_Summary",
                     full.names = TRUE, recursive = TRUE)
write.table(tmp_summ,"tmp_sumName.txt",sep="\t",quote=F,
             row.names=F,col.names=T)

# ***************************************************************************
# Samples from SSW-17-15822 - A2, A3, A6, and H3 have wrong ID
# It should be HTA3_8022_xxxx intead of HTA3_8023_xxxx
# --> add wrongID in name tag --> change ID when regenerating data
# in the previous anaysis in 2022, changed ID and re-saved --> the file cannot
# be read by R 4.3.3 
# 11/8/24
# **************************************************************************

# **********************************************************************************
# Regenerate meta data for Panel 2 data--> include PCGA2 ID
# Scond, list files since each batch is a folder --> multiple folders
# Using output from version 1.2
# re-save the summary files using CSV-UTF-8 (comma delimitted)
# Panel2_SU2C/Human Panel 2 Set 4 extra v1.2/ObjectData/Total_Summary_Results.csv
# Panel2_SU2C/Human Panel 2 Set 4 v1.2/Total_Summary_Results.csv
# 11/8/2024
# **********************************************************************************
# Run on R version 4.3.3
library(dplyr)

workFolder = c("D:/MIF/Human_SU2C/Ranalysis/")   #*****
setwd(workFolder)

# **************non-GGO, multiple folders before 12/2023 ********************
##inputGrandParentFolder <- paste("E:/",
##                        "MIF_UCLA/HALO_out/Panel2_SU2C/",sep="")

## input text file include full path to the folder (col 1), and summary file (basename - col2)
name_directories <- c("listHaloOutDirs_panel2_nonGGO_pre2023.txt")


# ------------------------------
listDirs <- read.delim(name_directories,sep="\t",
               header=T,  #***
               stringsAsFactors=FALSE)
nfolders <- nrow(listDirs)

for (k in c(4:nfolders)){
   # get folder and summary file name
   tmpin <- listDirs[k,1]   #*****
   name_summary_file <- listDirs[k,2]   #*****
   fout <- listDirs[k,3]   #*****

   # get object file names
   tmpf <- list.files(path = tmpin, pattern = "object_results", 
           full.names = TRUE)
   tmps <- list.files(path = tmpin, pattern = "object_results", 
           full.names = FALSE)

   # get summary file
   fin_summary_file <- paste(tmpin,"/",name_summary_file,sep="")
   data_summary <- read.csv(fin_summary_file,stringsAsFactors=FALSE,header=T)
 
   # extract meta data from file names
   tmps_cond1 <- t(sapply(tmps,function(x) 
                   unlist(strsplit(x,split="\\."))[1:2]))


   # extract slide ID: different format from first two batches and others
   tmps_slideID <- matrix(NA,nrow(tmps_cond1),2)
   for (i in c(1:nrow(tmps_cond1))){
      tmp_str <- unlist(strsplit(tmps_cond1[i,1],split=" ")) 
      if (length(tmp_str)>2){
         tmps_slideID[i,] <- tmp_str[2:3]   #*****
      }else{
         tmps_slideID[i,] <- tmp_str
      }
   }
   rownames(tmps_slideID) <- NULL

   tmps_jobID <- sapply(tmps_cond1[,2],function(x)
                      substr(x,5,11))
   names(tmps_jobID) <- NULL

   # extract ROI from each file
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
         tmp4sum <- ifelse(is_sum_avai==TRUE,name_summary_file,"NA")
         tmp_info <- cbind(rep(tmpf[j],n_roi),
                        rep(tmps_slideID[j,1],n_roi),
                        rep(tmps_slideID[j,2],n_roi),
                        rep(tmps_jobID[j],n_roi),
                        names(roi_info),as.numeric(roi_info),
                        is_sum_avai,tmp4sum)
         listFiles <-  rbind(listFiles,tmp_info)
      }else{
         stop("No column associated with roi_column_interest found in image")

     }
   }

   colnames(listFiles) <- c("originalFileName","PathologicalID","BlockID",
                         "jobID","originalROI","cellNumber","IsSummaryAvai","SumFileName")

   write.table(listFiles,fout,sep="\t",quote=F,
             row.names=F,col.names=T)
}


# ************************************************************************
# Panel 2: Inspect level 3 data for batches: 12/20/23, 2/23/2024
# Multiple Summary data file --> single file/ROI --> skip this column
# ************************************************************************
# Run on R version 4.3.3
library(dplyr)

workFolder = c("D:/MIF/Human_SU2C/Ranalysis/")   #*****
setwd(workFolder)

# --------------------------batch 12/20/23: >1 summary file
inputParentFolder <- paste("E:/",
                        "MIF_UCLA/HALO_out/Panel2_SU2C/",
                         "Human Panel 2 Re-do of 1 case 12-20-23 version 1.2/",sep="")

fout <- c("listHaloOut_panel2_batch20231223.txt")  #****

tmpin <- inputParentFolder   #*****

# --------------------------batch 02/23/24: along with panel 1 and >1 summary file
inputParentFolder <- paste("E:/",
                        "MIF_UCLA/HALO_out//",
                        "20240202 Adeno Cases Panel 1+2 few slides Halo analysis/Panel2/",sep="")

fout <- c("listHaloOut_panel2_batch022324.txt")  #****

tmpin <- inputParentFolder   #*****

# ------------------------------
# get file names
tmpf <- list.files(path = tmpin, pattern = "object_results.csv", 
            full.names = TRUE)
tmps <- list.files(path = tmpin, pattern = "object_results.csv", 
            full.names = FALSE)

# extract meta data from file names
tmps_cond1 <- t(sapply(tmps,function(x) 
                unlist(strsplit(x,split="\\."))[1:2]))

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
 
# extract ROI from each file
nfiles <- length(tmpf)
listFiles <- NULL
roi_column_interest <- c("Analysis.Region")
 
for (j in c(1:nfiles)){
   tmp_data_in <- read.csv(tmpf[j],stringsAsFactors=FALSE,header=T)
   image_colnames <- colnames(tmp_data_in)
   if (is.element(roi_column_interest,image_colnames)){
      roi_info <- table(tmp_data_in[,roi_column_interest])
      n_roi <- length(roi_info)
      tmp_info <- cbind(rep(tmpf[j],n_roi),
                        rep(tmps_slideID[j,1],n_roi),
                        rep(tmps_slideID[j,2],n_roi),
                        rep(tmps_jobID[j],n_roi),
                        names(roi_info),as.numeric(roi_info))
      listFiles <-  rbind(listFiles,tmp_info)
 
   } else {
      stop("roi_column_interest not found in image")
   }
}

colnames(listFiles) <- c("originalFileName","PathologicalID","BlockID",
                         "jobID","originalROI","cellNumber")

write.table(listFiles,fout,sep="\t",quote=F,
              row.names=F,col.names=T)

# ************************************************************************
# Panel 2: Inspect level 3 data for batche Set 3 and GGO:tag is "object_Data"
# Multiple Summary data file --> single file/ROI --> skip this column
# ************************************************************************
library(dplyr)

workFolder = c("D:/MIF/Human_SU2C/Ranalysis/")   #*****
setwd(workFolder)

# -----------------------  PCGA2 Set 3 ---------------
inputParentFolder <- paste("E:/",
                        "MIF_UCLA/HALO_out/Panel2_SU2C/",
                        "Human Panel 2 Set 3 v1.2/ObjectData/",sep="")

name_summary_file <- c("Total_Summary_Results.csv")
tag_object <- c("object_Data.csv")

fout <- c("listHaloOut_panel2_batchSet3.txt")  #****

# -----------------------  GGO cohort - different headers ---------------
inputParentFolder <- paste("E:/",
                        "MIF_UCLA/HALO_out/Panel2_GGO/",
                        "GGO/",sep="")

name_summary_file <- c("Total_Summary_Results.csv")
tag_object <- c("object_Data.csv")

fout <- c("listHaloOut_panel2_batchGGO.txt")  #****


# ------------------------------
tmpin <- inputParentFolder   #*****

# get file names
tmpf <- list.files(path = tmpin, pattern = tag_object, 
           full.names = TRUE)
tmps <- list.files(path = tmpin, pattern = tag_object, 
           full.names = FALSE)

# get summary file
fin_summary_file <- paste(inputParentFolder,name_summary_file,sep="")
data_summary <- read.csv(fin_summary_file,stringsAsFactors=FALSE,header=T)
##colnames(data_summary)

# extract meta data from file names
tmps_cond1 <- t(sapply(tmps,function(x) 
                   unlist(strsplit(x,split="\\."))[1:2]))


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

# extract ROI from each file
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
      tmp4sum <- ifelse(is_sum_avai==TRUE,name_summary_file,"NA")
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

write.table(listFiles,fout,sep="\t",quote=F,
             row.names=F,col.names=T)

# **************************************************************************
# Panel 2: data generated after December 2024 has wrong column order for 
# the summary data --> concatenate and arrange columns for the summary data in 
# 20240202 Adeno Cases Panel 1+2 few slides Halo analysis/Panel2
# manually edit header --> Total_Summary_Results_newColOrder2.csv
# 11/11/2024
# ****************************************************************************
# Run on R version 4.3.3
library(dplyr)

workFolder = c("D:/MIF/Human_SU2C/Ranalysis/")   #*****
setwd(workFolder)

# --------------------------batch 02/23/24: along with panel 1 and >1 summary file
inputParentFolder <- paste("E:/",
                        "MIF_UCLA/HALO_out/",
                        "20240202 Adeno Cases Panel 1+2 few slides Halo analysis/Panel2/",sep="")

fout <- paste(inputParentFolder,
         "Total_Summary_Results_newColOrder.csv",sep="")  #****

tmpin <- inputParentFolder   #*****

# ---------------
# get list of summary files
list4sum <-  list.files(path = tmpin, pattern = "Total_Summary_Results.csv", 
           full.names = FALSE)

nfiles <- length(list4sum)

tmpSumData <- NULL
for (i in c(1:nfiles)){
   fin_summary_file <- paste(tmpin,list4sum[i],sep="")
   data_summary <- read.csv(fin_summary_file,stringsAsFactors=FALSE,header=T)

   if (i==1){
     sum_headers <- colnames(data_summary)
     sum_nocols <- ncol(data_summary)
     tmpSumData <- bind_rows(tmpSumData,data_summary)
   }else{
     isMatch <- all(colnames(data_summary)==sum_headers)
     if (isMatch==TRUE){
        tmpSumData <- bind_rows(tmpSumData,data_summary)
     }else{
        stop("columns in summary files do not match")
     }
   }
}

newSumData <- tmpSumData[,c(1:6,8:158,7)]
write.csv(newSumData,fout,quote=FALSE,row.names=F)

# *************************************************************************************
# Regenerate level 3 data for Panel 2
# 1. The nonGGO and GGO object and summary data files have different headers
# 2. The summary data of nonGGO after December 2024 has "Analysis Input" column moved
# from end the the fist 10 column --> redo summary data with name tag "_newColOrd"
# 3. Set 4 batch was splitted by correct and wrong IDs (i.e. wrong patientID and lesionID order
# 4. The summary data of wrong ID lesions modified in 2022 was corrupted 
# --> use the original file and tagged as "wrongID"
# 5. all summary file must be located at the same folder of their object files
# 6. Using csv files for headers of object and summary (kept special characters)
# 11/11/2024
# re-write csv as "UTF-8" format - 11/29/24
# **************************************************************************************
# Run on R version 4.3.3
library(dplyr)

workFolder = c("D:/MIF/Human_SU2C/Ranalysis/")   #*****
setwd(workFolder)

# -------------------------
# Step 0: Specify the input and output file/folder

# *****Below is for the GGO data
outputParentFolder <- c("D:/MIF/Panel2_SU2C4Sub_110924/")   #********
fout_subfolder <- c("GGO")   #******
fin_conversion <- c("MIF_panel2hu_GGOBatch_fileConversion_110924.txt")   #*****
fout_conv <- c("mif_P2_GGOBatch_oldnew_QuantifiedFile.txt")

suffix_obj <- c("_panel2_v1.2_object_results.csv")   #*********phenotype version different bet. MCL and non-MCL for Panel1, not P2
suffix_summary <- c("_panel2_v1.2_Summary_Results.csv")  #******

# ---------------------------
# ****set4 wrong ID

outputParentFolder <- c("D:/MIF/Panel2_SU2C4Sub_110924/")   #********
fout_subfolder <- c("nonGGO")   #******
fin_conversion <- c("MIF_panel2hu_pcgaSet4wrongID_fileConversion_110924.txt")   #*****
fout_conv <- c("mif_P2_pcgaSet4wrongID_oldnew_QuantifiedFile.txt")

suffix_obj <- c("_panel2_v1.2_object_results.csv")   #*********phenotype version different bet. MCL and non-MCL for Panel1, not P2
suffix_summary <- c("_panel2_v1.2_Summary_Results.csv")  #******

# ---------------------------
# ****nonSet4

outputParentFolder <- c("D:/MIF/Panel2_SU2C4Sub_110924/")   #********
fout_subfolder <- c("nonGGO")   #******
fin_conversion <- c("MIF_panel2hu_pcgaBatchNonSet4_fileConversion_110924.txt")   #*****
fout_conv <- c("mif_P2_pcgaBatchNonSet4_oldnew_QuantifiedFile.txt")

suffix_obj <- c("_panel2_v1.2_object_results.csv")   #*********phenotype version different bet. MCL and non-MCL for Panel1, not P2
suffix_summary <- c("_panel2_v1.2_Summary_Results.csv")  #******

# ---------------------------
# ****set4 correct ID

outputParentFolder <- c("D:/MIF/Panel2_SU2C4Sub_110924/")   #********
fout_subfolder <- c("nonGGO")   #******
fin_conversion <- c("MIF_panel2hu_pcgaSet4_fileConversion_110924.txt")   #*****
fout_conv <- c("mif_P2_pcgaSet4_oldnew_QuantifiedFile.txt")

suffix_obj <- c("_panel2_v1.2_object_results.csv")   #*********phenotype version different bet. MCL and non-MCL for Panel1, not P2
suffix_summary <- c("_panel2_v1.2_Summary_Results.csv")  #******


# --------------------------------------
# Step 1: import and check existing files
# Step 1.0: import list of object files and their information
sampleInfo <- read.delim(fin_conversion,sep="\t",header=T,
         stringsAsFactors=FALSE)
sampleInfo$JobIDnew <- sapply(sampleInfo$JobID, function(x) substr(x,4,nchar(x)),simplify=TRUE)

# ***for regenaration UTF-8 on 11/29/24 on laptop --> samsung HD
sampleInfo$oldFileName <- sampleInfo$fileName
sampleInfo$fileName <- sub("E:/","D:/",sampleInfo$oldFileName)

# Step 1.1: check existence of object data/quantified data
nfiles <- nrow(sampleInfo)
Ftab <- data.frame(fileName=sampleInfo[,1], is_exist=rep(NA,nfiles))
for (i in c(1:nfiles)){
   Ftab$is_exist[i] <- file.exists(sampleInfo[i,1])
}

# Step 1.2: check existence of summary data
tmp_sumfile <- paste(dirname(sampleInfo$fileName),"/",sampleInfo$OriginalSummaryFile,sep="")
sampleInfo$OriginalSummaryFile_fullPath <- tmp_sumfile
list4sum <- sampleInfo$OriginalSummaryFile_fullPath[!duplicated(sampleInfo$OriginalSummaryFile_fullPath)]
nsum <- length(list4sum)
Fsum <- data.frame(fileName=list4sum, is_exist=rep(NA,nsum))
for (i in c(1:nsum)){
   Fsum$is_exist[i] <- file.exists(list4sum[i])
}

# Step 1.3: inspect summary data for the ORIGINAL ROIs and job IDs
Fsum_roi <- NULL
for (i in c(1:nsum)){
   tmp_sum <- read.csv(list4sum[i],stringsAsFactors=FALSE,header=T)
   tmp <- cbind(rep(list4sum[i],nrow(tmp_sum)),tmp_sum[,c("Analysis.Region","Job.Id")])
   Fsum_roi <- rbind(Fsum_roi,tmp)
}
colnames(Fsum_roi) <- c("fileName","Analysis.Region","Job.Id")

# Step 1.4: Check overlap between summary files and object_results
nrow(Fsum_roi)
nrow(sampleInfo)
table(duplicated(Fsum_roi$Analysis.Region))

table(is.element(Fsum_roi$Analysis.Region,sampleInfo$Original_HTAN_ROI_ID))
table(is.element(sampleInfo$Original_HTAN_ROI_ID,Fsum_roi$Analysis.Region))


table(is.element(Fsum_roi$Job.Id,sampleInfo$JobIDnew))
table(is.element(sampleInfo$JobIDnew,Fsum_roi$Job.Id))
Fsum_roi[!is.element(Fsum_roi$Job.Id,sampleInfo$JobIDnew),]

#****STOP TO EVALUE THE OVERLAP ***


# Step 2: Upload header files (summary and object)and prepare for level 3 regeneration
# Step 2.1: get headers for ggo
# ggo header files
fin_header4obj_ggo <-  paste("D:/MIF/Human_SU2C/Ranalysis/",
                 "object_header_panel2_ggo_111124.csv",sep="") #*****new file with a new column defined keeping variables
ggo.obj_headers <-  read.csv(fin_header4obj_ggo,header=T,stringsAsFactors=FALSE)

fin_header4sum_ggo <-  paste("D:/MIF/Human_SU2C/Ranalysis/",
                 "summary_header_panel2_ggo_111124.csv",sep="") #*****new file with a new column defined keeping variables
ggo.sum_headers <-  read.csv(fin_header4sum_ggo,header=T,stringsAsFactors=FALSE)


# Step 2.2: get headers for pcga
fin_header4obj_pcga <-  paste("D:/MIF/Human_SU2C/Ranalysis/",
                 "object_header_panel2_pcga_111124.csv",sep="") #*****new file with a new column defined keeping variables
pcga.obj_headers <-  read.csv(fin_header4obj_pcga,header=T,stringsAsFactors=FALSE)


fin_header4sum_pcga <-  paste("D:/MIF/Human_SU2C/Ranalysis/",
                 "summary_header_panel2_pcga_111124.csv",sep="") #*****new file with a new column defined keeping variables
pcga.sum_headers <-  read.csv(fin_header4sum_pcga,header=T,stringsAsFactors=FALSE)

# Step 3: Rewrite the object_result files and change names
nfiles <- nrow(sampleInfo)
Ftab_new <- data.frame(fileName_Old=sampleInfo$fileName, fileName_New=rep(NA,nfiles),
                       summary_new=rep(NA,nfiles))

for (i in c(1:nfiles)){
   fin_tmp <- sampleInfo$fileName[i]
   fin_sum <- sampleInfo$OriginalSummaryFile_fullPath[i]

   sampleID_new <- sampleInfo$PCGA2_BiospecimenID[i]
   sampleID_old <- sampleInfo$Original_HTAN_ROI_ID[i]
   sampleID_jobid <- sampleInfo$JobIDnew[i]
   sampleID_local <- paste(sampleInfo$localPatID[i],"_",sampleInfo$BlockID[i],
                       ".tif",sep="")  ##slide_ID/parentalID

   fout_tmp <- paste(outputParentFolder,fout_subfolder,"/",sampleID_new,suffix_obj,sep="")
   fout_sum <- paste(outputParentFolder,fout_subfolder,"/",sampleID_new,suffix_summary,sep="")
   Ftab_new$fileName_New[i] <- fout_tmp
   Ftab_new$summary_new[i] <- fout_sum

   # get object data
   data_in <- read.csv(fin_tmp,stringsAsFactors=FALSE,header=T)
   tmp_selR <- which(data_in$Analysis.Region==sampleID_old)

   data_out <- data_in[tmp_selR,]
   data_out$Analysis.Region <- sampleID_new
   data_out[,1] <- sampleID_local

   # get summary data: time consumed but save
   tmp_sum <- read.csv(fin_sum ,stringsAsFactors=FALSE,header=T)
   isMatchID <- tmp_sum$Analysis.Region==sampleID_old & tmp_sum$Job.Id==sampleID_jobid
   data_sum <- tmp_sum[which(isMatchID==TRUE),]
   data_sum[,1] <- sampleID_local   #"Image Location"
   data_sum[,2] <- sampleID_local   #"Image Tag"
   data_sum$Analysis.Region <- sampleID_new

   # get header and modify classification class for ggo
   if (grepl("GGO",sampleInfo$cohort[i])==TRUE){
    
      tmp_obj_headers <- ggo.obj_headers$header
      keptCol_obj <- which(ggo.obj_headers$isKept==1) #****  
      data_out$Classifier.Label <- c("Tissue")   #***origincal labeled as "Tumor"

      tmp_sum_headers <- ggo.sum_headers$new_header  ##**replace old eader (col 1) by new header(col 2)
      keptCol_sum <- which(ggo.sum_headers$isKept==1)  #****
      orderKeptCol_sum <- as.numeric(ggo.sum_headers$ggo_new_order[keptCol_sum])  #switching Glass Area and Tumor Area
     
   }else{
      tmp_obj_headers <- pcga.obj_headers$header
      keptCol_obj <- which(pcga.obj_headers$isKept==1) #****
     
      tmp_sum_headers <- pcga.sum_headers$new_header
      keptCol_sum <- which(pcga.sum_headers$isKept==1) #***
      orderKeptCol_sum <- as.numeric(pcga.sum_headers$pcga_new_order[keptCol_sum])  #switching Glass Area and Tumor Area
   
   }
   
    colnames(data_out) <- tmp_obj_headers
    new_data_out <- data_out[,keptCol_obj]

    colnames(data_sum) <- tmp_sum_headers
    data_sum <- data_sum[,keptCol_sum]
    new_data_sum <- data_sum[,orderKeptCol_sum]    

   ## On November 11, 2024
   ##write.csv(new_data_out,fout_tmp,quote = FALSE,row.names=F)
   ##write.csv(new_data_sum,fout_sum,quote = FALSE,row.names=F)

   ## On 11/29/24
   write.csv(new_data_out,fout_tmp,quote = FALSE,row.names=F, fileEncoding="UTF-8")
   write.csv(new_data_sum,fout_sum,quote = FALSE,row.names=F, fileEncoding="UTF-8")
}

write.table(Ftab_new,fout_conv,sep="\t",
            quote=F,row.names=F)

# **************************************************************************
# Double check combination markers
# See MIF_Panel1_regnerateLevel3dara_080323.R for checking ROI area between 
# two panels
# 1. Clincal data : The final "xxxxx_sampleInfo_xxx.txt" files are combination between 
#     "xxx_oldnew_QuantifiedFile.txt" : dictionary/file names and 
#     "xxxx_fileConversion_xxx" : original files with sample information
# 2. Channel info: copy from "K:/MIF_UCLA/Panel2_SU2C/Panel2_Human_ChannelInformation.csv"
#     to local inputFolder
# There are 6 ROIs shared same HALO object files with others --> double check
# clinical files for replication of biospecimen IDs
# 12/23/24
# **************************************************************************
# Run on R version 4.3.3
library(SingleCellExperiment) 
library(reshape2)
library(dplyr)
library(apcluster)
library(scater)
library(dbscan)   #for identify)cell_communities
library(pheatmap)

source(paste("D:/PharmG/","Rsamples/SPIAT_modified_LT_010422.R",sep=""))

workFolder = c("D:/MIF/Human_SU2C/Ranalysis/")   #*****
setwd(workFolder)

# Step 0: Setting workplace and input folders
# Step 0.1: check if input files and those in sampleInfo are matched
inputParentFolder <- c("D:/MIF/Panel2_SU2C4Sub_110924/")  #****
subFolders <- c("nonGGO","GGO")

fin_clin <- c("MIF_panel2hu_GGOpcga_sampleInfo_110924.txt")
clin.data <- read.delim(fin_clin,header=T,
               stringsAsFactors=FALSE)

tmp_sc_files <- NULL
for (subinput in subFolders){
  tmpin <- paste(inputParentFolder,subinput,sep="")
  tmpf <- list.files(path = tmpin, pattern = "object_results.csv", 
           full.names = TRUE)
  tmp_sc_files <- c(tmp_sc_files,tmpf)
}

sum(clin.data$fileName_New==tmp_sc_files)

# Step 0.2: setting input, using files listed in fin_clin
# Step 0.2.1: setting chanel infor and phenotype
inputParentFolder <- c("D:/MIF/Panel2_SU2C4Sub_110924/")  #****
subFolders <- c("nonGGO","GGO")

channelInfo_fpath <- paste(inputParentFolder,
        "Panel2_Human_ChannelInformation.csv",sep="")   #****update here

# Step 0.2.2
fin_clin <- c("MIF_panel2hu_GGOpcga_sampleInfo_110924.txt")
clin.data <- read.delim(fin_clin,header=T,
               stringsAsFactors=FALSE)

# Step 1. import single-cell 
# Step 1.1: import and check all positive combination
allPhenotypes_mks <- NULL
nsamples <- nrow(clin.data)
for (i in c(1:nsamples)){   #**
  image_fpath <- clin.data$fileName_New[i]     #****
  sce <- format_halo_to_sce_4DL(image_fpath,
                      channelInfo_fpath,
                      haloPhenotype_fpath=NULL)

  # remove duplicated rows in some files
  tmpdup <- duplicated(colnames(sce))
  if (sum(tmpdup)>0){
     sce <- sce[,!tmpdup]
  }
  list_phenotype <- table(colData(sce)$Phenotype)
  mat_phenotype <- as.data.frame(list_phenotype)
  colnames(mat_phenotype) <- c("phenotype",clin.data$PCGA2_BiospecimenID[i]) ##**incorrect here
  rownames(mat_phenotype) <- NULL
 
  if (i==1){
     allPhenotypes_mks <- mat_phenotype
  }else{
     tmp_old <- allPhenotypes_mks
     allPhenotypes_mks <- full_join(tmp_old,mat_phenotype,
                 by = c("phenotype"="phenotype"))
  }
}

# create heatmap summary

tmp2 <- t(t(allPhenotypes_mks[,-1])/colSums(allPhenotypes_mks[,-1],na.rm=T))
rownames(tmp2) <- allPhenotypes_mks[,1]
tmp2 <- rbind(tmp2,colSums(allPhenotypes_mks[,-1],na.rm=T))

tmp_out <- cbind(c(as.character(allPhenotypes_mks[,1]),"totalCells"),
                 tmp2)
colnames(tmp_out) <- c("markers",colnames(tmp2))

write.table(tmp_out,"panel2_possibleCombination_122624.txt",sep="\t",
   quote=F,row.names=F)


tmp2[is.na(tmp2)] <- 0
pheatmap(tmp2)

# **************************************************************************
# Import data and define phenotypes. Input are
# 1. Clincal data : The final "xxxxx_sampleInfo_xxx.txt" files are combination between 
#     "xxx_oldnew_QuantifiedFile.txt" : dictionary/file names and 
#     "xxxx_fileConversion_xxx" : original files with sample information
# 2. Channel info: copy from "K:/MIF_UCLA/Panel2_SU2C/Panel2_Human_ChannelInformation.csv"
#     to local inputFolder
# 3. cell type dictionary: using "Markers_CT_conversion_122324.txt" whihc is almost
#     identical to "Markers_CT_conversion_010422.txt", except col 2 which
# "none: is replaced by "others"
# Note: Panel2 files dos NOT include Tissue Classification --> Using
#       format_halo_to_sce_4DL_DAPIfilter, specify $Phenotype as combnation markers
# 12/23/24
# **************************************************************************
# Run on R version 4.3.3
library(SingleCellExperiment) 
library(reshape2)
library(dplyr)
library(apcluster)
library(scater)
library(dbscan)   #for identify)cell_communities
library(pheatmap)
library(tibble)   # for rownames_to_column


source(paste("D:/PharmG/","Rsamples/SPIAT_modified_LT_010422.R",sep=""))

workFolder = c("D:/MIF/Human_SU2C/Ranalysis/")   #*****
setwd(workFolder)

# Step 0: Setting workplace and input folders
inputParentFolder <- c("D:/MIF/Panel2_SU2C4Sub_110924/")  #****
subFolders <- c("nonGGO","GGO")

# *******specify select cols for sammry files*******************
# ****see MIF_Panel1_regeneratedLevel3data_080324************
sumFile.selCol <- c(1:9,54,49,94)
sumFile.ColID <- c("Image.Location","Image.Tag","Algorithm.Name","Job.Id",
                       "Analysis.Region","Classified.Area","Glass.Area",
                       "Tissue.Area","Total.Cells","Tissue.Total.Cells",
                       "Avg.Cell.Area","Tissue.Avg.Cell.Area")

# import channel info
channelInfo_fpath <- paste(inputParentFolder,
        "Panel2_Human_ChannelInformation.csv",sep="")   #****update here

# phenotype definition:***** check script for final cols used in cell type definition ******
# ***there are two options --> cell_type have two options: tier 1 and tier 2
fin_dictionary  <- paste(inputParentFolder,
         "Markers_CT_conversion_122324.txt",sep="")   #*******
cellType_def <- read.delim(fin_dictionary,sep="\t",header=T,
               stringsAsFactors=FALSE)  

# import samples/ROIs info
fin_clin <- c("MIF_panel2hu_GGOpcga_sampleInfo_110924.txt")   #*******
clin.data <- read.delim(fin_clin,header=T,
               stringsAsFactors=FALSE)

# Step 1. import single-cell and summary data (1 ROI/file)
# Step 1.1: import and check all positive combination
selVar <- c("Phenotype")   #***update here
nct <- nrow(cellType_def)
stage_markers <- c("GranzymeB","FoxP3","Ki67")   #******
#stage_markers <- c("GranzymeB","FoxP3","Ki67","CD3","CD8","PanCK","DAPI")   #******


allPhenotypes <- NULL
allSummary <- NULL
nsamples <- nrow(clin.data)
for (i in c(1:nsamples)){   #**
  # import summary data
  fin4sum <- clin.data$summary_new[i]  #***check col names of clin.data
  tmpsum <- read.csv(fin4sum,header=T,stringsAsFactors=FALSE)
   
  if (i>2){
    colnames(tmpsum) <- colnames(allSummary)
  }
  allSummary <- bind_rows(allSummary,tmpsum)

  # import summary data
  image_fpath <- clin.data$fileName_New[i]     #****

  sce <- format_halo_to_sce_4DL_DAPIfilter(image_fpath,
                      channelInfo_fpath,
                      haloPhenotype_fpath=NULL)

  # remove duplicated rows in some files
  tmpdup <- duplicated(colnames(sce))
  if (sum(tmpdup)>0){
     sce <- sce[,!tmpdup]
  }

  # ****define cell_type based on column defined in the dictionary file***********
  # *****tier 1
  tmp_old <- colData(sce)[,selVar]
  tmp_new <- tmp_old
  for (jj in c(1:nct)){
     tmp_new[tmp_old==cellType_def[jj,1]] <- cellType_def[jj,2]   #****
  }
  colData(sce)$pheno_tier1 <- tmp_new   #****

  # *****tier 2
  tmp_old <- colData(sce)[,selVar]
  tmp_new <- tmp_old
  for (jj in c(1:nct)){
    tmp_new[tmp_old==cellType_def[jj,1]] <- cellType_def[jj,3]   #*****
  }
  colData(sce)$pheno_tier2 <- tmp_new  #****

  tmp_meta <- data.frame(colData(sce)) %>% rownames_to_column("Cell.ID")
  tmp_expr <- data.frame(t(assay(sce)[stage_markers,]))

  # add on 01/10/22 to skip stacking files with multiple regions (set 3) multiple times
  # include import expression of Gzmb, Foxp3 and Ki67
  if (i==1){ 
     allPhenotypes <- rbind(allPhenotypes,cbind(tmp_meta,tmp_expr))
  }else if (i>1 & image_fpath!=clin.data$fileName_New[(i-1)]){
     allPhenotypes <- rbind(allPhenotypes,cbind(tmp_meta,tmp_expr))
  }
}
rownames(allPhenotypes) <- NULL

tmp <- allSummary[,sumFile.selCol]
colnames(tmp) <- sumFile.ColID

allSummary_data <- tmp

# check duplicates
sum(duplicated(allPhenotypes$Cell.ID))

# -----------------------------
## store data
channelInfo <- read.csv(channelInfo_fpath,stringsAsFactors=FALSE,header=T)
keptL <- c("allPhenotypes","allSummary_data","clin.data","cellType_def",
           "channelInfo")
rm(list=setdiff(ls(),keptL))

# using format_halo_to_sce_4DL_DAPIfilter function: 3904168
fout <- c("Panel2_GGOpcga_colDataExpr_DAPIfilter_122324.rds")   #****update here
saveRDS(allPhenotypes,file=fout)

# keep intensity of all markers
fout <- c("Panel2_GGOpcga_colDataExpr_DAPIfilter_122324QC.rds")   #****update here
saveRDS(allPhenotypes,file=fout)

# save workplace
fout <- c("Panel2_GGOpcga_workplace_DAPIfilter_122324.RData")   #****update here
save.image(file="Panel2_GGOpcga_workplace_DAPIfilter_122324.RData")

fout <- c("Panel2_GGOpcga_summary_DAPIfilter_122324.rds")   #****update here
saveRDS(allSummary_data,file=fout)

# *****************************************************************************
# Check intensity threshold of each marker
# 1/2/25
# *****************************************************************************
# Run in 4.3.3
library(SingleCellExperiment) 
library(reshape2)
library(dplyr)
library(apcluster)
library(scater)
library(dbscan)   #for identify)cell_communities
library(pheatmap)
library(tibble)   # for rownames_to_column
library(lmerTest)
library(ggridges)

source(paste("D:/PharmG/","Rsamples/SPIAT_modified_LT_010422.R",sep=""))

workFolder = c("D:/MIF/Human_SU2C/Ranalysis/")   #*****
setwd(workFolder)

# Step 0: load QC data
fin_objData <- c("Panel2_GGOpcga_colDataExpr_DAPIfilter_122324QC.rds")   #****update here
allPhenotypes_all <- readRDS(fin_objData)

# import clinical data
fin_clin <- c("MIF_panel2hu_GGOpcga_sampleInfo_110924.txt")   #*******
clin.data <- read.delim(fin_clin,header=T,
               stringsAsFactors=FALSE)

clin.data$His3 <- as.character(clin.data$tissueType)
clin.data$His3 <- ifelse(is.element(clin.data$His3,c("MIA","ADC")),c("MIA_ADC"),clin.data$His3)
clin.data$His3 <- factor(clin.data$His3,
         levels=c("Normal","AAH","AIS","MIA_ADC"))
clin.data$tissueType <- factor(clin.data$tissueType,
         levels=c("Normal","AAH","AIS","MIA","ADC"))

# Step 1: inspect intensity data
# Step 1.1: check if marker intensity are correlated (i.e. batch effect)
stage_markers <- c("GranzymeB","FoxP3","Ki67","CD3","CD8","PanCK","DAPI")
selCol <- c("ROI","CellArea","Phenotype") 


# whole data = Tissue for Panel 2
tmpdata <- allPhenotypes_all[,c(selCol,stage_markers)]   #******
tmp_mean <- aggregate(tmpdata[,stage_markers],by=list(tmpdata$ROI),mean)
tmp_med <- aggregate(tmpdata[,stage_markers],by=list(tmpdata$ROI),median)

cor(tmp_mean[,-1])
cor(tmp_med[,-1])

# check cutoff
tmpdata <- allPhenotypes_all[,c(selCol,stage_markers)]   #******
tmp_thres <- matrix(NA,length(stage_markers),2)
for (i in c(1:length(stage_markers))){
  tmp_pos <- grepl(stage_markers[i],tmpdata$Phenotype)
  tmp_thres[i,] <- as.numeric(aggregate(tmpdata[,stage_markers[i]],
                         by=list(tmp_pos),min)[,2])
}
colnames(tmp_thres) <- c("Negative","Positive")
rownames(tmp_thres) <- stage_markers

write.table(cbind(tmp_thres,rownames(tmp_thres)),
      "tmp.txt",sep="\t",quote=F,row.names=F)

# ***********************************************************************************# Run on R version 4.3.3
# Calculate the density
# 12/26/24
# ****************************************************************************
# Run in 4.3.3
library(SingleCellExperiment) 
library(reshape2)
library(dplyr)
library(apcluster)
library(scater)
library(dbscan)   #for identify)cell_communities
library(pheatmap)
library(tibble)   # for rownames_to_column
library(lmerTest)
library(ggridges)

source(paste("D:/PharmG/","Rsamples/SPIAT_modified_LT_010422.R",sep=""))

workFolder = c("D:/MIF/Human_SU2C/Ranalysis/")   #*****
setwd(workFolder)

# Step 1: import data and clean up
# Step 1.1: import data
fin_p2 <- c("Panel2_GGOpcga_workplace_DAPIfilter_122324.RData")   #****update here
load(fin_p2)

table(allPhenotypes$Phenotype)

# Step 1.2: clean up single-cell objects
# Filter out "Removed" cells (non-biological combinations)
allPhenotypes_filter <- allPhenotypes
tmpi <- which(allPhenotypes_filter$pheno_tier1=="PanCK,GranzymeB,FoxP3,Ki67,CD8")
allPhenotypes_filter$pheno_tier1[tmpi] <- c("Removed")

keptR <- allPhenotypes_filter$pheno_tier1 != c("Removed")
allPhenotypes_filter <- allPhenotypes_filter[keptR,]

# Filter out "doublets" by cell area
area_thres <- 100 #*****
keptR <- allPhenotypes_filter$CellArea<area_thres  #****
allPhenotypes_filter <- allPhenotypes_filter[keptR,]

# Step 1.3: clean up summary data
# Re-count cells in Tissue.Total.Cells
tmp_count <- table(allPhenotypes_filter$ROI)

tmp_count_new <- left_join(allSummary_data[,c("Analysis.Region","Total.Cells")],
             data.frame(ROI=names(tmp_count),Tissue.Total.Cells=as.numeric(tmp_count)),
             by=c("Analysis.Region"="ROI"))

# check ROIs with large % of cells removed by area conditions
tmp <- tmp_count_new$Tissue.Total.Cells/tmp_count_new$Total.Cells
summary(tmp)
tmp_bad_idx <- which(tmp<0.9)
cbind(tmp_count_new[tmp_bad_idx,],tmp[tmp_bad_idx])
allSummary_data[tmp_bad_idx,]

# Update summary data
allSummary_data_filter <- allSummary_data
allSummary_data_filter$Tissue.Total.Cells <- tmp_count_new$Tissue.Total.Cells

# Step 2: calculate density
selTier = c("org_tier1") #***including sub-lineages
tmp_cellNum <- table(allPhenotypes_filter$ROI,
                     allPhenotypes_filter$pheno_tier1)

tmp_cellNum_short <- data.frame(rbind(tmp_cellNum)) %>%
     rownames_to_column("ROI")

# merge data: data4cellNum in wide-form 
selCol4Sum <- c("Analysis.Region","Tissue.Area","Tissue.Total.Cells")
data4cellNum <- clin.data[,-c(1:3)] %>%
      left_join(allSummary_data_filter[,selCol4Sum],
           by=c("PCGA2_BiospecimenID"="Analysis.Region")) %>%
      left_join(tmp_cellNum_short,by=c("PCGA2_BiospecimenID"="ROI"))


data4cellNum$Tissue.Total.Cells<- NULL

# -----
fout <- c("Panel2_GGOpcga_cellNumber4BU_122324.rds") 
saveRDS(data4cellNum,file=fout)
# ----


# merge data: data4density in long-form
data4density <- as.data.frame(tmp_cellNum)   #**long/stacked form
colnames(data4density) <- c("PCGA2_BiospecimenID","cellType","cellNumber")

data4density <- data4density %>%
      left_join(allSummary_data_filter[,c("Analysis.Region","Tissue.Area")],
            by=c("PCGA2_BiospecimenID"="Analysis.Region")) %>%
      left_join(clin.data[,-c(1:3)],join_by(PCGA2_BiospecimenID))

data4density$density <- data4density$cellNumber/data4density$Tissue.Area
data4density$tissueType <- factor(data4density$tissueType,
                          levels=c("Normal","AAH","AIS","MIA","ADC"),
                          labels=c("Normal","AAH","AIS","MIA","ADC"))  #*****
data4density$His3 <- as.character(data4density$tissueType)
data4density$His3 <- ifelse(is.element(data4density$His3,c("MIA","ADC")),c("MIA_ADC"),data4density$His3)
data4density$His3 <- factor(data4density$His3,
                          levels=c("Normal","AAH","AIS","MIA_ADC"),
                          labels=c("Normal","AAH","AIS","MIA_ADC"))  #*****

# Step 3: Statistical test: using delta
ct_list <- names(table(data4density$cellType))
noCT <- length(ct_list)
Ftest <- matrix(NA,nrow=noCT,12)
for (i in c(1:noCT)){
   selCT <- ct_list[i]
   tmpsel <- is.element(data4density$cellType,selCT)
   tmpdata <- subset(data4density,tmpsel)
 
   ## kruskal test all groups
   Ftest[i,1] <- kruskal.test(density~His3,tmpdata)$p.value

   ## mixed linear all groups
   fm1 = lmer(density~His3+(1|PCGA2_PatID),tmpdata)    ##mixed linear          
   Ftest[i,2] <- anova(fm1)$"Pr(>F)"
  
   ## pairwise: AAH vs. Nor
   tmpdata2 <- subset(tmpdata,His3 %in% c("Normal","AAH"))
   tmpdata2$His3 <- droplevels(tmpdata2$His3)
   Ftest[i,3] <- wilcox.test(density~His3,tmpdata2)$p.value
   fm1 = lmer(density~His3+(1|PCGA2_PatID),tmpdata2)    ##mixed linear          
   Ftest[i,4] <- anova(fm1)$"Pr(>F)"

   ## pairwise: AIS vs. Nor
   tmpdata2 <- subset(tmpdata,His3 %in% c("Normal","AIS"))
   tmpdata2$His3 <- droplevels(tmpdata2$His3)
   Ftest[i,5] <- wilcox.test(density~His3,tmpdata2)$p.value
   fm1 = lmer(density~His3+(1|PCGA2_PatID),tmpdata2)    ##mixed linear          
   Ftest[i,6] <- anova(fm1)$"Pr(>F)"

   ## pairwise: Invasive vs. Nor
   tmpdata2 <- subset(tmpdata,His3 %in% c("Normal","MIA_ADC"))
   tmpdata2$His3 <- droplevels(tmpdata2$His3)
   Ftest[i,7] <- wilcox.test(density~His3,tmpdata2)$p.value
   fm1 = lmer(density~His3+(1|PCGA2_PatID),tmpdata2)    ##mixed linear          
   Ftest[i,8] <- anova(fm1)$"Pr(>F)"
  
   ## pairwise: AIS vs. AAH
   tmpdata2 <- subset(tmpdata,His3 %in% c("AAH","AIS"))
   tmpdata2$His3 <- droplevels(tmpdata2$His3)
   Ftest[i,9] <- wilcox.test(density~His3,tmpdata2)$p.value
   fm1 = lmer(density~His3+(1|PCGA2_PatID),tmpdata2)    ##mixed linear          
   Ftest[i,10] <- anova(fm1)$"Pr(>F)"

   ## pairwise: MIA_ADC vs. AAH
   tmpdata2 <- subset(tmpdata,His3 %in% c("AAH","MIA_ADC"))
   tmpdata2$His3 <- droplevels(tmpdata2$His3)
   Ftest[i,11] <- wilcox.test(density~His3,tmpdata2)$p.value
   fm1 = lmer(density~His3+(1|PCGA2_PatID),tmpdata2)    ##mixed linear          
   Ftest[i,12] <- anova(fm1)$"Pr(>F)"

}
colnames(Ftest) <- c("all_kruskal_pv","all_mixed_pv",
                     "AAHvNor_wilcox_pv","AAHvNor_mixed_pv",
                     "AISvNor_wilcox_pv","AISvNor_mixed_pv",
                     "InvvNor_wilcox_pv","InvvNor_mixed_pv",
                     "AISvAAH_wilcox_pv","AISvAAH_mixed_pv",
                     "InvvAAH_wilcox_pv","InvvAAH_mixed_pv")                    
rownames(Ftest) <- ct_list

# ---
Ftab_out <- as.data.frame(Ftest) %>% rownames_to_column("cellType")
write.table(Ftab_out,"tmp.txt",sep="\t",quote=F,row.names=F)


## Step 4: plotting
colCode4Hist <- c("#999999","#009E73","#56B4E9","#D55E00")
names(colCode4Hist) <- c("Normal","AAH","AIS","MIA_ADC")

tmpSeclType <- c("CD8")   #*****
tmpSeclType_label <- c("CD8")   #*****

ymax <- max(data4density$density[is.element(data4density$cellType,tmpSeclType)])
p1 <- ggplot(subset(data4density,is.element(data4density$cellType,tmpSeclType)),
     aes(x=His3,y=density,fill=His3))+
     geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
     geom_jitter(position=position_jitter(width=.1, height=0))+
     ylim(y=c(0,signif(ymax,2)))+
     ggtitle(paste(tmpSeclType_label,":subsolid cohort",sep=""))+
     scale_fill_manual(values=colCode4Hist)+theme_bw()

