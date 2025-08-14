# *******************************************************************************
# Check HALO output for panel1 LUAD using powershell
# 12/28/23
cd 'K:\MIF_UCLA\Panel1_SU2C\20210607 Adeno Cases Panel 1\Halo Analysis'

Get-ChildItem -Path .\ -Recurse -File | 
    Get-FileHash -Algorithm MD5 | 
    Select-Object  @{n='Hash';e={$_.Hash.ToLower()}},@{N='Path';E={$_.Path | Resolve-Path -Relative}} |
    Export-Csv -Path 'K:\MIF_UCLA\Panel1_SU2C\20210607 Adeno Cases Panel 1\panel1_20210607_md5checksum.csv' -UseCulture -NoTypeInformation

# *******************************************************************************
# Note of data history
# 1. First data batch in 12/2023 has tumor, stromal, glass calssification
# including:
# a. 20210607 Adeno Case Panel 1: MCL, 51 ROIs
# b. 20231211 Adeno Case Panel 1 Halo Analysis: PCGA2 with 93 ROIs
# c. 20231220 Adeno Case Panel 1 - 1 case from initial Panel 1 analysis: missing case from MCL, 2 ROIs
# --> high % of glass --> requantified and received in 1/2024 as 
# "20240109 Adeno Cases Panel 1 Halo Analysis"
# --> remove a-c folders to tobedel
# 2. Re-quantified data "20240109 Adeno Cases Panel 1 Halo Analysis" with 95 ROIs
# Inspect quantified data in 
# a. PCGA: "20240109 Adeno Cases Panel 1 Halo Analysis"
# b. MCL: "20240131 Adeno initial Cases Panel 1 redo of classifier"
# 3. Redo/add missing slides
# 4. MCL: "20241105 Adeno Cases Panel 1 redo with new classifier ObjectData":
#    increase cells in tissue area
# *********************************************************************************

# **********************************************************************************
# Regenerate meta data for batches: 
# a. PCGA: "20240109 Adeno Cases Panel 1 Halo Analysis"
# b. MCL: "20240131 Adeno initial Cases Panel 1 redo of classifier"
# c. batch 02/28/2024,
# d. batch 05/03/2024
# Other small batches have its own script because of (1) file name format and 
# (2) muilple sumary files
# 10/25/24
# redo MCL batch 11/05/24 with new tissue classifier --> run on desktop
# **********************************************************************************
# Run on R version 4.3.3
library(dplyr)

workFolder = c("D:/MIF/Human_SU2C/Ranalysis/")   #*****
setwd(workFolder)

# ----********PCGA2 ******************************
inputParentFolder <- paste("E:/",
                        "MIF_UCLA/HALO_out/Panel1_SU2C/",
                         "20240109 Adeno Cases Panel 1 Halo Analysis/",sep="")
name_summary_file <- c("20240109 Adeno Cases Pane 1 export.csv")

fout <- c("listHaloOut_panel1_batch20240109.txt")  #****

# ----********MCL - Camelia to re-export data on 1/30/24 ******************************
inputParentFolder <- paste("E:/",
                        "MIF_UCLA/HALO_out/Panel1_SU2C/",
                        "20240131 Adeno initial Cases Panel 1 redo of classifier/",sep="")
name_summary_file <- c("Total_Summary_Results.csv")
fout <- c("listHaloOut_panel1_mcl20240131.txt")  #****

#----********batch 02/28/2024 ********************************************
inputParentFolder <- paste("E:/",
                        "MIF_UCLA/HALO_out/Panel1_SU2C/",
                        "20240228 Adeno Cases Panel 1 a few cases/",sep="")
name_summary_file <- c("Total_Summary_Results.csv")
fout <- c("listHaloOut_panel1_batch20240224.txt")  #****

#----********batch 05/03/2024 ********************************************
inputParentFolder <- paste("E:/",
                        "MIF_UCLA/HALO_out/Panel1_SU2C/",
                        "20240503 Adeno Cases Panel 1 a couple of cases/",sep="")
name_summary_file <- c("Total_Summary_Results.csv")
fout <- c("listHaloOut_panel1_batch20240503.txt")  #****

# ***********bach 11/05/24 *************************************************
inputParentFolder <- paste("D:/",
                    "MIF_UCLA/HALO_out/Panel1/",
                    "20241104 Adeno Cases Panel 1 analysis redo with different classifier/",
                    "MCL/",sep="")

name_summary_file <- c("20241105 Adeno Cases Panel 1 redo with new classifier Total_Summary_Results.csv")
fout <- c("listHaloOut_panel1_mcl20241104.txt")  #****

# ------------------------------
tmpin <- inputParentFolder   #*****

# get file names
tmpf <- list.files(path = tmpin, pattern = "object_results.csv", 
           full.names = TRUE)
tmps <- list.files(path = tmpin, pattern = "object_results.csv", 
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

# ************************************************************************
# Panel 1: Inspect level 3 data for batch 02/23/2024, received with Panel2
# Multiple Summary data file --> single file/ROI --> skip this column
# ************************************************************************
# Run on R version 4.3.3
library(dplyr)

workFolder = c("D:/MIF/Human_SU2C/Ranalysis/")   #*****
setwd(workFolder)

# list the original file names of qunatified data
inputParentFolder <- paste("E:/MIF_UCLA/HALO_out/",
       "20240202 Adeno Cases Panel 1+2 few slides Halo analysis/Panel1/",sep="")
fout <- c("listHaloOut_panel1_batch022324.txt")  #***
 
tmpin <- inputParentFolder   #*****

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


# ***************************************************************
# Regenerate level 3 data for Panel 1 MCL data
# Checked headers of obj and summary files and created text files which specify 
# kept and removed columns. 
# If the headers are not the same among batches --> specify kept and removed
# column in each batch --> rearranged column to a unit order
# MCL tissue classification: using v2 algorithm --> tagged as v2
# 10/25/24
# 11/04/24: regnerate the MCL data using v3 algorithm for tissue classification
# --> tagged as v3 on laptop and samsung HD
# 11/29/24: correct output for write.csv with fileEncoding="UTF-8"
# ***************************************************************
# Run on R version 4.3.3
library(dplyr)

workFolder = c("D:/MIF/Human_SU2C/Ranalysis/")   #*****
setwd(workFolder)


# *****Below is for the new quantified by v2 in 1/2024 --> rerun on 10/25/24
# Step 0: Specify the input and output file/folder

#-----
# ***fileConversion includes FULL PATH to the input
outputParentFolder <- c("D:/MIF/Panel1_SU2C4Sub_031224/")   #********
fin_conversion <- c("MIF_panel1hu_MCLBatch013124_fileConversion_031424.txt")   #*****
fin_sumFile <- paste("E:/",
                "MIF_UCLA/HALO_out/Panel1_SU2C/",
                "20240131 Adeno initial Cases Panel 1 redo of classifier/",
                "Total_Summary_Results.csv",sep="")
fout_conv <- c("mif_P1_MCLBatch013124_oldnew_QuantifiedFile.txt")

suffix_obj <- c("_panel1_v2_object_results.csv")   #*********phenotype version different bet. MCL and non-MCL for Panel1, not P2
suffix_summary <- c("_panel1_v2_Summary_Results.csv")  #******

#-------------
# singe MCL file in batch 2/23/24
# ***fileConversion includes FULL PATH to the input
outputParentFolder <- c("D:/MIF/Panel1_SU2C4Sub_031224/")   #********
fin_conversion <- c("MIF_panel1hu_MCLBatch022324_fileConversion_042324.txt")   #*****
fin_sumFile <- paste("E:/MIF_UCLA/HALO_out/",
                    "20240202 Adeno Cases Panel 1+2 few slides Halo analysis/Panel1/",
                    "SSW-17-13377 A10 Panel 1 Total_Summary_Results.csv",sep="") #*********
fout_conv <- c("mif_P1_MCLBatch022324_oldnew_QuantifiedFile.txt")

suffix_obj <- c("_panel1_v2_object_results.csv")   #*********phenotype version different bet. MCL and non-MCL for Panel1, not P2
suffix_summary <- c("_panel1_v2_Summary_Results.csv")  #******

# -------
# ***fileConversion includes FULL PATH to the input for 11/04/24 batch
outputParentFolder <- c("D:/MIF/Panel1_SU2C4Sub_110424/")   #********
fin_conversion <- c("MIF_panel1hu_MCL_fileConversion_110424.txt")   #*****ROIs in 01/31/24 & 02/23/24
fin_sumFile <- paste("D:/",
                    "MIF_UCLA/HALO_out/Panel1/",
                    "20241104 Adeno Cases Panel 1 analysis redo with different classifier/",
                    "MCL/",
                    "20241105 Adeno Cases Panel 1 redo with new classifier Total_Summary_Results.csv",sep="")

fout_conv <- c("mif_P1_MCLBatch110424_oldnew_QuantifiedFile.txt")

suffix_obj <- c("_panel1_v3_object_results.csv")   #*********phenotype version different bet. MCL and non-MCL for Panel1, not P2
suffix_summary <- c("_panel1_v3_Summary_Results.csv")  #******
#-------------

# Step 1: Check existing files
# Step 1.0: import samples/ROIs
sampleInfo <- read.delim(fin_conversion,sep="\t",header=T,
         stringsAsFactors=FALSE)

# Step 1.1: check quantified data
nfiles <- nrow(sampleInfo)
Ftab <- data.frame(fileName=sampleInfo[,1], is_exist=rep(NA,nfiles))
for (i in c(1:nfiles)){
   Ftab$is_exist[i] <- file.exists(sampleInfo[i,1])
}

# Step 1.2: check summary data
tmp_sum <- read.csv(fin_sumFile,stringsAsFactors=FALSE,header=T)
dim(tmp_sum)

# Step 1.3: Check overlap between summary files and object_results
table(is.element(tmp_sum$Analysis.Region,sampleInfo$Original_HTAN_ROI_ID))
table(is.element(sampleInfo$Original_HTAN_ROI_ID,tmp_sum$Analysis.Region))

#****STOP IF EXISTING FILE(S) IS NOT OVERLAPPED AMONG TWO FILE CATEGORIES ********

# Step 2: Rewrite the object_result files and change names
# Step 2.1: get headers
fin_header_pcgamcl <-  paste("D:/MIF/Human_SU2C/Ranalysis/",
                 "object_header_panel1_pcgamcl_031424.txt",sep="") #*****new file with a new column defined keeping variables
pcgamcl.obj_headers <-  read.delim(fin_header_pcgamcl,sep="\t",header=T,  #***
               stringsAsFactors=FALSE)

fin_header_summary <-  paste("D:/MIF/Human_SU2C/Ranalysis/",
                 "summary_header_panel1_pcgamcl_031424.txt",sep="") #*****new file with a new column defined keeping variables
pcgamcl.sum_headers <-  read.delim(fin_header_summary,sep="\t",header=T,  #***
               stringsAsFactors=FALSE)

# Step 2.2: reload summary file
mcl.allSummary <- read.csv(fin_sumFile,stringsAsFactors=FALSE,header=T)
ncol(mcl.allSummary)
nrow(pcgamcl.sum_headers)

nfiles <- nrow(sampleInfo)
Ftab_new <- data.frame(fileName_Old=sampleInfo[,1], fileName_New=rep(NA,nfiles),
                       summary_new=rep(NA,nfiles))

for (i in c(1:nfiles)){
   fin_tmp <- sampleInfo$fileName[i]

   sampleID_new <- sampleInfo$PCGA2_BiospecimenID[i]
   sampleID_old <- sampleInfo$Original_HTAN_ROI_ID[i]
   sampleID_local <- paste(sampleInfo$localPatID[i],"_",sampleInfo$BlockID[i],
                       ".tif",sep="")  ##slide_ID/parentalID                     

   fout_subfolder <- c("MCL")   #******
   fout_tmp <- paste(outputParentFolder,fout_subfolder,"/",sampleID_new,
                      suffix_obj,sep="")       #**diff from P2
   fout_sum <- paste(outputParentFolder,fout_subfolder,"/",sampleID_new,
                       suffix_summary,sep="")   #**diff from P2

   Ftab_new$fileName_New[i] <- fout_tmp
   Ftab_new$summary_new[i] <- fout_sum

   # get object data
   data_in <- read.csv(fin_tmp,stringsAsFactors=FALSE,header=T)
   tmp_selR <- which(data_in$Analysis.Region==sampleID_old)

   data_out <- data_in[tmp_selR,]
   data_out$Analysis.Region <- sampleID_new
   data_out[,1] <- sampleID_local

   # get summary data, information of kept variables
   tmp_obj_headers <- pcgamcl.obj_headers[,1]
   keptCol_obj <- which(pcgamcl.obj_headers[,2]==1) #****

   tmp_sum_headers <- pcgamcl.sum_headers[,1]
   keptCol_sum <- which(pcgamcl.sum_headers[,2]==1) #***

   data_sum <- mcl.allSummary[which(mcl.allSummary$Analysis.Region==sampleID_old),]

   data_sum[,1] <- sampleID_local   #"Image Location"
   data_sum[,2] <- sampleID_local   #"Image Tag"
   data_sum$Analysis.Region <- sampleID_new

   colnames(data_out) <- tmp_obj_headers
   colnames(data_sum) <- tmp_sum_headers

   ## ***remove column associated with HALO-defined phenotype
   new_data_out <- data_out[,keptCol_obj]
   new_data_sum <- data_sum[,keptCol_sum]

   write.csv(new_data_out,fout_tmp,quote = FALSE,row.names=F)
   write.csv(new_data_sum,fout_sum,quote = FALSE,row.names=F)
}

write.table(Ftab_new,fout_conv,sep="\t",quote=F,row.names=F)


# ***************************************************************
# Regenerate level 3 data for Panel 1 non-MCL data
# Checked headers of obj and summary files and created text files which specify 
# kept and removed columns. 
# If the headers are not the same among batches --> specify kept and removed
# column in each batch --> rearranged column to a unit order
# version 3: before May 2024: have Glass > tissue
# version 4: re-classification of glass vs tisues of 5 ROIs (patID = 30006 and 30034
# 10/23/24
# ***************************************************************
# Run on R version 4.3.3
library(dplyr)

workFolder = c("D:/MIF/Human_SU2C/Ranalysis/")   #*****
setwd(workFolder)
# -------------------------
# Step 0: Specify the input and output file/folder

# *****Below is for the new quantified by v3 in 1/2024
outputParentFolder <- c("D:/MIF/Panel1_SU2C4Sub_031224/")   #********
fin_conversion <- c("MIF_panel1hu_nonMCLv3Batches_fileConversion_102324.txt")   #*****
fout_conv <- c("mif_P1_nonMCLv3Batches_oldnew_QuantifiedFile.txt")

suffix_obj <- c("_panel1_v3_object_results.csv")   #*********phenotype version different bet. MCL and non-MCL for Panel1, not P2
suffix_summary <- c("_panel1_v3_Summary_Results.csv")  #******

# ---------------------
# *****Below is for the new quantified by v4 in 5/2024

# ***fileConversion includes FULL PATH to the input
outputParentFolder <- c("D:/MIF/Panel1_SU2C4Sub_031224/")   #********
fin_conversion <- c("MIF_panel1hu_nonMCLv4Batch_fileConversion_102324.txt") 
fout_conv <- c("mif_P1_nonMCLv4Batches_oldnew_QuantifiedFile.txt")

suffix_obj <- c("_panel1_v4_object_results.csv")   #*********phenotype version different bet. MCL and non-MCL for Panel1, not P2
suffix_summary <- c("_panel1_v4_Summary_Results.csv")  #******  #*****

# -----------------
sampleInfo <- read.delim(fin_conversion,sep="\t",header=T,
         stringsAsFactors=FALSE)

# Step 1: Check existing files
# Step 1.1: check quantified data
nfiles <- nrow(sampleInfo)
Ftab <- data.frame(fileName=sampleInfo[,1], is_exist=rep(NA,nfiles))
for (i in c(1:nfiles)){
   Ftab$is_exist[i] <- file.exists(sampleInfo[i,1])
}

# Step 1.2: check summary data
tmp_sumfile <- paste(dirname(sampleInfo$fileName),"/",sampleInfo$OriginalSummaryFile,sep="")
sampleInfo$OriginalSummaryFile_fullPath=tmp_sumfile
list4sum <- sampleInfo$OriginalSummaryFile_fullPath[!duplicated(sampleInfo$OriginalSummaryFile_fullPath)]
nsum <- length(list4sum)
Fsum <- data.frame(fileName=list4sum, is_exist=rep(NA,nsum))
for (i in c(1:nsum)){
   Fsum$is_exist[i] <- file.exists(list4sum[i])
}

Fsum_roi <- NULL
for (i in c(1:nsum)){
   tmp_sum <- read.csv(list4sum[i],stringsAsFactors=FALSE,header=T)
   tmp <- cbind(rep(list4sum[i],nrow(tmp_sum)),tmp_sum[,c("Analysis.Region","Job.Id")])
   Fsum_roi <- rbind(Fsum_roi,tmp)
}
colnames(Fsum_roi) <- c("fileName","Analysis.Region","Job.Id")

# Step 1.3: Check overlap between summary files and object_results
nrow(Fsum_roi)
nrow(sampleInfo)
table(duplicated(Fsum_roi$Analysis.Region))
table(is.element(Fsum_roi$Analysis.Region,sampleInfo$Original_HTAN_ROI_ID))
table(is.element(sampleInfo$Original_HTAN_ROI_ID,Fsum_roi$Analysis.Region))

## sampleInfo does not list quantified data  for GGO7, and requalified data (GGO6 and 80024)
Fsum_roi[!is.element(Fsum_roi$Analysis.Region,sampleInfo$Original_HTAN_ROI_ID),]

#****STOP IF EXISTING FILE(S) IS NOT OVERLAPPED AMONG TWO FILE CATEGORIES ********

# Step 3: Rewrite the object_result files and change names
# Step 3.1: get headers
fin_header_pcgamcl <-  paste("D:/MIF/Human_SU2C/Ranalysis/",
                 "object_header_panel1_pcgamcl_031424.txt",sep="") #*****new file with a new column defined keeping variables
pcgamcl.obj_headers <-  read.delim(fin_header_pcgamcl,sep="\t",header=T,  #***
               stringsAsFactors=FALSE)

fin_header_summary <-  paste("D:/MIF/Human_SU2C/Ranalysis/",
                 "summary_header_panel1_pcgamcl_031424.txt",sep="") #*****new file with a new column defined keeping variables
pcgamcl.sum_headers <-  read.delim(fin_header_summary,sep="\t",header=T,  #***
               stringsAsFactors=FALSE)

# Step 3.2: reload summary file
list4sum <- sampleInfo$OriginalSummaryFile_fullPath[!duplicated(sampleInfo$OriginalSummaryFile_fullPath)]
nsum <- length(list4sum)
non_mcl.allSummary  <- NULL
for (i in c(1:nsum)){
   tmp_sum <- read.csv(list4sum[i],stringsAsFactors=FALSE,header=T)
   non_mcl.allSummary <- bind_rows(non_mcl.allSummary,tmp_sum)
}

ncol(non_mcl.allSummary)
nrow(pcgamcl.sum_headers)

nfiles <- nrow(sampleInfo)
Ftab_new <- data.frame(fileName_Old=sampleInfo[,1], fileName_New=rep(NA,nfiles),
                       summary_new=rep(NA,nfiles))

for (i in c(1:nfiles)){
   fin_tmp <- sampleInfo$fileName[i]

   sampleID_new <- sampleInfo$PCGA2_BiospecimenID[i]
   sampleID_old <- sampleInfo$Original_HTAN_ROI_ID[i]
   sampleID_local <- paste(sampleInfo$localPatID[i],"_",sampleInfo$BlockID[i],
                       ".tif",sep="")  ##slide_ID/parentalID                     

   fout_subfolder <- c("nonMCL")   #******
   fout_tmp <- paste(outputParentFolder,fout_subfolder,"/",sampleID_new,
                      suffix_obj,sep="")       #**diff from P2
   fout_sum <- paste(outputParentFolder,fout_subfolder,"/",sampleID_new,
                       suffix_summary,sep="")   #**diff from P2

   Ftab_new$fileName_New[i] <- fout_tmp
   Ftab_new$summary_new[i] <- fout_sum

   # get object data
   data_in <- read.csv(fin_tmp,stringsAsFactors=FALSE,header=T)
   tmp_selR <- which(data_in$Analysis.Region==sampleID_old)

   data_out <- data_in[tmp_selR,]
   data_out$Analysis.Region <- sampleID_new
   data_out[,1] <- sampleID_local

   # get summary data, information of kept variables
   tmp_obj_headers <- pcgamcl.obj_headers[,1]
   keptCol_obj <- which(pcgamcl.obj_headers[,2]==1) #****

   tmp_sum_headers <- pcgamcl.sum_headers[,1]
   keptCol_sum <- which(pcgamcl.sum_headers[,2]==1) #***

   data_sum <- non_mcl.allSummary[which(non_mcl.allSummary$Analysis.Region==sampleID_old),]

   data_sum[,1] <- sampleID_local   #"Image Location"
   data_sum[,2] <- sampleID_local   #"Image Tag"
   data_sum$Analysis.Region <- sampleID_new

   colnames(data_out) <- tmp_obj_headers
   colnames(data_sum) <- tmp_sum_headers

   ## ***remove column associated with HALO-defined phenotype
   new_data_out <- data_out[,keptCol_obj]
   new_data_sum <- data_sum[,keptCol_sum]

   write.csv(new_data_out,fout_tmp,quote = FALSE,row.names=F)
   write.csv(new_data_sum,fout_sum,quote = FALSE,row.names=F)
}

write.table(Ftab_new,fout_conv,sep="\t",quote=F,row.names=F)


# ******************************************************************************************
# Rewrite level 3 with "UTF-8"
# 11/29/24: correct output for write.csv with fileEncoding="UTF-8"
# ****************************************************************************
# Run on R version 4.3.3
library(dplyr)

workFolder = c("D:/MIF/Human_SU2C/Ranalysis/")   #*****
setwd(workFolder)

# header files: object and summary
fin_header_pcgamcl <-  paste("D:/MIF/Human_SU2C/Ranalysis/",
                 "object_header_panel1_pcgamcl_110924.csv",sep="") #*****new file with a new column defined keeping variables
pcgamcl.obj_headers <-  read.csv(fin_header_pcgamcl,header=T,  #***
               stringsAsFactors=FALSE)

fin_header_summary <-  paste("D:/MIF/Human_SU2C/Ranalysis/",
                 "summary_header_panel1_pcgamcl_110924.csv",sep="") #*****new file with a new column defined keeping variables
pcgamcl.sum_headers <-  read.csv(fin_header_summary,header=T,  #***
               stringsAsFactors=FALSE)

# ----------------
fin_listFiles <- c("mif_P1_MCLBatch110424_oldnew_QuantifiedFile.txt")   #******

fin_listFiles <- c("mif_P1_nonMCLv3Batches_oldnew_QuantifiedFile.txt")   #******

fin_listFiles <- c("mif_P1_nonMCLv4Batches_oldnew_QuantifiedFile.txt")   #******

fin_listFiles <- c("mif_P1_MCLBatch013124_oldnew_QuantifiedFile.txt")   #******

fin_listFiles <- c("mif_P1_MCLBatch022324_oldnew_QuantifiedFile.txt")   #******

listFiles_p1 <- c("mif_P1_MCLBatch110424_oldnew_QuantifiedFile.txt",
                  "mif_P1_nonMCLv3Batches_oldnew_QuantifiedFile.txt",
                  "mif_P1_nonMCLv4Batches_oldnew_QuantifiedFile.txt")
                  

 tmps
[1] "mif_P1_MCLBatch013124_oldnew_QuantifiedFile.txt"  
[2] "mif_P1_MCLBatch022324_oldnew_QuantifiedFile.txt"  
[3] "mif_P1_MCLBatch110424_oldnew_QuantifiedFile.txt"  
[4] "mif_P1_nonMCLv3Batches_oldnew_QuantifiedFile.txt" 
[5] "mif_P1_nonMCLv4Batches_oldnew_QuantifiedFile.txt" 
[6] "mif_P2_GGOBatch_oldnew_QuantifiedFile.txt"        
[7] "mif_P2_pcgaBatchNonSet4_oldnew_QuantifiedFile.txt"
[8] "mif_P2_pcgaSet4_oldnew_QuantifiedFile.txt"        
[9] "mif_P2_pcgaSet4wrongID_oldnew_QuantifiedFile.txt"

# -----------------------------
listFiles  <- read.delim(fin_listFiles,sep="\t",header=T,
         stringsAsFactors=FALSE)


nfiles <- nrow(listFiles)
for (i in c(1:nfiles)){
   fin_obj <- listFiles[i,2]
   fin_sum <- listFiles[i,3]

   fout_obj_new <- paste(dirname(fin_obj),"_new/",basename(fin_obj),sep="")
   fout_sum_new <- paste(dirname(fin_sum),"_new/",basename(fin_sum),sep="")

   # get object data
   data_obj_in <- read.csv(fin_obj,stringsAsFactors=FALSE,header=T,
                 check.names = F)
   data_obj_out <- data_obj_in
   colnames(data_obj_out) <- pcgamcl.obj_headers[which(pcgamcl.obj_headers[,2]==1),1]

   # get sum data
   data_sum_in <- read.csv(fin_sum,stringsAsFactors=FALSE,header=T,
                 check.names = F)
   data_sum_out <- data_sum_in
   colnames(data_sum_out) <- pcgamcl.sum_headers[which(pcgamcl.sum_headers[,2]==1),1]
 
   write.csv(data_obj_out,fout_obj_new,quote = FALSE,row.names=F, fileEncoding="UTF-8")  #******
   write.csv(data_sum_out,fout_sum_new,quote = FALSE,row.names=F, fileEncoding="UTF-8")  #******
}


##tmp1 <- read.csv(fout_obj_new,header=T,stringsAsFactors=FALSE)
##tmp2 <- read.csv(fout_sum_new,header=T,stringsAsFactors=FALSE)

# *********************************************************************************
# Inspect ROI- tissue area of panel 1 and 2
# Concatenate panel 2 "xx_oldnew_QuantifiedFile.txt" files
# 11/29/24 - on laptop
# *********************************************************************************
# Run on R version 4.3.3
library(dplyr)

workFolder = c("D:/MIF/Human_SU2C/Ranalysis/")   #*****
setwd(workFolder)

# list all "oldnew_QuantifiedFile.txt"
inputFolder <- c("D:/MIF/Human_SU2C/Ranalysis/")
tag_suffix <- c("oldnew_QuantifiedFile.txt")
tmps <- list.files(path = inputFolder, pattern = tag_suffix, 
           full.names = FALSE)
tmps

# using MCL generated on 11/03/24
listFiles_p1 <- c("mif_P1_MCLBatch110424_oldnew_QuantifiedFile.txt",
                  "mif_P1_nonMCLv3Batches_oldnew_QuantifiedFile.txt",
                  "mif_P1_nonMCLv4Batches_oldnew_QuantifiedFile.txt")
                  

listFiles_p2 <- c("mif_P2_pcgaBatchNonSet4_oldnew_QuantifiedFile.txt",
               "mif_P2_pcgaSet4_oldnew_QuantifiedFile.txt",
               "mif_P2_pcgaSet4wrongID_oldnew_QuantifiedFile.txt",
               "mif_P2_GGOBatch_oldnew_QuantifiedFile.txt")

# panel 1
nfiles_p1 <- length(listFiles_p1)
all_files_p1 <- NULL
for (i in c(1:nfiles_p1)){
   tmp_data_in <- read.delim(listFiles_p1[i],sep="\t",header=T,
         stringsAsFactors=FALSE)
   all_files_p1 <- bind_rows(all_files_p1,tmp_data_in)
}

# panel 2
nfiles_p2 <- length(listFiles_p2)
all_files_p2 <- NULL
for (i in c(1:nfiles_p2)){
   tmp_data_in <- read.delim(listFiles_p2[i],sep="\t",header=T,
         stringsAsFactors=FALSE)
   all_files_p2 <- bind_rows(all_files_p2,tmp_data_in)
}

all_files_p2 <- tmpall  

# import sumary data
nrois_p1 <- nrow(all_files_p1)
allROI_sum_p1 <- NULL
for (i in c(1:nrois_p1)){
   tmp_data_in <- read.csv(all_files_p1[i,3],header=T,stringsAsFactors=FALSE)
   allROI_sum_p1 <- bind_rows(allROI_sum_p1,tmp_data_in)
}

tmp_data_in <- read.table(all_files_p1[i,3],sep=",",header=T,stringsAsFactors=FALSE,
                  fileEncoding = "Latin1",check.names = F)



# ******************************************************************************
# Inspect ROI- tissue area of panel 1 and 2
# using "check.names=F" option in read.csv when import all files witn symbols in the header
# 11/29/24 (on new Dell Desktop)
# ******************************************************************************
# Run on R version 4.3.3
library(dplyr)

workFolder = c("D:/MIF/Human_SU2C/Ranalysis/")   #*****
setwd(workFolder)

##panel1_inputFolder <- c("D:/MIF/Panel1_SU2C4Sub_031224/MCL_old/",
##                        "D:/MIF/Panel1_SU2C4Sub_031224/nonMCL_old/")

panel1_inputFolder <- c("D:/MIF/Panel1_SU2C4Sub_110424/MCL/",
                        "D:/MIF/Panel1_SU2C4Sub_031224/nonMCL/")

panel2_inputFolder <- c("D:/MIF/Panel2_SU2C4Sub_110924/nonGGO/",
                        "D:/MIF/Panel2_SU2C4Sub_110924/GGO/")

tag_suffix <- c("Summary_Results.csv")

# get all summary file names instead of using those listed 
# in "xxx_oldnew_QuantifiedFile.txt"
listFiles_p1  <- NULL
for (inputFolder in panel1_inputFolder){
   tmp_sum_list <- list.files(path = inputFolder, pattern = tag_suffix, 
           full.names = TRUE)
   listFiles_p1 <- c(listFiles_p1,tmp_sum_list)
}

listFiles_p2  <- NULL
for (inputFolder in panel2_inputFolder){
   tmp_sum_list <- list.files(path = inputFolder, pattern = tag_suffix, 
           full.names = TRUE)
   listFiles_p2 <- c(listFiles_p2,tmp_sum_list)
}

# import summary data
# panel 1
p1_summ_data <- NULL
for (i in c(1:length(listFiles_p1))){
   tmpsum <- read.csv(listFiles_p1[i],header=T,
               stringsAsFactors=FALSE)   ###check.names=F

   p1_summ_data <- bind_rows(p1_summ_data,tmpsum)
}
p1.selSumCol <- c(1:9,99,54,49,139,94)
colnames(p1_summ_data)[p1.selSumCol]

p1_summ <- p1_summ_data[,p1.selSumCol]
colnames(p1_summ) <- c("Image.Location","Image.Tag","Algorithm.Name","Job.Id",
                       "Analysis.Region","Classified.Area","Glass.Area",
                       "Tissue.Area","Total.Cells","Tissue.Total.Cells",
                       "Glass.Total.Cells","Avg.Cell.Area",
                       "Tissue.Avg.Cell.Area","Glass.Avg.Cell.Area")
p1_summ$area_T2total <- p1_summ$Tissue.Area/p1_summ$Classified.Area
p1_summ$cell_T2total <- p1_summ$Tissue.Total.Cells/p1_summ$Total.Cells

# panel 2
p2_summ_data <- NULL
p2.selSumCol <- c(1:9,54,49,94)   #no export for "Glass.Total.Cells","Glass.Avg.Cell.Area"
for (i in c(1:length(listFiles_p2))){
   tmpsum <- read.csv(listFiles_p2[i],header=T,
               stringsAsFactors=FALSE)

   p2_summ_data <- bind_rows(p2_summ_data,tmpsum)
}
p2.selSumCol <- c(1:9,54,49,94)   #no export for "Glass.Total.Cells","Glass.Avg.Cell.Area"
colnames(p2_summ_data)[p2.selSumCol]

p2_summ <- p2_summ_data[,p2.selSumCol]
colnames(p2_summ) <- c("Image.Location","Image.Tag","Algorithm.Name","Job.Id",
                       "Analysis.Region","Classified.Area","Glass.Area",
                       "Tissue.Area","Total.Cells","Tissue.Total.Cells",
                       "Avg.Cell.Area","Tissue.Avg.Cell.Area")
p2_summ$area_T2total <- p2_summ$Tissue.Area/p2_summ$Classified.Area
p2_summ$cell_T2total <- p2_summ$Tissue.Total.Cells/p2_summ$Total.Cells


pall_summ <- inner_join(p1_summ,p2_summ,
     by=join_by("Analysis.Region"))

# Explore data
plot(p1_summ$Tissue.Area,p1_summ$Glass.Area)
plot(p1_summ$Tissue.Avg.Cell.Area,p1_summ$Glass.Avg.Cell.Area)
plot(log10(p1_summ$Tissue.Total.Cells),log10(p1_summ$Glass.Total.Cells))
abline(a=0,b=1)

plot(pall_summ$area_T2total.x,pall_summ$area_T2total.y)
cor(pall_summ$area_T2total.x,pall_summ$area_T2total.y)
summary(lm(pall_summ$area_T2total.x~pall_summ$area_T2total.y))

plot(pall_summ$cell_T2total.x,pall_summ$cell_T2total.y)
plot(log10(pall_summ$Tissue.Total.Cells.x),log10(pall_summ$Tissue.Total.Cells.y))

# ******************** MCL data *******************************
# *****There are two classification associated with tissue (Region) and cell type (Phenotype)
# *****Old files do not have a column of "Analysis.Region"
# *****Ask Camelia to re-export data on 1/30/24
# *****Below is for the new quantified in 1/2024

inputParentFolder <- paste("E:/",
                          "MIF_UCLA/HALO_out/Panel1_SU2C/",
                           "20240131 Adeno initial Cases Panel 1 redo of classifier/",sep="")
fout <- c("listHaloOut_panel1_mcl20240131.txt")  #****

tmpin <- inputParentFolder   #*****

# get file names
tmpf <- list.files(path = tmpin, pattern = "object_results.csv", 
           full.names = TRUE)
tmps <- list.files(path = tmpin, pattern = "object_results.csv", 
           full.names = FALSE)

# extract meta data from file names
tmps_cond1 <- t(sapply(tmps,function(x) 
                   unlist(strsplit(x,split="\\."))[1:2]))
tmps_slideID <- t(sapply(tmps_cond1[,1],function(x)
                   unlist(strsplit(x,split=" "))[1:2]))
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
             row.names=F,col.names=F)

 
# **************************************************************************
# **************************************************************************
# Double check combination markers for panel 1 MCL and non MCL
# See MIF_Panel1_regnerateLevel3dara_122823.R for checking ROI area between 
# two panels
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


# ****MCL ************************
inputParentFolder <- c("D:/MIF/Panel1_SU2C4Sub_110424/")  #****
subFolders <- c("MCL")
channelInfo_fpath <- paste(inputParentFolder,
        "Panel1MCL_Human_ChannelInformation.csv",sep="")   #****update here

fin_clin <- c("MIF_panel1hu_MCL_sampleInfo_110424.txt")
fout_combine <- c("panel1MCL_possibleCombination_122624.txt")

# ****nonMCL ************************
inputParentFolder <- c("D:/MIF/Panel1_SU2C4Sub_031224/")  #****
subFolders <- c("nonMCL")
channelInfo_fpath <- paste(inputParentFolder,
        "Panel1nonMCL_Human_ChannelInformation.csv",sep="")   #****update here

fin_clin <- c("MIF_panel1hu_nonMCL_sampleInfo_102324.txt")
fout_combine <- c("panel1nonMCL_possibleCombination_122624.txt")

#---------------------------------
# Step 1: Check if input files and those in sampleInfo are matched
clin.data <- read.delim(fin_clin,header=T,
               stringsAsFactors=FALSE)


tmp_sc_files <- NULL
for (subinput in subFolders){
  tmpin <- paste(inputParentFolder,subinput,sep="")
  tmpf <- list.files(path = tmpin, pattern = "object_results.csv", 
           full.names = TRUE)
  tmp_sc_files <- c(tmp_sc_files,tmpf)
}

dim(clin.data)
length(tmp_sc_files)
nrow(clin.data)-sum(is.element(clin.data$fileName_New,tmp_sc_files))

# --------------------
# Step 2: check combination of markers
clin.data <- read.delim(fin_clin,header=T,
               stringsAsFactors=FALSE)
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

tmp2 <- t(t(allPhenotypes_mks[,-1])/colSums(allPhenotypes_mks[,-1],na.rm=T))
rownames(tmp2) <- allPhenotypes_mks[,1]
tmp2 <- rbind(tmp2,colSums(allPhenotypes_mks[,-1],na.rm=T))

tmp_out <- cbind(c(as.character(allPhenotypes_mks[,1]),"totalCells"),
                 tmp2)
colnames(tmp_out) <- c("markers",colnames(tmp2))

write.table(tmp_out,fout_combine,sep="\t",
   quote=F,row.names=F)

tmp2[is.na(tmp2)] <- 0
pheatmap(tmp2[-65,])


# -----
# Concatenate combination in MCL and nonMCL P1
fin_mcl <- c("panel1MCL_possibleCombination_122624.txt")
mcl_data <- read.delim(fin_mcl,sep="\t",header=T,
         stringsAsFactors=FALSE)
mcl_data$median <- apply(mcl_data[,-1],1,median,na.rm=T)
mcl_data$mean<- apply(mcl_data[,-1],1,mean,na.rm=T)

fin_nonmcl <- c("panel1nonMCL_possibleCombination_122624.txt")
nonmcl_data <- read.delim(fin_nonmcl,sep="\t",header=T,
         stringsAsFactors=FALSE)
nonmcl_data$median <- apply(nonmcl_data[,-1],1,median,na.rm=T)
nonmcl_data$mean<- apply(nonmcl_data[,-1],1,mean,na.rm=T)


tmp_mcl_header <- mcl_data$markers
new_mcl_com <- rep(NA,length(tmp_mcl_header))
for (i in c(1:length(tmp_mcl_header))){
  tmp_ele <- sort(unlist(strsplit(tmp_mcl_header[i],split=",")))
  new_mcl_com[i] <- paste0(tmp_ele,collapse=",")
}
mcl_data$marker_order <- new_mcl_com

tmp_nonmcl_header <- nonmcl_data$markers
new_nonmcl_com <- rep(NA,length(tmp_nonmcl_header))
for (i in c(1:length(tmp_nonmcl_header))){
  tmp_ele <- sort(unlist(strsplit(tmp_nonmcl_header[i],split=",")))
  new_nonmcl_com[i] <- paste0(tmp_ele,collapse=",")
}
nonmcl_data$marker_order <- new_nonmcl_com

length(new_nonmcl_com)  #64
length(new_mcl_com)   #64

tmp4mcl <- mcl_data[,c("marker_order","markers","median","mean")] 
colnames(tmp4mcl) <- paste("mcl_",colnames(tmp4mcl),sep="")

tmp4nonmcl <- nonmcl_data[,c("marker_order","markers","median","mean")] 
colnames(tmp4nonmcl) <- paste("nonmcl_",colnames(tmp4nonmcl),sep="")

tmpall <- full_join(tmp4mcl,tmp4nonmcl,
              by=c("mcl_marker_order"="nonmcl_marker_order"))
write.table(tmpall,"tmp.txt",sep="\t",quote=F,row.names=F)


# ******************************************************************************************
# Import data and define phenotypes. Input are
# 1. Clincal data : The final "xxxxx_sampleInfo_xxx.txt" files are combination between 
#     "xxx_oldnew_QuantifiedFile.txt" : dictionary/file names and 
#     "xxxx_fileConversion_xxx" : original files with sample information
# 2. Channel info: There are two versions for MCL and nonMCL cohort, copy from HALO output
#     to local inputFolder
# 3. cell type dictionary: using "Markers_CT_conversion_P1_122324.txt" 
# Note: Panel1 files DOES include Tissue Classification --> Using
#       format_haloTissueClass_to_sce_4DL_DAPIfilter, 
#       where specify $Phenotype as combnation markers (in the alphabet order)
# 12/23/24
# re-run with exporting all intensitive --> "QC.rds" 
# *******************************************************************************************
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

# Step 0: setting worksplace and input folders
# ****MCL ************************
inputParentFolder <- c("D:/MIF/Panel1_SU2C4Sub_110424/")  #****
subFolders <- c("MCL")
channelInfo_fpath <- paste(inputParentFolder,
        "Panel1MCL_Human_ChannelInformation.csv",sep="")   #****update here

fin_clin <- c("MIF_panel1hu_MCL_sampleInfo_110424.txt")
fout_objData <- c("Panel1_MCL_colDataExpr_DAPIfilter_122324.rds")
fout_sumData <- c("Panel1_MCL_summary_DAPIfilter_122324.rds")

# ****nonMCL ************************
inputParentFolder <- c("D:/MIF/Panel1_SU2C4Sub_031224/")  #****
subFolders <- c("nonMCL")
channelInfo_fpath <- paste(inputParentFolder,
        "Panel1nonMCL_Human_ChannelInformation.csv",sep="")   #****update here

fin_clin <- c("MIF_panel1hu_nonMCL_sampleInfo_102324.txt")
fout_objData <- c("Panel1_nonMCL_colDataExpr_DAPIfilter_122324.rds")
fout_sumData <- c("Panel1_nonMCL_summary_DAPIfilter_122324.rds")

#---------------------------------
# *******specify select cols for sammry files*******************
# ****see MIF_Panel1_regeneratedLevel3data_080324************
sumFile.selCol <- c(1:9,99,54,49,139,94)
sumFile.ColID <- c("Image.Location","Image.Tag","Algorithm.Name","Job.Id",
                       "Analysis.Region","Classified.Area","Glass.Area",
                       "Tissue.Area","Total.Cells","Tissue.Total.Cells",
                       "Glass.Total.Cells","Avg.Cell.Area",
                       "Tissue.Avg.Cell.Area","Glass.Avg.Cell.Area")

# phenotype definition:***** check script for final cols used in cell type definition ******
# ***there are two options --> cell_type have two options: tier 1 and tier 2
fin_dictionary  <- paste(inputParentFolder,
         "Markers_CT_conversion_P1_122324.txt",sep="")   #*******
cellType_def <- read.delim(fin_dictionary,sep="\t",header=T,
               stringsAsFactors=FALSE)  

# import samples/ROIs info
clin.data <- read.delim(fin_clin,header=T,
               stringsAsFactors=FALSE)

# Step 1. import single-cell and summary data (1 ROI/file)
# Step 1.1: import and check all positive combination
selVar <- c("Phenotype")   #***update here
nct <- nrow(cellType_def)
stage_markers <- c("PDL1","PD1")   #******
#stage_markers <- c("PDL1","PD1","CD3","CD8","CD68","PanCK","DAPI")   #******

allPhenotypes <- NULL
allSummary <- NULL
nsamples <- nrow(clin.data)
for (i in c(1:nsamples)){   #**
  # import summary data
  fin4sum <- clin.data$summary_file[i]  #***check col names of clin.data
  tmpsum <- read.csv(fin4sum,header=T,stringsAsFactors=FALSE)
   
  if (i>2){
    colnames(tmpsum) <- colnames(allSummary)
  }
  allSummary <- bind_rows(allSummary,tmpsum)

  # import summary data
  image_fpath <- clin.data$fileName_New[i]     #****

  sce <- format_haloTissueClass_to_sce_4DL_DAPIfilter(image_fpath,
                      channelInfo_fpath,
                      haloPhenotype_fpath=NULL)

  # remove duplicated rows in some files
  tmpdup <- duplicated(colnames(sce))
  if (sum(tmpdup)>0){
     sce <- sce[,!tmpdup]
  }

  # re-order list of markers due to different fluorophore-antibody pairs
  tmp_old_pheno <- colData(sce)[,selVar]
  tmp_old_list <- names(table(tmp_old_pheno))
  tmp_oldNew_list <- reorder_markerList(tmp_old_list)
  tmp_alphabet_phenotype <- tmp_old_pheno
  for (jj in c(1:nrow(tmp_oldNew_list))){
     tmp_alphabet_phenotype[tmp_old_pheno==tmp_oldNew_list[jj,1]] <- tmp_oldNew_list[jj,2]   #****
  }
  colData(sce)$PhenotypeOrdered <- tmp_alphabet_phenotype

  # ****define cell_type based on column defined in the dictionary file***********
  # *****lenient
  tmp_old <-  tmp_alphabet_phenotype
  tmp_new <- tmp_old
  for (jj in c(1:nct)){
     tmp_new[tmp_old==cellType_def[jj,1]] <- cellType_def[jj,4]   #****
  }
  colData(sce)$pheno_lenient <- tmp_new   #****

  # *****strict
  tmp_old <- tmp_alphabet_phenotype
  tmp_new <- tmp_old
  for (jj in c(1:nct)){
    tmp_new[tmp_old==cellType_def[jj,1]] <- cellType_def[jj,5]   #*****
  }
  colData(sce)$pheno_strict <- tmp_new  #****

  # *****tier 1 (major cell types)
  tmp_old <- tmp_alphabet_phenotype
  tmp_new <- tmp_old
  for (jj in c(1:nct)){
    tmp_new[tmp_old==cellType_def[jj,1]] <- cellType_def[jj,6]   #*****
  }
  colData(sce)$pheno_tier1 <- tmp_new  #****

  tmp_meta <- data.frame(colData(sce)) %>% rownames_to_column("Cell.ID")
  tmp_expr <- data.frame(t(assay(sce)[stage_markers,]))

  # add on 01/10/22 to skip stacking files with multiple regions (set 3) multiple times
  # include import expression of Gzmb, Foxp3 and Ki67
  if (i==1){ 
     allPhenotypes <- rbind(allPhenotypes,cbind(tmp_meta,tmp_expr))
  }else if (i>1 & image_fpath!=clin.data$fileName[(i-1)]){
     allPhenotypes <- rbind(allPhenotypes,cbind(tmp_meta,tmp_expr))
  }
}
rownames(allPhenotypes) <- NULL

tmp <- allSummary[,sumFile.selCol]
colnames(tmp) <- sumFile.ColID

allSummary_data <- tmp

# using format_halo_to_sce_4DL_DAPIfilter function: 3904168
saveRDS(allPhenotypes,file=fout_objData)
saveRDS(allSummary_data ,file=fout_sumData)


# ---------------
allPhenotypes_mcl <- allPhenotypes
allSummary_mcl <-allSummary_data

allPhenotypes_nonmcl <- allPhenotypes
allSummary_nonmcl <-allSummary_data


allPhenotypes_p1 <- bind_rows(allPhenotypes_mcl,allPhenotypes_nonmcl)
allSummary_p1 <- bind_rows(allSummary_mcl,allSummary_nonmcl)

fout_objData <- c("Panel1_all_colDataExpr_DAPIfilter_122324.rds")
fout_sumData <- c("Panel1_all_summary_DAPIfilter_122324.rds")
saveRDS(allPhenotypes_p1,file=fout_objData)
saveRDS(allSummary_p1 ,file=fout_sumData)

fin_clin_mcl <- c("MIF_panel1hu_MCL_sampleInfo_110424.txt")
fin_clin_nonMCL <- c("MIF_panel1hu_nonMCL_sampleInfo_102324.txt")
clin.data_mcl <- read.delim(fin_clin_mcl,header=T,stringsAsFactors=FALSE)
clin.data_nonMCL <- read.delim(fin_clin_nonMCL,header=T,stringsAsFactors=FALSE)

sum(colnames(clin.data_mcl)==colnames(clin.data_nonMCL))
dim(clin.data_mcl)
dim(clin.data_nonMCL)

clin.data_p1 <- bind_rows(clin.data_mcl,clin.data_nonMCL)
fout_clin <- c("MIF_panel1hu_all_sampleInfo_122324.txt")
write.table(clin.data_p1,fout_clin,sep="\t",quote=F,row.names=F)


saveRDS(allPhenotypes_p1,file=fout_objData)

# **************************************************************************
# Inspect intensity of panel 1 data
# Note: 
# 1. Intensities of markers are not correlated to each other for both MCL and nonMCL,
# 2. mean has higher correlation than median --> outlines is batch sensitive    
# 
# 12/30/24
# **************************************************************************
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
# load MCL data
fin_objData <- c("Panel1_MCL_colDataExpr_DAPIfilter_122324QC.rds")  #****
fin_sumData <- c("Panel1_MCL_summary_DAPIfilter_122324.rds") #*****

allPhenotypes_mcl <- readRDS(fin_objData)
allSummary_mcl <- readRDS(fin_sumData)

# load nonMCL data
fin_objData <- c("Panel1_nonMCL_colDataExpr_DAPIfilter_122324QC.rds")  #****
fin_sumData <- c("Panel1_nonMCL_summary_DAPIfilter_122324.rds")   #****

allPhenotypes_nonMCL <- readRDS(fin_objData)
allSummary_nonMCL <- readRDS(fin_sumData)

# import clinical data
fin_clin <- c("MIF_panel1hu_all_sampleInfo_122324.txt")
clin.data <- read.delim(fin_clin,sep="\t",header=T,stringsAsFactors=FALSE)  

clin.data$His3 <- as.character(clin.data$tissueType)
clin.data$His3 <- ifelse(is.element(clin.data$His3,c("MIA","ADC")),c("MIA_ADC"),clin.data$His3)
clin.data$His3 <- factor(clin.data$His3,
         levels=c("Normal","AAH","AIS","MIA_ADC"))
clin.data$tissueType <- factor(clin.data$tissueType,
         levels=c("Normal","AAH","AIS","MIA","ADC"))

colCode4Hist <- c("#999999","#009E73","#56B4E9","#0072B2","#D55E00")
names(colCode4Hist) <- c("Normal","AAH","AIS","MIA_ADC")
 

# Step 1: inspect intensity data
# Step 1.1: check if marker intensity are correlated (i.e. batch effect)
# --> found there is no correlation for either whole or "Tissue" cells
stage_markers <- c("PDL1","PD1","CD3","CD8","CD68","PanCK","DAPI")
selCol <- c("ROI","CellArea","ClassifierLabel","PhenotypeOrdered")

# whole data
tmpdata <- allPhenotypes_nonMCL[,c(selCol,stage_markers)]   #******
tmp_mean <- aggregate(tmpdata[,stage_markers],by=list(tmpdata$ROI),mean)
tmp_med <- aggregate(tmpdata[,stage_markers],by=list(tmpdata$ROI),median)

cor(tmp_mean[,-1])

# "Tissue" or "Glass" associated cells  
selTissue <- c("Tissue")  #*****      
tmpdata <- allPhenotypes_mcl[,c(selCol,stage_markers)] %>%   #******
           dplyr::filter(ClassifierLabel==selTissue)
tmp_mean <- aggregate(tmpdata[,stage_markers],by=list(tmpdata$ROI),mean)
tmp_med <- aggregate(tmpdata[,stage_markers],by=list(tmpdata$ROI),median)

cor(tmp_mean[,-1])
cor(tmp_med[,-1])

# check cutoff
selTissue <- c("Tissue")  #***** 
# ---
tmpdata <- allPhenotypes_mcl[,c(selCol,stage_markers)]   %>%   #******
           dplyr::filter(ClassifierLabel==selTissue)

tmpdata <- allPhenotypes_nonMCL[,c(selCol,stage_markers)] %>%   #******
           dplyr::filter(ClassifierLabel==selTissue)

# ---
tmp_thres <- matrix(NA,length(stage_markers),2)
for (i in c(1:length(stage_markers))){
  tmp_pos <- grepl(stage_markers[i],tmpdata$PhenotypeOrdered)
  tmp_thres[i,] <- as.numeric(aggregate(tmpdata[,stage_markers[i]],
                         by=list(tmp_pos),min)[,2])
}
colnames(tmp_thres) <- c("Negative","Positive")
rownames(tmp_thres) <- stage_markers

write.table(cbind(tmp_thres,rownames(tmp_thres)),
      "tmp.txt",sep="\t",quote=F,row.names=F)


# Step 1.2: check distribution of data
selCol <- c("ROI","CellArea","ClassifierLabel","PhenotypeOrdered")
selMarker <- c("PanCK")
selTissue <- c("Tissue")  #*****
selCohort <- c("nonMCL")   #****** 

tmpdata <- allPhenotypes_nonMCL[,selCol]  #******
tmpdata$Intensity <- allPhenotypes_nonMCL[,selMarker] #******
tmpdata$positive <- grepl(selMarker,tmpdata$PhenotypeOrdered)
tmp4plot <- tmpdata %>%   
               dplyr::filter(ClassifierLabel==selTissue) %>%
               left_join(clin.data[,c("PCGA2_BiospecimenID","PCGA2_PatID","tissueType","His3")],
                  by=c("ROI"="PCGA2_BiospecimenID"))

tmp4title <- paste(selCohort,selTissue,selMarker,sep="+")  
ggplot(tmp4plot, aes(x=ROI,y=Intensity))+
  geom_violin()+ 
  scale_y_continuous(trans='log10')+ggtitle(tmp4title)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

# -------
tmp4title <- paste(selCohort,selTissue,sep="+")  
ggplot(tmp4plot, aes(x=ROI,y=CellArea))+
  geom_violin()+ 
  scale_y_continuous(trans='log10')+ggtitle(tmp4title)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

#----
# check cutoff
selCol <- c("ROI","CellArea","ClassifierLabel","PhenotypeOrdered")
selMarker <- c("PD1")
selTissue <- c("Tissue")  #*****
selCohort <- c("nonMCL")   #****** 

# plot Tissue
tmpdata <- allPhenotypes_nonMCL[,selCol]  #******
tmpdata$Intensity <- allPhenotypes_nonMCL[,selMarker] #******
tmpdata$positive <- grepl(selMarker,tmpdata$PhenotypeOrdered)
tmp4plot <- tmpdata %>%   
               dplyr::filter(ClassifierLabel==selTissue) %>%
               left_join(clin.data[,c("PCGA2_BiospecimenID","PCGA2_PatID","tissueType","His3")],
                  by=c("ROI"="PCGA2_BiospecimenID"))

tmp4title <- paste(selCohort,selTissue,selMarker,sep="+")  
ggplot(tmp4plot, aes(x=positive,y=Intensity))+
  geom_violin()+ 
  scale_y_continuous(trans='log10')+ggtitle(tmp4title)+
  theme_bw()

aggregate(tmp4plot$Intensity,by=list(tmp4plot$positive),min)

# plot tissue and glass
tmpdata <- allPhenotypes_nonMCL[,selCol]  #******
tmpdata$Intensity <- allPhenotypes_nonMCL[,selMarker] #******
tmpdata$positive <- grepl(selMarker,tmpdata$PhenotypeOrdered)
tmp4plot <- tmpdata %>%   
               left_join(clin.data[,c("PCGA2_BiospecimenID","PCGA2_PatID","tissueType","His3")],
                  by=c("ROI"="PCGA2_BiospecimenID"))

tmp4title <- paste(selCohort,selMarker,sep="+")  
ggplot(tmp4plot, aes(x=ClassifierLabel,y=Intensity))+
  geom_violin()+ 
  scale_y_continuous(trans='log10')+ggtitle(tmp4title)+
  theme_bw()

# Step 2: Check combinations
selCol <- c("ROI","CellArea","ClassifierLabel","PhenotypeOrdered")
selTissue <- c("Tissue")  #*****
selCohort <- c("nonMCL")   #****** 

selMarker <- c("CD68")   #****
selComb_pos <- c("CD8")  #****
selComb_neg <- c("CD8,PanCK")  #****

#---
tmpdata <- allPhenotypes_nonMCL[,selCol]  #******
tmpdata$Intensity <- allPhenotypes_nonMCL[,selMarker] #******
tmp4plot <- tmpdata %>%   
               dplyr::filter(PhenotypeOrdered %in% selComb_pos) %>%
               left_join(clin.data[,c("PCGA2_BiospecimenID","PCGA2_PatID","tissueType","His3")],
                  by=c("ROI"="PCGA2_BiospecimenID"))

ggplot(tmp4plot, aes(x=ROI,y=Intensity))+
  geom_violin()+ 
  scale_y_continuous(trans='log10')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

ggplot(tmp4plot, aes(x=His3,y=Intensity))+
  geom_violin()+ 
  scale_y_continuous(trans='log10')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))


# ----
tmpdata <- allPhenotypes_mcl[,selCol]  #******
tmpdata$Intensity <- allPhenotypes_mcl[,selMarker] #******
tmp4plot <- tmpdata %>%  
               dplyr::filter(ClassifierLabel==selTissue) %>%
               dplyr::filter(PhenotypeOrdered %in% c(selComb_pos,selComb_neg)) %>%
               dplyr::mutate(PosCell = PhenotypeOrdered %in% selComb_pos) %>%
               left_join(clin.data[,c("PCGA2_BiospecimenID","PCGA2_PatID","tissueType","His3")],
                  by=c("ROI"="PCGA2_BiospecimenID"))

ggplot(tmp4plot, aes(x=PosCell,y=Intensity))+
  geom_violin()+ 
  scale_y_continuous(trans='log10')+
  theme_bw()




# *********************************************************************************************
# Re-define postive status of all markers > re-define cellTypes > calculate density
# > compared with panel 2 data
# Final thresholds:
stage_markers <- c("CD3","CD68","CD8","PanCK","PD1","PDL1","DAPI")
mcl_thres <- c(1.0,0.2,1.0,3.0,1.0,3.0,1.0)
nonMCL_thres <- c(1.0,0.2,0.5,5.0,1.0,1.0,1.0)
names(mcl_thres) <- stage_markers[1:7]
names(nonMCL_thres) <- stage_markers[1:7]
# Note: dictionary file: 01/07/25: 
#       tier1 has 1 row as CD8PD1--> CD8 on 2/5/25
#       tier1 has 1 row as NonMacMyeloud --> Myeloid on 2/5/25
#       strict_tier1: clean up
#       CD68,CD8,PD1: CD8PD1 -> Other
# 1/6/25
# *********************************************************************************************
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

# Step 0: load data
# load MCL data
fin_objData <- c("Panel1_MCL_colDataExpr_DAPIfilter_122324QC.rds")  #****
fin_sumData <- c("Panel1_MCL_summary_DAPIfilter_122324.rds") #*****

allPhenotypes_mcl <- readRDS(fin_objData)
allSummary_mcl <- readRDS(fin_sumData)

# load nonMCL data
fin_objData <- c("Panel1_nonMCL_colDataExpr_DAPIfilter_122324QC.rds")  #****
fin_sumData <- c("Panel1_nonMCL_summary_DAPIfilter_122324.rds")   #****

allPhenotypes_nonMCL <- readRDS(fin_objData)
allSummary_nonMCL <- readRDS(fin_sumData)

# import clinical data
fin_clin <- c("MIF_panel1hu_all_sampleInfo_122324.txt")
clin.data <- read.delim(fin_clin,sep="\t",header=T,stringsAsFactors=FALSE)  

clin.data$His3 <- as.character(clin.data$tissueType)
clin.data$His3 <- ifelse(is.element(clin.data$His3,c("MIA","ADC")),c("MIA_ADC"),clin.data$His3)
clin.data$His3 <- factor(clin.data$His3,
         levels=c("Normal","AAH","AIS","MIA_ADC"))
clin.data$tissueType <- factor(clin.data$tissueType,
         levels=c("Normal","AAH","AIS","MIA","ADC"))

colCode4Hist <- c("#999999","#009E73","#56B4E9","#0072B2","#D55E00")
names(colCode4Hist) <- c("Normal","AAH","AIS","MIA_ADC")
 
# load summary data
fin_sumData <- c("Panel1_all_summary_DAPIfilter_122324.rds")
allSummary_data <- readRDS(fin_sumData)

# load dictionary data
##fin_dictionary  <- paste("D:/MIF/Panel1_SU2C4Sub_031224/",
##                      "Markers_CT_conversion_P1_122324.txt",sep="")   #*******

fin_dictionary  <- paste("D:/MIF/Panel1_SU2C4Sub_031224/",
                      "Markers_CT_conversion_P1_010725.txt",sep="")   #*******
cellType_allDef <- read.delim(fin_dictionary,sep="\t",header=T,
                      stringsAsFactors=FALSE) 
 
# Step 1: re-define positive threshold per marker
# Step 1.0: setup threshold based on quantile of tissue cells
stage_markers <- c("CD3","CD68","CD8","PanCK","PD1","PDL1","DAPI")

#---
#mcl_thres <- c(2.0,0.2,1.5,1.0,0.5,1.5,3.0)
#nonMCL_thres <- c(1.0,1.2,0.5,6,1.0,1.0,1.0)

mcl_thres <- c(5.0,0.5,1.5,1.0,0.5,1.5,3.0)
nonMCL_thres <- c(1.0,1.2,0.5,1.0,1.0,1.0,1.0)

# non much change
##mcl_thres <- c(5.0,0.5,1.5,1.0,0.5,1.5,3.0)
##nonMCL_thres <- c(1.0,1.2,0.5,5.0,1.0,1.0,1.0)

mcl_thres <- c(5.0,0.75,1.5,1.0,1.0,1.5,3.0)
nonMCL_thres <- c(1.0,1.2,0.5,1.0,1.0,0.05,1.0)

mcl_thres <- c(5.0,0.75,1.5,1.0,1.0,1.5,3.0)
nonMCL_thres <- c(1.0,0.2,0.5,1.0,1.0,1.0,1.0)


mcl_thres <- c(1.0,0.2,2.0,1.0,1.0,2.0,1.0)
nonMCL_thres <- c(1.0,0.2,0.5,1.0,1.0,1.0,1.0)


mcl_thres <- c(1.0,0.2,1.0,3.0,1.0,3.0,10.0)
nonMCL_thres <- c(1.0,0.2,0.5,5.0,1.0,1.0,1.0)

names(mcl_thres) <- stage_markers[1:7]
names(nonMCL_thres) <- stage_markers[1:7]

# Step 1.1: extract postive status based on the new threshold
mcl_new_status <- intensity_to_binary_phenotype(allPhenotypes_mcl[,stage_markers],
                          mcl_thres)
nonMCL_new_status <- intensity_to_binary_phenotype(allPhenotypes_nonMCL[,stage_markers],
                           nonMCL_thres)

# ------
tmp <- table(mcl_new_status$Phenotype)
write.table(tmp, "tmp.txt",sep="\t",row.names=T,quote=F)

tmp <- table(nonMCL_new_status$Phenotype)
write.table(tmp, "tmp.txt",sep="\t",row.names=T,quote=F)

# check "Negative_all" in tissue segmentation
tmpneg <- which(mcl_new_status$Phenotype=="Negative_all")
table(allPhenotypes_mcl$ClassifierLabel[tmpneg])

tmpneg <- which(nonMCL_new_status$Phenotype=="Negative_all")
table(allPhenotypes_nonMCL$ClassifierLabel[tmpneg])
# -----

# step 1.3: convert to new cell type
cellType_tmp <- cellType_allDef[,c("phenotype_order","Tier1")]   #******
mcl_new_status$pheno_tier1  <- new_cellType(mcl_new_status$Phenotype,
                                     cellType_tmp)   #****
nonMCL_new_status$pheno_tier1  <- new_cellType(nonMCL_new_status$Phenotype,
                                     cellType_tmp)   #******

tmpSelVars <- c("Cell.ID","ROI","CellArea","ClassifierLabel","PhenotypeOrdered","pheno_tier1",stage_markers)
allPhenotypes_mcl_new <- allPhenotypes_mcl[,tmpSelVars] %>%
                            dplyr::rename(oldPhenotypeOrdered=PhenotypeOrdered,old_tier1=pheno_tier1)
allPhenotypes_mcl_new <- bind_cols(allPhenotypes_mcl_new,
                              mcl_new_status[,c("Phenotype","pheno_tier1")])

allPhenotypes_nonMCL_new <- allPhenotypes_nonMCL[,tmpSelVars] %>%
                              dplyr::rename(oldPhenotypeOrdered=PhenotypeOrdered,old_tier1=pheno_tier1)
allPhenotypes_nonMCL_new <- bind_cols(allPhenotypes_nonMCL_new,
                              nonMCL_new_status[,c("Phenotype","pheno_tier1")])

# combine mcl and nonMCL data
allPhenotypes_new <- bind_rows(allPhenotypes_mcl_new,allPhenotypes_nonMCL_new)

# Step 3: clean up Panel 1
# Step 3.1: clean up single-cell objects for Panel 1
# Filter out "Removed" cells (non-biological combinations)
allPhenotypes_filter <- allPhenotypes_new
table(allPhenotypes_filter$pheno_tier1)

tmpi <- which(allPhenotypes_filter$pheno_tier1=="Negative_all")
allPhenotypes_filter$pheno_tier1[tmpi] <- c("Removed")

keptR <- allPhenotypes_filter$pheno_tier1 != c("Removed")
allPhenotypes_filter <- allPhenotypes_filter[keptR,]

# Filter out "doublets" by cell area
quantile(allPhenotypes_filter$CellArea,c(0.9,0.95,0.99))

area_thres <- 100 #*****
keptR <- allPhenotypes_filter$CellArea<area_thres  #****
allPhenotypes_filter <- allPhenotypes_filter[keptR,]

# *****Filter out "Glass" cells for PANEL 1
selTissue <- c("Tissue")  #***** 
keptR <- allPhenotypes_filter$ClassifierLabel==selTissue  #****
allPhenotypes_filter <- allPhenotypes_filter[keptR,]

# # Step 3.1: clean up summary data for Panel 1
# Re-count cells in Tissue.Total.Cells
tmp_count <- table(allPhenotypes_filter$ROI)

tmp_count_new <- left_join(allSummary_data[,c("Analysis.Region","Total.Cells")],
             data.frame(ROI=names(tmp_count),Tissue.Total.Cells=as.numeric(tmp_count)),
             by=c("Analysis.Region"="ROI"))

cor(allSummary_data$Tissue.Total.Cells,tmp_count_new$Tissue.Total.Cells)

# Update summary data
allSummary_data_filter <- allSummary_data
allSummary_data_filter$Tissue.Total.Cells <- tmp_count_new$Tissue.Total.Cells

p1.allSummary_filter <- allSummary_data_filter
p1.allPhenotypes_filter <- allPhenotypes_filter
p1.clin_data <- clin.data

# Step 4:import Panel 2 data and clean up
# Step 4.1: import data
fin_objData <- c("Panel2_GGOpcga_colDataExpr_DAPIfilter_122324.rds")
fin_sumData <- c("Panel2_GGOpcga_summary_DAPIfilter_122324.rds")

allPhenotypes <- readRDS(fin_objData)
allSummary_data <- readRDS(fin_sumData)

fin_clin <- c("MIF_panel2hu_GGOpcga_sampleInfo_110924.txt")
clin.data <- read.delim(fin_clin,sep="\t",header=T,stringsAsFactors=FALSE) 

# update pheno_tier1
tmpi <- which(allPhenotypes$pheno_tier1=="PanCK,GranzymeB,FoxP3,Ki67,CD8")
allPhenotypes$pheno_tier1[tmpi] <- c("Removed")
allPhenotypes$pheno_tier2[tmpi] <- c("Removed")

# Step 4.2: clean up single-cell objects
# Filter out "Removed" cells (non-biological combinations)
allPhenotypes_filter <- allPhenotypes

table(allPhenotypes_filter$pheno_tier2)

keptR <- allPhenotypes_filter$pheno_tier2 != c("Removed")
allPhenotypes_filter <- allPhenotypes_filter[keptR,]

# Filter out "doublets" by cell area
quantile(allPhenotypes_filter$CellArea,c(0.9,0.95,0.99))

area_thres <- 100 #*****
keptR <- allPhenotypes_filter$CellArea<area_thres  #****
allPhenotypes_filter <- allPhenotypes_filter[keptR,]

# Step 4.3: clean up summary data
# Re-count cells in Tissue.Total.Cells
tmp_count <- table(allPhenotypes_filter$ROI)

tmp_count_new <- left_join(allSummary_data[,c("Analysis.Region","Total.Cells")],
             data.frame(ROI=names(tmp_count),Tissue.Total.Cells=as.numeric(tmp_count)),
             by=c("Analysis.Region"="ROI"))

cor(allSummary_data$Tissue.Total.Cells,tmp_count_new$Tissue.Total.Cells)

# Update summary data
allSummary_data_filter <- allSummary_data
allSummary_data_filter$Tissue.Total.Cells <- tmp_count_new$Tissue.Total.Cells

p2.allSummary_filter <- allSummary_data_filter
p2.allPhenotypes_filter <- allPhenotypes_filter
p2.clin_data <- clin.data
#------

# Step 5: calculate density
selCol4Sum <- c("Analysis.Region","Tissue.Area","Tissue.Total.Cells")

p1_cellNum <- table(p1.allPhenotypes_filter$ROI,
                    p1.allPhenotypes_filter$pheno_tier1)   #*****
p2_cellNum <- table(p2.allPhenotypes_filter$ROI,
                    p2.allPhenotypes_filter$pheno_tier2)   #*****

selCellType <- c("CD4")
p1_tmp_density <- data.frame(ROI=rownames(p1_cellNum),cellNo = p1_cellNum[,selCellType]) %>% 
             left_join(p1.allSummary_filter[,selCol4Sum],by=c("ROI"="Analysis.Region")) %>%
             dplyr::mutate(density = cellNo/Tissue.Area)
p2_tmp_density <- data.frame(ROI=rownames(p2_cellNum),cellNo = p2_cellNum[,selCellType]) %>% 
             left_join(p2.allSummary_filter[,selCol4Sum],by=c("ROI"="Analysis.Region")) %>%
             dplyr::mutate(density = cellNo/Tissue.Area) 

keptR <- is.element(p1.clin_data$PCGA2_BiospecimenID,p2.clin_data$PCGA2_BiospecimenID)
all_tmp_density <- p1.clin_data[keptR,-c(1:3)] %>%
             left_join(p1_tmp_density[,c("ROI","density")],
                   by=c("PCGA2_BiospecimenID"="ROI")) %>%
             dplyr::rename(p1_density=density) %>%
             left_join(p2_tmp_density[,c("ROI","density")],
                   by=c("PCGA2_BiospecimenID"="ROI")) %>%
             dplyr::rename(p2_density=density) 
all_tmp_density$cohortMain <- ifelse(all_tmp_density$cohort=="MCL","MCL","nonMCL")

cor(all_tmp_density$p1_density,all_tmp_density$p2_density)

tmpsel <- which(all_tmp_density$cohort=="MCL")
cor(all_tmp_density$p1_density[tmpsel],all_tmp_density$p2_density[tmpsel])

tmpsel2 <- which(all_tmp_density$cohort!="MCL")
cor(all_tmp_density$p1_density[tmpsel2],all_tmp_density$p2_density[tmpsel2])

tmp4title <- selCellType
ggplot(all_tmp_density, aes(x=p1_density,y=p2_density,color=cohortMain))+
  geom_point()+ 
  geom_smooth(method=lm)+
  geom_abline()+
  ggtitle(tmp4title)+
  theme_bw()


# ------
# check distribution of markers on Tissue vs. Glass based on the new cutoff values
selCol <- c("ROI","CellArea","ClassifierLabel","Phenotype")  #*****
selMarker <- c("CD3")
selTissue <- c("Tissue")  #*****
selCohort <- c("MCL")   #****** 

# plot Tissue
tmpdata <- allPhenotypes_mcl_new[,selCol]  #******
tmpdata$Intensity <- allPhenotypes_mcl_new[,selMarker] #******
tmpdata$positive <- grepl(selMarker,tmpdata$Phenotype)
tmp4plot <- tmpdata %>%   
               dplyr::filter(ClassifierLabel==selTissue)

tmp4title <- paste(selCohort,selTissue,selMarker,sep="+")  
ggplot(tmp4plot, aes(x=positive,y=Intensity))+
  geom_violin()+ 
  scale_y_continuous(trans='log10')+ggtitle(tmp4title)+
  theme_bw()

# plot tissue and glass
tmpdata <- allPhenotypes_mcl_new[,selCol]  #******
tmpdata$Intensity <- allPhenotypes_mcl_new[,selMarker] #******
tmpdata$positive <- grepl(selMarker,tmpdata$PhenotypeOrdered)
tmp4plot <- tmpdata 

tmp4title <- paste(selCohort,selMarker,sep="+")  
ggplot(tmp4plot, aes(x=ClassifierLabel,y=Intensity))+
  geom_violin()+ 
  scale_y_continuous(trans='log10')+ggtitle(tmp4title)+
  theme_bw()


tmp <- table(allPhenotypes_nonMCL_new$oldPhenotypeOrdered,allPhenotypes_nonMCL_new$Phenotype)
pheatmap(tmp,scale="row",cluster_rows=FALSE,cluster_cols=FALSE)

tmp <- table(allPhenotypes_mcl_new$oldPhenotypeOrdered,allPhenotypes_mcl_new$Phenotype)
pheatmap(tmp,scale="row",cluster_rows=FALSE,cluster_cols=FALSE)

# ***********************************************************************************
# Calculate the density for panel 1 - using adjusted cell identiy based on
# new cutoff
# 2/5/2025
# ***********************************************************************************
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

# Step 0: load data
# load MCL data
fin_objData <- c("Panel1_MCL_colDataExpr_DAPIfilter_122324QC.rds")  #****
fin_sumData <- c("Panel1_MCL_summary_DAPIfilter_122324.rds") #*****

allPhenotypes_mcl <- readRDS(fin_objData)
allSummary_mcl <- readRDS(fin_sumData)

# load nonMCL data
fin_objData <- c("Panel1_nonMCL_colDataExpr_DAPIfilter_122324QC.rds")  #****
fin_sumData <- c("Panel1_nonMCL_summary_DAPIfilter_122324.rds")   #****

allPhenotypes_nonMCL <- readRDS(fin_objData)
allSummary_nonMCL <- readRDS(fin_sumData)

# import clinical data
fin_clin <- c("MIF_panel1hu_all_sampleInfo_122324.txt")
clin.data <- read.delim(fin_clin,sep="\t",header=T,stringsAsFactors=FALSE)  

clin.data$His3 <- as.character(clin.data$tissueType)
clin.data$His3 <- ifelse(is.element(clin.data$His3,c("MIA","ADC")),c("MIA_ADC"),clin.data$His3)
clin.data$His3 <- factor(clin.data$His3,
         levels=c("Normal","AAH","AIS","MIA_ADC"))
clin.data$tissueType <- factor(clin.data$tissueType,
         levels=c("Normal","AAH","AIS","MIA","ADC"))

colCode4Hist <- c("#999999","#009E73","#56B4E9","#0072B2","#D55E00")
names(colCode4Hist) <- c("Normal","AAH","AIS","MIA_ADC")
 
# load summary data
fin_sumData <- c("Panel1_all_summary_DAPIfilter_122324.rds")
allSummary_data <- readRDS(fin_sumData)

# load dictionary data
##fin_dictionary  <- paste("D:/MIF/Panel1_SU2C4Sub_031224/",
##                      "Markers_CT_conversion_P1_122324.txt",sep="")   #*******

fin_dictionary  <- paste("D:/MIF/Panel1_SU2C4Sub_031224/",
                      "Markers_CT_conversion_P1_010725.txt",sep="")   #*******
cellType_allDef <- read.delim(fin_dictionary,sep="\t",header=T,
                      stringsAsFactors=FALSE) 
 
# Step 1: re-define positive threshold per marker
# Step 1.0: setup threshold based on quantile of tissue cells
stage_markers <- c("CD3","CD68","CD8","PanCK","PD1","PDL1","DAPI")

mcl_thres <- c(1.0,0.2,1.0,3.0,1.0,3.0,1.0)
nonMCL_thres <- c(1.0,0.2,0.5,5.0,1.0,1.0,1.0)

names(mcl_thres) <- stage_markers[1:7]
names(nonMCL_thres) <- stage_markers[1:7]

# Step 1.1: extract postive status based on the new threshold
mcl_new_status <- intensity_to_binary_phenotype(allPhenotypes_mcl[,stage_markers],
                          mcl_thres)
nonMCL_new_status <- intensity_to_binary_phenotype(allPhenotypes_nonMCL[,stage_markers],
                           nonMCL_thres)


# step 1.2: convert to new cell type
# ****Using tier 1
cellType_tmp <- cellType_allDef[,c("phenotype_order","Tier1")]   #******
mcl_new_status$pheno_tier1  <- new_cellType(mcl_new_status$Phenotype,
                                     cellType_tmp)   #****
nonMCL_new_status$pheno_tier1  <- new_cellType(nonMCL_new_status$Phenotype,
                                     cellType_tmp)   #******

# ****Using tier 2 - subtype
cellType_tmp <- cellType_allDef[,c("phenotype_order","lenient_phenotype")]   #******
mcl_new_status$pheno_lenientSub <- new_cellType(mcl_new_status$Phenotype,
                                     cellType_tmp)   #****
nonMCL_new_status$pheno_lenientSub  <- new_cellType(nonMCL_new_status$Phenotype,
                                     cellType_tmp)   #******


tmpSelVars <- c("Cell.ID","ROI","CellArea","ClassifierLabel","PhenotypeOrdered","pheno_tier1",stage_markers)
allPhenotypes_mcl_new <- allPhenotypes_mcl[,tmpSelVars] %>%
                            dplyr::rename(oldPhenotypeOrdered=PhenotypeOrdered,old_tier1=pheno_tier1)
allPhenotypes_mcl_new <- bind_cols(allPhenotypes_mcl_new,
                              mcl_new_status[,c("Phenotype","pheno_tier1","pheno_lenientSub")])

allPhenotypes_nonMCL_new <- allPhenotypes_nonMCL[,tmpSelVars] %>%
                              dplyr::rename(oldPhenotypeOrdered=PhenotypeOrdered,old_tier1=pheno_tier1)
allPhenotypes_nonMCL_new <- bind_cols(allPhenotypes_nonMCL_new,
                              nonMCL_new_status[,c("Phenotype","pheno_tier1","pheno_lenientSub")])

# combine mcl and nonMCL data
allPhenotypes_new <- bind_rows(allPhenotypes_mcl_new,allPhenotypes_nonMCL_new)

# Step 2: clean up Panel 1
# Step 2.1: clean up single-cell objects for Panel 1
# Filter out "Removed" cells (non-biological combinations)
allPhenotypes_filter <- allPhenotypes_new
table(allPhenotypes_filter$pheno_tier1)

tmpi <- which(allPhenotypes_filter$pheno_tier1=="Negative_all")
allPhenotypes_filter$pheno_tier1[tmpi] <- c("Removed")

keptR <- allPhenotypes_filter$pheno_tier1 != c("Removed")
allPhenotypes_filter <- allPhenotypes_filter[keptR,]

# Filter out "doublets" by cell area
quantile(allPhenotypes_filter$CellArea,c(0.9,0.95,0.99))

area_thres <- 100 #*****
keptR <- allPhenotypes_filter$CellArea<area_thres  #****
allPhenotypes_filter <- allPhenotypes_filter[keptR,]

# *****Filter out "Glass" cells for PANEL 1
selTissue <- c("Tissue")  #***** 
keptR <- allPhenotypes_filter$ClassifierLabel==selTissue  #****
allPhenotypes_filter <- allPhenotypes_filter[keptR,]

# Step 2.2: clean up summary data for Panel 1
# Re-count cells in Tissue.Total.Cells
tmp_count <- table(allPhenotypes_filter$ROI)

tmp_count_new <- left_join(allSummary_data[,c("Analysis.Region","Total.Cells")],
             data.frame(ROI=names(tmp_count),Tissue.Total.Cells=as.numeric(tmp_count)),
             by=c("Analysis.Region"="ROI"))

cor(allSummary_data$Tissue.Total.Cells,tmp_count_new$Tissue.Total.Cells)

# Update summary data
allSummary_data_filter <- allSummary_data
allSummary_data_filter$Tissue.Total.Cells <- tmp_count_new$Tissue.Total.Cells

p1.allSummary_filter <- allSummary_data_filter
p1.allPhenotypes_filter <- allPhenotypes_filter
p1.clin_data <- clin.data

# Step 3: calculate density
# calculate cel number per samples --> define subtype level
selTier = c("pheno_lenientSub") #***including sub-lineages
tmp_cellNum <- table(allPhenotypes_filter$ROI,
                     allPhenotypes_filter$pheno_lenientSub)    #********

tmp_cellNum_short <- data.frame(rbind(tmp_cellNum)) %>%
     rownames_to_column("ROI")


selCol4Sum <- c("Analysis.Region","Tissue.Area","Tissue.Total.Cells")

data4cellNum <- p1.clin_data[,-c(1:3)] %>%
      left_join(p1.allSummary_filter[,selCol4Sum],
           by=c("PCGA2_BiospecimenID"="Analysis.Region")) %>%
      left_join(tmp_cellNum_short,by=c("PCGA2_BiospecimenID"="ROI"))


# remove total cell to avoid confusion --> encourage using rowSums
data4cellNum$Tissue.Total.Cells<- NULL

# -----
fout <- c("Panel1_GGOpcga_cellNumber4BU_020525.rds") 
saveRDS(data4cellNum,file=fout)
# ----




p1_cellNum <- table(p1.allPhenotypes_filter$ROI,
                    p1.allPhenotypes_filter$pheno_tier1)   #*****






















# ***********************************************************************************
# Calculate the density
# 12/26/24
# ***********************************************************************************
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

# Step 1: import Panel 1 data and clean up
# Step 1.1: import data
fin_objData <- c("Panel1_all_colDataExpr_DAPIfilter_122324.rds")
fin_sumData <- c("Panel1_all_summary_DAPIfilter_122324.rds")

allPhenotypes <- readRDS(fin_objData)
allSummary_data <- readRDS(fin_sumData)

# update pheno_tier1
tmpidx <- which(allPhenotypes$pheno_tier1==c("CD8PD1"))
allPhenotypes$pheno_tier1[tmpidx] <- c("CD8")

fin_clin <- c("MIF_panel1hu_all_sampleInfo_122324.txt")
clin.data <- read.delim(fin_clin,sep="\t",header=T,stringsAsFactors=FALSE) 

# Step 1.2: clean up single-cell objects
# Filter out "Removed" cells (non-biological combinations)
allPhenotypes_filter <- allPhenotypes

table(allPhenotypes_filter$pheno_lenient)

keptR <- allPhenotypes_filter$pheno_lenient != c("Removed")
allPhenotypes_filter <- allPhenotypes_filter[keptR,]

# Filter out "doublets" by cell area
quantile(allPhenotypes_filter$CellArea,c(0.9,0.95,0.99))

area_thres <- 100 #*****
keptR <- allPhenotypes_filter$CellArea<area_thres  #****
allPhenotypes_filter <- allPhenotypes_filter[keptR,]

# *****Filter out "Glass" cells for PANEL 1
selTissue <- c("Tissue")  #***** 
keptR <- allPhenotypes_filter$ClassifierLabel==selTissue  #****
allPhenotypes_filter <- allPhenotypes_filter[keptR,]

# update pheno_tier1
##tmpidx <- which(allPhenotypes_filter$pheno_tier1==c("CD8PD1"))
##allPhenotypes_filter$pheno_tier1[tmpidx] <- c("CD8")


# Step 1.3: clean up summary data
# Re-count cells in Tissue.Total.Cells
tmp_count <- table(allPhenotypes_filter$ROI)

tmp_count_new <- left_join(allSummary_data[,c("Analysis.Region","Total.Cells")],
             data.frame(ROI=names(tmp_count),Tissue.Total.Cells=as.numeric(tmp_count)),
             by=c("Analysis.Region"="ROI"))

cor(allSummary_data$Tissue.Total.Cells,tmp_count_new$Tissue.Total.Cells)

# Update summary data
allSummary_data_filter <- allSummary_data
allSummary_data_filter$Tissue.Total.Cells <- tmp_count_new$Tissue.Total.Cells

p1.allSummary_filter <- allSummary_data_filter
p1.allPhenotypes_filter <- allPhenotypes_filter
p1.clin_data <- clin.data

# Step 2:import Panel 2 data and clean up
# Step 2.1: import data
fin_objData <- c("Panel2_GGOpcga_colDataExpr_DAPIfilter_122324.rds")
fin_sumData <- c("Panel2_GGOpcga_summary_DAPIfilter_122324.rds")

allPhenotypes <- readRDS(fin_objData)
allSummary_data <- readRDS(fin_sumData)

fin_clin <- c("MIF_panel2hu_GGOpcga_sampleInfo_110924.txt")
clin.data <- read.delim(fin_clin,sep="\t",header=T,stringsAsFactors=FALSE) 

# update pheno_tier1
tmpi <- which(allPhenotypes$pheno_tier1=="PanCK,GranzymeB,FoxP3,Ki67,CD8")
allPhenotypes$pheno_tier1[tmpi] <- c("Removed")
allPhenotypes$pheno_tier2[tmpi] <- c("Removed")

# Step 2.2: clean up single-cell objects
# Filter out "Removed" cells (non-biological combinations)
allPhenotypes_filter <- allPhenotypes

table(allPhenotypes_filter$pheno_tier2)

keptR <- allPhenotypes_filter$pheno_tier2 != c("Removed")
allPhenotypes_filter <- allPhenotypes_filter[keptR,]

# Filter out "doublets" by cell area
quantile(allPhenotypes_filter$CellArea,c(0.9,0.95,0.99))

area_thres <- 100 #*****
keptR <- allPhenotypes_filter$CellArea<area_thres  #****
allPhenotypes_filter <- allPhenotypes_filter[keptR,]

# Step 2.3: clean up summary data
# Re-count cells in Tissue.Total.Cells
tmp_count <- table(allPhenotypes_filter$ROI)

tmp_count_new <- left_join(allSummary_data[,c("Analysis.Region","Total.Cells")],
             data.frame(ROI=names(tmp_count),Tissue.Total.Cells=as.numeric(tmp_count)),
             by=c("Analysis.Region"="ROI"))

cor(allSummary_data$Tissue.Total.Cells,tmp_count_new$Tissue.Total.Cells)

# Update summary data
allSummary_data_filter <- allSummary_data
allSummary_data_filter$Tissue.Total.Cells <- tmp_count_new$Tissue.Total.Cells

p2.allSummary_filter <- allSummary_data_filter
p2.allPhenotypes_filter <- allPhenotypes_filter
p2.clin_data <- clin.data
#------

# Step 3: calculate density
selCol4Sum <- c("Analysis.Region","Tissue.Area","Tissue.Total.Cells")

p1_cellNum <- table(p1.allPhenotypes_filter$ROI,
                    p1.allPhenotypes_filter$pheno_tier1)   #*****
p2_cellNum <- table(p2.allPhenotypes_filter$ROI,
                    p2.allPhenotypes_filter$pheno_tier2)   #*****

selCellType <- c("CD4")
p1_tmp_density <- data.frame(ROI=rownames(p1_cellNum),cellNo = p1_cellNum[,selCellType]) %>% 
             left_join(p1.allSummary_filter[,selCol4Sum],by=c("ROI"="Analysis.Region")) %>%
             dplyr::mutate(density = cellNo/Tissue.Area)
p2_tmp_density <- data.frame(ROI=rownames(p2_cellNum),cellNo = p2_cellNum[,selCellType]) %>% 
             left_join(p2.allSummary_filter[,selCol4Sum],by=c("ROI"="Analysis.Region")) %>%
             dplyr::mutate(density = cellNo/Tissue.Area) 

keptR <- is.element(p1.clin_data$PCGA2_BiospecimenID,p2.clin_data$PCGA2_BiospecimenID)
all_tmp_density <- p1.clin_data[keptR,-c(1:3)] %>%
             left_join(p1_tmp_density[,c("ROI","density")],
                   by=c("PCGA2_BiospecimenID"="ROI")) %>%
             dplyr::rename(p1_density=density) %>%
             left_join(p2_tmp_density[,c("ROI","density")],
                   by=c("PCGA2_BiospecimenID"="ROI")) %>%
             dplyr::rename(p2_density=density) 
all_tmp_density$cohortMain <- ifelse(all_tmp_density$cohort=="MCL","MCL","nonMCL")

cor(all_tmp_density$p1_density,all_tmp_density$p2_density)

tmpsel <- which(all_tmp_density$cohort=="MCL")
cor(all_tmp_density$p1_density[tmpsel],all_tmp_density$p2_density[tmpsel])
cor(all_tmp_density$p1_density[tmpsel],all_tmp_density$p2_density[tmpsel],method="spearman")

tmpsel2 <- which(all_tmp_density$cohort!="MCL")
cor(all_tmp_density$p1_density[tmpsel2],all_tmp_density$p2_density[tmpsel2])
cor(all_tmp_density$p1_density[tmpsel2],all_tmp_density$p2_density[tmpsel2],method="spearman")

tmp4title <- selCellType
ggplot(all_tmp_density, aes(x=p1_density,y=p2_density,color=cohortMain))+
  geom_point()+ 
  geom_smooth(method=lm)+
  geom_abline()+
  ggtitle(tmp4title)+
  theme_bw()



# reference codes
# *******************************************************************************
# redefine positive status of all markers, including DAPI, and create
# Phenotype, which is combination of positive status for each cell
# modified from format_halo_to_sce_4DL_DAPIfilter
# 1/6/24
intensity_to_binary_phenotype <- function(intensity_of_markers,markers_thres){

    markers <- colnames(intensity_of_markers)
    expression_status_cols <- matrix(0,nrow(intensity_of_markers),length(markers))
    for (k in c(1:length(stage_markers))){
        tmpi <- which(intensity_of_markers[,k] > markers_thres[k])
        expression_status_cols[tmpi,k] <- 1
    }
    colnames(expression_status_cols) <- markers
    expression_status_cols <- as.data.frame(expression_status_cols)

    # start reading in the Phenotypes of every cell based on binary status 
    # ****Phenotype as concatenate of positive markers******
    expression_status_cols$Phenotype <- ""
    for (marker in markers) {
        if (marker == "DAPI") {
            phenotype <- "OTHER,"
        } else {
            phenotype <- paste(marker, ",", sep = "")
        }

        # get the row idx of the cells that express the specific marker, and paste the phenotype
        rows_true_exp <- which(expression_status_cols[,marker] != 0)
        if (length(rows_true_exp) != 0) {
            expression_status_cols[rows_true_exp,]$Phenotype <- paste(expression_status_cols[rows_true_exp,]$Phenotype, 
                    phenotype, sep="")
        }
    }

    # now clean the phenotype column - remove DAPI ("OTHER") marker in the concatenated list
    if (nrow(expression_status_cols[expression_status_cols$Phenotype == "", ]) != 0) {
        expression_status_cols[expression_status_cols$Phenotype == "", ]$Phenotype <- "Negative_all"
    }
    if (nrow(expression_status_cols[expression_status_cols$Phenotype == "OTHER,", ]) != 0) {
        expression_status_cols[expression_status_cols$Phenotype == "OTHER,", ]$Phenotype <- "OTHER"
    }
    expression_status_cols$Phenotype <- gsub("OTHER,", "", expression_status_cols$Phenotype)
    expression_status_cols$Phenotype <- gsub(",OTHER", "", expression_status_cols$Phenotype)
    expression_status_cols$Phenotype <- gsub(",$", "", expression_status_cols$Phenotype)
    return(expression_status_cols)
}

# ************************************************************************
# redefine cell types based on dictionary file
# PositivePhenotype: Phenotype vector listing positive marker in each cell
# cellType_def: matrix with 1) col 1: list combination of postive markers, 
# and 2) col 2: specify equivalent cell type
# combination in PositivePhenotype must be matched with those in col 1 of cellType_def
new_cellType <- function(PositivePhenotype,cellType_def){
   nct <- nrow(cellType_def)
   new_cellType <-  PositivePhenotype
   for (jj in c(1:nct)){
      new_cellType[PositivePhenotype==cellType_def[jj,1]] <- cellType_def[jj,2]   #****
   }
   return(new_cellType)
}

