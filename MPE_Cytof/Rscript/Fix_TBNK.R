wrongID <- c("TBNK_PID19-PBMC_01_batch5")

table(sub_sce_tier1$parental_sample_id)
table(sub_sce_tier1$tissue_type)


# update phenotype matrix for the sce object, named sub_sce_tier1
tmpi <- which(sub_sce_tier1@metadata$experiment_info$sample_id==wrongID)
tmp_tissue_new <- sub_sce_tier1@metadata$experiment_info$tissue_type   #**use temporary var to avoid error of modified variable adding new level
tmp_tissue_new[tmpi] <- c("MPE")
sub_sce_tier1@metadata$experiment_info$tissue_type <- factor(tmp_tissue_new, levels = c("PBMC", "MPE"))

# check if updates are correct
ei(sub_sce_tier1) # not updated



# update colDtaa for the sce object, named sub_sce_tier1
tmpi <- which(colData(sub_sce_tier1)$sample_id==wrongID)
tmp_tissue_new <- sub_sce_tier1$tissue_type
tmp_tissue_new[tmpi] <- c("MPE")
sub_sce_tier1$tissue_type <- factor(tmp_tissue_new, levels = c("PBMC", "MPE"))


tmp_psID_new <- as.character(sub_sce_tier1$parental_sample_id)
tmp_psID_new[tmpi] <- c("PID19-MPE")
sub_sce_tier1$parental_sample_id <- factor(tmp_psID_new)


# check if updates are correct
table(sub_sce_tier1$parental_sample_id)
table(sub_sce_tier1$tissue_type)



# save the updated sce object with associated file name/
