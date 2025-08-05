sel_panel <- "Myeloid"

rds_subclust_dir <- paste0("D:/CHRISTOPHER_BUI/MPE_CYTOF_RDS/", sel_panel, "/Results_Subcluster/RDS_Subcluster")

rds_files <- list.files(rds_subclust_dir)

# sce_c1 <- readRDS(file.path(rds_subclust_dir, rds_files[1]))
sce_c234 <- readRDS(file.path(rds_subclust_dir, rds_files[2]))


sce_main <- readRDS(r'(D:\CHRISTOPHER_BUI\MPE_CYTOF_RDS\Myeloid\sce_main.rds)')
sce_fl <- readRDS(r'(D:\CHRISTOPHER_BUI\MPE_CYTOF_RDS\Myeloid\sce_full_lineage.rds)')

# ------------------------------------------------------------------------------
# NOTE:
#   - Not much batch effect observed, so we just analyze tissue type
#   - Use ALL markers
# ------------------------------------------------------------------------------


################################## C1 ##########################################
metak <- "meta9"

# filter for selected cluster
sce_main_tmp <- filterSCE(sce_main, meta6_full_lin %in% c(1))

# select SCE
sce_sub <- sce_c1

# transfer metadata
sce_main_tmp$cluster_id <- sce_sub$cluster_id
sce_main_tmp@metadata <- metadata(sce_sub)
sce_main_tmp$meta9 <- cluster_ids(sce_sub, "meta9")  # adjust as needed

# set PBMC as reference
colData(sce_main_tmp)$tissue_type <- factor(colData(sce_main_tmp)$tissue_type, levels = c("PBMC", "MPE"))

sce_tmp_all <- sce_main_tmp

# filter out refs
sce_tmp <- filterSCE(sce_tmp_all, patient_id != "Ref")


# OPTION 1: TISSUE ONLY --------------------------------------------------------
# y ~ 1 + tissue_type
metak <- "meta9"
tmp_metadata <- ei(sce_tmp)

design <- createDesignMatrix(tmp_metadata, cols_design = c("tissue_type"))
# colnames(design)
contrast <- createContrast(c(0,1))
# 0=intercept, 1=tissue_type

# all markers
sel_markers <- rownames(sce_tmp) %in% rownames(sce_tmp)

C1_res_DS_1 <- diffcyt(sce_tmp,
                       clustering_to_use = metak,
                       analysis_type = "DS",
                       method_DS = "diffcyt-DS-limma",
                       design = design,
                       contrast = contrast,
                       markers_to_test = sel_markers,
                       verbose = FALSE,
                       transform = FALSE)



################################## C2,C3,C4 ####################################
metak <- "meta6"

# filter for selected cluster
sce_main_tmp <- filterSCE(sce_main, meta6_full_lin %in% c(2,3,4))

# select SCE
sce_sub <- sce_c234

# transfer metadata
sce_main_tmp$cluster_id <- sce_sub$cluster_id
sce_main_tmp@metadata <- metadata(sce_sub)
sce_main_tmp$meta9 <- cluster_ids(sce_sub, "meta6")  # adjust as needed

# set PBMC as reference
colData(sce_main_tmp)$tissue_type <- factor(colData(sce_main_tmp)$tissue_type, levels = c("PBMC", "MPE"))

sce_tmp_all <- sce_main_tmp

# filter out refs
sce_tmp <- filterSCE(sce_tmp_all, patient_id != "Ref")


# OPTION 1: TISSUE ONLY --------------------------------------------------------
# y ~ 1 + tissue_type
metak <- "meta6"
tmp_metadata <- ei(sce_tmp)

design <- createDesignMatrix(tmp_metadata, cols_design = c("tissue_type"))
# colnames(design)
contrast <- createContrast(c(0,1))
# 0=intercept, 1=tissue_type

# all markers
sel_markers <- rownames(sce_tmp) %in% rownames(sce_tmp)

C234_res_DS_1 <- diffcyt(sce_tmp,
                         clustering_to_use = metak,
                         analysis_type = "DS",
                         method_DS = "diffcyt-DS-limma",
                         design = design,
                         contrast = contrast,
                         markers_to_test = sel_markers,
                         verbose = FALSE,
                         transform = FALSE)



# COMBINE RESULTS  -------------------------------------------------------------

C1_DS_1 <- as.data.frame(rowData(C1_res_DS_1$res))
C234_DS_1 <- as.data.frame(rowData(C234_res_DS_1$res))

# option 1
C1_DS_1 <- C1_DS_1 %>%
  select(cluster_id, marker_id, logFC, p_val, p_adj) %>%
  dplyr::rename(
    cluster = cluster_id,
    marker = marker_id,
    logFC_DS = logFC,
    p_val_DS = p_val,
    p_adj_DS = p_adj
  ) %>% mutate(
    model = "1: Tissue",
    tier1_cluster = "C1"
  )
# option 1
C234_DS_1 <- C234_DS_1 %>%
  select(cluster_id, marker_id, logFC, p_val, p_adj) %>%
  dplyr::rename(
    cluster = cluster_id,
    marker = marker_id,
    logFC_DS = logFC,
    p_val_DS = p_val,
    p_adj_DS = p_adj
  ) %>% mutate(
    model = "1: Tissue",
    tier1_cluster = "C234"
  )


# COMBINE ALL
DS_FINAL <- bind_rows(C1_DS_1, C234_DS_1)
rownames(DS_FINAL) <- NULL

file_name <- paste("Myeloid_subclust_DS.txt")
write.table(DS_FINAL, file.path(res_dir, file_name), sep = "\t", row.names = FALSE, quote = FALSE)




############################### TIER 1 #########################################

metak <- "meta6"

# select SCE
sce_fl$meta6 <- cluster_ids(sce_fl, "meta6")

# transfer metadata
sce_main$cluster_id <- sce_fl$cluster_id
sce_main@metadata <- metadata(sce_fl)
sce_main$meta6 <- sce_fl$meta6  # adjust as needed

# set PBMC as reference
colData(sce_main)$tissue_type <- factor(colData(sce_main)$tissue_type, levels = c("PBMC", "MPE"))

# filter out refs
sce_tmp <- filterSCE(sce_main, patient_id != "Ref")


# OPTION 1: TISSUE ONLY --------------------------------------------------------
# y ~ 1 + tissue_type
metak <- "meta6"
tmp_metadata <- ei(sce_tmp)

design <- createDesignMatrix(tmp_metadata, cols_design = c("tissue_type"))
# colnames(design)
contrast <- createContrast(c(0,1))
# 0=intercept, 1=tissue_type

# all markers
sel_markers <- rownames(sce_tmp) %in% rownames(sce_tmp)

TIER1_META6_res_DS_1 <- diffcyt(sce_tmp,
                                clustering_to_use = metak,
                                analysis_type = "DS",
                                method_DS = "diffcyt-DS-limma",
                                design = design,
                                contrast = contrast,
                                markers_to_test = sel_markers,
                                verbose = FALSE,
                                transform = FALSE)

TIER1_FINAL_DS <- as.data.frame(rowData(TIER1_META6_res_DS_1$res)) %>%
  select(cluster_id, marker_id, logFC, p_val, p_adj) %>%
  dplyr::rename(
    cluster = cluster_id,
    marker = marker_id,
    logFC_DS = logFC,
    p_val_DS = p_val,
    p_adj_DS = p_adj
  ) %>% mutate(
    model = "1: Tissue",
    tier1_k = "meta6"
  )


file_name <- paste("Myeloid_tier1_meta6_DS.txt")
write.table(TIER1_FINAL_DS, file.path(res_dir, file_name), sep = "\t", row.names = FALSE, quote = FALSE)














