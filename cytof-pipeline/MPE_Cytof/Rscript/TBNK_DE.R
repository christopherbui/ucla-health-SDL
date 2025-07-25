sel_panel <- "TBNK"

rds_subclust_dir <- paste0("D:/CHRISTOPHER_BUI/MPE_CYTOF_RDS/", sel_panel, "/Results_Subcluster/RDS_Subcluster")

rds_files <- list.files(rds_subclust_dir)

sce_c1 <- readRDS(file.path(rds_subclust_dir, rds_files[1]))
# sce_c2 <- readRDS(file.path(rds_subclust_dir, rds_files[2]))
# sce_c3c4 <- readRDS(file.path(rds_subclust_dir, rds_files[3]))
# sce_c5c9 <- readRDS(file.path(rds_subclust_dir, rds_files[4]))
sce_c6c7 <- readRDS(file.path(rds_subclust_dir, rds_files[5]))
# sce_c8 <- readRDS(file.path(rds_subclust_dir, rds_files[6]))

# subset lineage
sce_main <- readRDS("D:/CHRISTOPHER_BUI/MPE_CYTOF_RDS/TBNK/sce_subset_lineage.rds")

# sce with all markers
# 1. remove outliers from meta9 C1
sce_main_c1 <- sce_main[, sce_main$meta9 == 1] # this has outliers; process below to remove outliers

# identify outliers via UMAP
umap_matrix <- reducedDim(sce_main_c1, "UMAP")
tmp_sel <- umap_matrix[, 1] < (-10) | umap_matrix[, 2] > 10
tmp_sel[is.na(tmp_sel)] <- FALSE

# initialize column marking outliers
colData(sce_main)$c1_outlier <- FALSE
# assign outlier selection to only C1
c1_idx <- which(sce_main$meta9 == 1)
colData(sce_main)$c1_outlier[c1_idx] <- tmp_sel

# view outliers
colData(sce_main)[sce_main$c1_outlier == TRUE, ]
# remove outliers
sce_main <- sce_main[, sce_main$c1_outlier != TRUE]
sce_main_c1_clean <- sce_main[, colData(sce_main)$meta9 == 1]




rds_subclust_dir <- c("D:/CHRISTOPHER_BUI/MPE_CYTOF_RDS/TBNK/Results_Subcluster/RDS_Subcluster")

# load sce
sce_main <- readRDS("D:/CHRISTOPHER_BUI/MPE_CYTOF_RDS/TBNK/sce_main_0.rds") # all markers
sce_sublin <- readRDS("D:/CHRISTOPHER_BUI/MPE_CYTOF_RDS/TBNK/sce_subset_lineage_no_outliers.rds") # subset lineage markers

sce_main <- sce_main[, sce_main$c1_outlier != TRUE] # remove outliers

sce_main$meta9_tier1 <- sce_sublin$meta9  # add tier1 meta9 to sce_main



################################## C1 ##########################################
metak <- "meta5"

# filter for selected cluster
sce_main_tmp <- filterSCE(sce_main, meta9_tier1 %in% c(1))

# select SCE
sce_sub <- sce_c1

# transfer metadata
sce_main_tmp$cluster_id <- sce_sub$cluster_id
sce_main_tmp@metadata <- metadata(sce_sub)
sce_main_tmp$meta5 <- cluster_ids(sce_sub, "meta5")  # adjust as needed

# set PBMC as reference
colData(sce_main_tmp)$tissue_type <- factor(colData(sce_main_tmp)$tissue_type, levels = c("PBMC", "MPE"))

sce_tmp_all <- sce_main_tmp

# filter out refs
sce_tmp <- filterSCE(sce_tmp_all, patient_id != "Ref")


# OPTION 1: TISSUE ONLY --------------------------------------------------------
# y ~ 1 + tissue_type
metak <- "meta5"
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

# OPTION 2: BATCH AS FIXED -----------------------------------------------------
# y ~ 1 + tissue_type + BATCH
metak <- "meta5"
tmp_metadata <- ei(sce_tmp)

design <- createDesignMatrix(tmp_metadata, cols_design = c("tissue_type", "BATCH"))
# colnames(design)
contrast <- createContrast(c(0,1,rep(0,6)))
# 0=intercept, 1=tissue_type, 0...= batch 2 to batch 7

# all markers
sel_markers <- rownames(sce_tmp) %in% rownames(sce_tmp)

C1_res_DS_2 <- diffcyt(sce_tmp,
                       clustering_to_use = metak,
                       analysis_type = "DS",
                       method_DS = "diffcyt-DS-limma",
                       design = design,
                       contrast = contrast,
                       markers_to_test = sel_markers,
                       verbose = FALSE,
                       transform = FALSE)

# OPTION 3: MIXED MODEL  -------------------------------------------------------
# y ~ tissue_type + BATCH + (1 | patient_id)
metak <- "meta5"
tmp_metadata <- ei(sce_tmp)

designF <- createFormula(tmp_metadata, cols_fixed = c("tissue_type", "BATCH"), cols_random = c("patient_id"))
# designF$formula
contrastF <- createContrast(c(1,0,0))
# 1=tissue_type, 0=batch, 0=random effect

# all markers
sel_markers <- rownames(sce_tmp) %in% rownames(sce_tmp)

C1_res_DS_3 <- diffcyt(sce_tmp,
                       clustering_to_use = metak,
                       analysis_type = "DS",
                       method_DS = "diffcyt-DS-LMM", # diffcyt-DS-limma doesn't work; use diffcyt-DS_LMM
                       design = designF,
                       contrast = contrastF,
                       markers_to_test = sel_markers,
                       verbose = FALSE,
                       transform = FALSE)


################################## C2 ##########################################
metak <- "meta6"

# filter for selected cluster
sce_main_tmp <- filterSCE(sce_main, meta9_tier1 %in% c(2))

# select SCE
sce_sub <- sce_c2

# transfer metadata
sce_main_tmp$cluster_id <- sce_sub$cluster_id
sce_main_tmp@metadata <- metadata(sce_sub)
sce_main_tmp$meta6 <- cluster_ids(sce_sub, "meta6")  # adjust as needed

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

C2_res_DS_1 <- diffcyt(sce_tmp,
                       clustering_to_use = metak,
                       analysis_type = "DS",
                       method_DS = "diffcyt-DS-limma",
                       design = design,
                       contrast = contrast,
                       markers_to_test = sel_markers,
                       verbose = FALSE,
                       transform = FALSE)

# OPTION 2: BATCH AS FIXED -----------------------------------------------------
# y ~ 1 + tissue_type + BATCH
metak <- "meta6"
tmp_metadata <- ei(sce_tmp)

design <- createDesignMatrix(tmp_metadata, cols_design = c("tissue_type", "BATCH"))
# colnames(design)
contrast <- createContrast(c(0,1,rep(0,6)))
# 0=intercept, 1=tissue_type, 0...= batch 2 to batch 7

# all markers
sel_markers <- rownames(sce_tmp) %in% rownames(sce_tmp)

C2_res_DS_2 <- diffcyt(sce_tmp,
                       clustering_to_use = metak,
                       analysis_type = "DS",
                       method_DS = "diffcyt-DS-limma",
                       design = design,
                       contrast = contrast,
                       markers_to_test = sel_markers,
                       verbose = FALSE,
                       transform = FALSE)

# OPTION 3: MIXED MODEL  -------------------------------------------------------
# y ~ tissue_type + BATCH + (1 | patient_id)
metak <- "meta6"
tmp_metadata <- ei(sce_tmp)

designF <- createFormula(tmp_metadata, cols_fixed = c("tissue_type", "BATCH"), cols_random = c("patient_id"))
# designF$formula
contrastF <- createContrast(c(1,0,0))
# 1=tissue_type, 0=batch, 0=random effect

# all markers
sel_markers <- rownames(sce_tmp) %in% rownames(sce_tmp)

C2_res_DS_3 <- diffcyt(sce_tmp,
                       clustering_to_use = metak,
                       analysis_type = "DS",
                       method_DS = "diffcyt-DS-LMM",  # diffcyt-DS-limma doesn't work; use diffcyt-DS_LMM
                       design = designF,
                       contrast = contrastF,
                       markers_to_test = sel_markers,
                       verbose = FALSE,
                       transform = FALSE)



################################## C3,C4 #######################################
metak <- "meta6"

sce_main_tmp <- filterSCE(sce_main, meta9_tier1 %in% c(3,4))  # filter for selected cluster

# select SCE
sce_sub <- sce_c3c4

# transfer metadata
sce_main_tmp$cluster_id <- sce_sub$cluster_id
sce_main_tmp@metadata <- metadata(sce_sub)
sce_main_tmp$meta6 <- cluster_ids(sce_sub, "meta6")  # adjust as needed

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

C3C4_res_DS_1 <- diffcyt(sce_tmp,
                       clustering_to_use = metak,
                       analysis_type = "DS",
                       method_DS = "diffcyt-DS-limma",
                       design = design,
                       contrast = contrast,
                       markers_to_test = sel_markers,
                       verbose = FALSE,
                       transform = FALSE)

# OPTION 2: BATCH AS FIXED -----------------------------------------------------
# y ~ 1 + tissue_type + BATCH
metak <- "meta6"
tmp_metadata <- ei(sce_tmp)

design <- createDesignMatrix(tmp_metadata, cols_design = c("tissue_type", "BATCH"))
# colnames(design)
contrast <- createContrast(c(0,1,rep(0,6)))
# 0=intercept, 1=tissue_type, 0...= batch 2 to batch 7

# all markers
sel_markers <- rownames(sce_tmp) %in% rownames(sce_tmp)

C3C4_res_DS_2 <- diffcyt(sce_tmp,
                       clustering_to_use = metak,
                       analysis_type = "DS",
                       method_DS = "diffcyt-DS-limma",
                       design = design,
                       contrast = contrast,
                       markers_to_test = sel_markers,
                       verbose = FALSE,
                       transform = FALSE)

# OPTION 3: MIXED MODEL  -------------------------------------------------------
# y ~ tissue_type + BATCH + (1 | patient_id)
metak <- "meta6"
tmp_metadata <- ei(sce_tmp)

designF <- createFormula(tmp_metadata, cols_fixed = c("tissue_type", "BATCH"), cols_random = c("patient_id"))
# designF$formula
contrastF <- createContrast(c(1,0,0))
# 1=tissue_type, 0=batch, 0=random effect

# all markers
sel_markers <- rownames(sce_tmp) %in% rownames(sce_tmp)

C3C4_res_DS_3 <- diffcyt(sce_tmp,
                       clustering_to_use = metak,
                       analysis_type = "DS",
                       method_DS = "diffcyt-DS-LMM",  # diffcyt-DS-limma doesn't work; use diffcyt-DS_LMM
                       design = designF,
                       contrast = contrastF,
                       markers_to_test = sel_markers,
                       verbose = FALSE,
                       transform = FALSE)



################################## C5,C9 #######################################
metak <- "meta6"

# filter for selected cluster
sce_main_tmp <- filterSCE(sce_main, meta9_tier1 %in% c(5,9))

# select SCE
sce_sub <- sce_c5c9

# transfer metadata
sce_main_tmp$cluster_id <- sce_sub$cluster_id
sce_main_tmp@metadata <- metadata(sce_sub)
sce_main_tmp$meta6 <- cluster_ids(sce_sub, "meta6")  # adjust as needed

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

C5C9_res_DS_1 <- diffcyt(sce_tmp,
                         clustering_to_use = metak,
                         analysis_type = "DS",
                         method_DS = "diffcyt-DS-limma",
                         design = design,
                         contrast = contrast,
                         markers_to_test = sel_markers,
                         verbose = FALSE,
                         transform = FALSE)

# OPTION 2: BATCH AS FIXED -----------------------------------------------------
# y ~ 1 + tissue_type + BATCH
metak <- "meta6"
tmp_metadata <- ei(sce_tmp)

design <- createDesignMatrix(tmp_metadata, cols_design = c("tissue_type", "BATCH"))
# colnames(design)
contrast <- createContrast(c(0,1,rep(0,6)))
# 0=intercept, 1=tissue_type, 0...= batch 2 to batch 7

# all markers
sel_markers <- rownames(sce_tmp) %in% rownames(sce_tmp)

C5C9_res_DS_2 <- diffcyt(sce_tmp,
                         clustering_to_use = metak,
                         analysis_type = "DS",
                         method_DS = "diffcyt-DS-limma",
                         design = design,
                         contrast = contrast,
                         markers_to_test = sel_markers,
                         verbose = FALSE,
                         transform = FALSE)

# OPTION 3: MIXED MODEL  -------------------------------------------------------
# y ~ tissue_type + BATCH + (1 | patient_id)
metak <- "meta6"
tmp_metadata <- ei(sce_tmp)

designF <- createFormula(tmp_metadata, cols_fixed = c("tissue_type", "BATCH"), cols_random = c("patient_id"))
# designF$formula
contrastF <- createContrast(c(1,0,0))
# 1=tissue_type, 0=batch, 0=random effect

# all markers
sel_markers <- rownames(sce_tmp) %in% rownames(sce_tmp)

C5C9_res_DS_3 <- diffcyt(sce_tmp,
                         clustering_to_use = metak,
                         analysis_type = "DS",
                         method_DS = "diffcyt-DS-LMM",  # diffcyt-DS-limma doesn't work; use diffcyt-DS_LMM
                         design = designF,
                         contrast = contrastF,
                         markers_to_test = sel_markers,
                         verbose = FALSE,
                         transform = FALSE)



################################## C6,C7 #######################################
metak <- "meta5"

# filter for selected cluster
sce_main_tmp <- filterSCE(sce_main, meta9_tier1 %in% c(6,7))

# select SCE
sce_sub <- sce_c6c7

# transfer metadata
sce_main_tmp$cluster_id <- sce_sub$cluster_id
sce_main_tmp@metadata <- metadata(sce_sub)
sce_main_tmp$meta5 <- cluster_ids(sce_sub, "meta5") # adjust as needed

# set PBMC as reference
colData(sce_main_tmp)$tissue_type <- factor(colData(sce_main_tmp)$tissue_type, levels = c("PBMC", "MPE"))

sce_tmp_all <- sce_main_tmp

# filter out refs
sce_tmp <- filterSCE(sce_tmp_all, patient_id != "Ref")


# OPTION 1: TISSUE ONLY --------------------------------------------------------
# y ~ 1 + tissue_type
metak <- "meta5"
tmp_metadata <- ei(sce_tmp)

design <- createDesignMatrix(tmp_metadata, cols_design = c("tissue_type"))
# colnames(design)
contrast <- createContrast(c(0,1))
# 0=intercept, 1=tissue_type

# all markers
sel_markers <- rownames(sce_tmp) %in% rownames(sce_tmp)

C6C7_res_DS_1 <- diffcyt(sce_tmp,
                         clustering_to_use = metak,
                         analysis_type = "DS",
                         method_DS = "diffcyt-DS-limma",
                         design = design,
                         contrast = contrast,
                         markers_to_test = sel_markers,
                         verbose = FALSE,
                         transform = FALSE)

# OPTION 2: BATCH AS FIXED -----------------------------------------------------
# y ~ 1 + tissue_type + BATCH
metak <- "meta5"
tmp_metadata <- ei(sce_tmp)

design <- createDesignMatrix(tmp_metadata, cols_design = c("tissue_type", "BATCH"))
# colnames(design)
contrast <- createContrast(c(0,1,rep(0,6)))
# 0=intercept, 1=tissue_type, 0...= batch 2 to batch 7

# all markers
sel_markers <- rownames(sce_tmp) %in% rownames(sce_tmp)

C6C7_res_DS_2 <- diffcyt(sce_tmp,
                         clustering_to_use = metak,
                         analysis_type = "DS",
                         method_DS = "diffcyt-DS-limma",
                         design = design,
                         contrast = contrast,
                         markers_to_test = sel_markers,
                         verbose = FALSE,
                         transform = FALSE)

# OPTION 3: MIXED MODEL  -------------------------------------------------------
# y ~ tissue_type + BATCH + (1 | patient_id)
metak <- "meta5"
tmp_metadata <- ei(sce_tmp)

designF <- createFormula(tmp_metadata, cols_fixed = c("tissue_type", "BATCH"), cols_random = c("patient_id"))
# designF$formula
contrastF <- createContrast(c(1,0,0))
# 1=tissue_type, 0=batch, 0=random effect

# all markers
sel_markers <- rownames(sce_tmp) %in% rownames(sce_tmp)

C6C7_res_DS_3 <- diffcyt(sce_tmp,
                         clustering_to_use = metak,
                         analysis_type = "DS",
                         method_DS = "diffcyt-DS-LMM",  # diffcyt-DS-limma doesn't work; use diffcyt-DS_LMM
                         design = designF,
                         contrast = contrastF,
                         markers_to_test = sel_markers,
                         verbose = FALSE,
                         transform = FALSE)



# COMBINE RESULTS  -------------------------------------------------------------


################################## C1 ##########################################
C1_DS_1 <- as.data.frame(rowData(C1_res_DS_1$res))
C1_DS_2 <- as.data.frame(rowData(C1_res_DS_2$res))
C1_DS_3 <- as.data.frame(rowData(C1_res_DS_3$res))

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
# option 2
C1_DS_2 <- C1_DS_2 %>%
  select(cluster_id, marker_id, logFC, p_val, p_adj) %>%
  dplyr::rename(
    cluster = cluster_id,
    marker = marker_id,
    logFC_DS = logFC,
    p_val_DS = p_val,
    p_adj_DS = p_adj
  ) %>% mutate(
    model = "2: Tissue + Batch",
    tier1_cluster = "C1"
  )
# option 3
C1_DS_3 <- C1_DS_3 %>%
  select(cluster_id, marker_id, p_val, p_adj) %>%
  dplyr::rename(
    cluster = cluster_id,
    marker = marker_id,
    p_val_DS = p_val,
    p_adj_DS = p_adj
  ) %>% mutate(
    model = "3: Mixed",
    tier1_cluster = "C1"
  )

# C1_DS_all <- bind_rows(C1_DS_1, C1_DS_2, C1_DS_3)
C1_DS_all <- bind_rows(C1_DS_1, C1_DS_2)


################################## C2 ##########################################
C2_DS_1 <- as.data.frame(rowData(C2_res_DS_1$res))
C2_DS_2 <- as.data.frame(rowData(C2_res_DS_2$res))
C2_DS_3 <- as.data.frame(rowData(C2_res_DS_3$res))

# option 1
C2_DS_1 <- C2_DS_1 %>%
  select(cluster_id, marker_id, logFC, p_val, p_adj) %>%
  dplyr::rename(
    cluster = cluster_id,
    marker = marker_id,
    logFC_DS = logFC,
    p_val_DS = p_val,
    p_adj_DS = p_adj
  ) %>% mutate(
    model = "1: Tissue",
    tier1_cluster = "C2"
  )
# option 2
C2_DS_2 <- C2_DS_2 %>%
  select(cluster_id, marker_id, logFC, p_val, p_adj) %>%
  dplyr::rename(
    cluster = cluster_id,
    marker = marker_id,
    logFC_DS = logFC,
    p_val_DS = p_val,
    p_adj_DS = p_adj
  ) %>% mutate(
    model = "2: Tissue + Batch",
    tier1_cluster = "C2"
  )
# option 3
C2_DS_3 <- C2_DS_3 %>%
  select(cluster_id, marker_id, p_val, p_adj) %>%
  dplyr::rename(
    cluster = cluster_id,
    marker = marker_id,
    p_val_DS = p_val,
    p_adj_DS = p_adj
  ) %>% mutate(
    model = "3: Mixed",
    tier1_cluster = "C2"
  )

# C2_DS_all <- bind_rows(C2_DS_1, C2_DS_2, C2_DS_3)
C2_DS_all <- bind_rows(C2_DS_1, C2_DS_2)


################################## C3,C4 #######################################
C3C4_DS_1 <- as.data.frame(rowData(C3C4_res_DS_1$res))
C3C4_DS_2 <- as.data.frame(rowData(C3C4_res_DS_2$res))
C3C4_DS_3 <- as.data.frame(rowData(C3C4_res_DS_3$res))

# option 1
C3C4_DS_1 <- C3C4_DS_1 %>%
  select(cluster_id, marker_id, logFC, p_val, p_adj) %>%
  dplyr::rename(
    cluster = cluster_id,
    marker = marker_id,
    logFC_DS = logFC,
    p_val_DS = p_val,
    p_adj_DS = p_adj
  ) %>% mutate(
    model = "1: Tissue",
    tier1_cluster = "C3C4"
  )
# option 2
C3C4_DS_2 <- C3C4_DS_2 %>%
  select(cluster_id, marker_id, logFC, p_val, p_adj) %>%
  dplyr::rename(
    cluster = cluster_id,
    marker = marker_id,
    logFC_DS = logFC,
    p_val_DS = p_val,
    p_adj_DS = p_adj
  ) %>% mutate(
    model = "2: Tissue + Batch",
    tier1_cluster = "C3C4"
  )
# option 3
C3C4_DS_3 <- C3C4_DS_3 %>%
  select(cluster_id, marker_id, p_val, p_adj) %>%
  dplyr::rename(
    cluster = cluster_id,
    marker = marker_id,
    p_val_DS = p_val,
    p_adj_DS = p_adj
  ) %>% mutate(
    model = "3: Mixed",
    tier1_cluster = "C3C4"
  )

# C3C4_DS_all <- bind_rows(C3C4_DS_1, C3C4_DS_2, C3C4_DS_3)
C3C4_DS_all <- bind_rows(C3C4_DS_1, C3C4_DS_2)


################################## C5,C9 #######################################
C5C9_DS_1 <- as.data.frame(rowData(C5C9_res_DS_1$res))
C5C9_DS_2 <- as.data.frame(rowData(C5C9_res_DS_2$res))
C5C9_DS_3 <- as.data.frame(rowData(C5C9_res_DS_3$res))

# option 1
C5C9_DS_1 <- C5C9_DS_1 %>%
  select(cluster_id, marker_id, logFC, p_val, p_adj) %>%
  dplyr::rename(
    cluster = cluster_id,
    marker = marker_id,
    logFC_DS = logFC,
    p_val_DS = p_val,
    p_adj_DS = p_adj
  ) %>% mutate(
    model = "1: Tissue",
    tier1_cluster = "C5C9"
  )
# option 2
C5C9_DS_2 <- C5C9_DS_2 %>%
  select(cluster_id, marker_id, logFC, p_val, p_adj) %>%
  dplyr::rename(
    cluster = cluster_id,
    marker = marker_id,
    logFC_DS = logFC,
    p_val_DS = p_val,
    p_adj_DS = p_adj
  ) %>% mutate(
    model = "2: Tissue + Batch",
    tier1_cluster = "C5C9"
  )
# option 3
C5C9_DS_3 <- C5C9_DS_3 %>%
  select(cluster_id, marker_id, p_val, p_adj) %>%
  dplyr::rename(
    cluster = cluster_id,
    marker = marker_id,
    p_val_DS = p_val,
    p_adj_DS = p_adj
  ) %>% mutate(
    model = "3: Mixed",
    tier1_cluster = "C5C9"
  )

# C5C9_DS_all <- bind_rows(C5C9_DS_1, C5C9_DS_2, C5C9_DS_3)
C5C9_DS_all <- bind_rows(C5C9_DS_1, C5C9_DS_2)


################################## C6,C7 #######################################
C6C7_DS_1 <- as.data.frame(rowData(C6C7_res_DS_1$res))
C6C7_DS_2 <- as.data.frame(rowData(C6C7_res_DS_2$res))
C6C7_DS_3 <- as.data.frame(rowData(C6C7_res_DS_3$res))

# option 1
C6C7_DS_1 <- C6C7_DS_1 %>%
  select(cluster_id, marker_id, logFC, p_val, p_adj) %>%
  dplyr::rename(
    cluster = cluster_id,
    marker = marker_id,
    logFC_DS = logFC,
    p_val_DS = p_val,
    p_adj_DS = p_adj
  ) %>% mutate(
    model = "1: Tissue",
    tier1_cluster = "C6C7"
  )
# option 2
C6C7_DS_2 <- C6C7_DS_2 %>%
  select(cluster_id, marker_id, logFC, p_val, p_adj) %>%
  dplyr::rename(
    cluster = cluster_id,
    marker = marker_id,
    logFC_DS = logFC,
    p_val_DS = p_val,
    p_adj_DS = p_adj
  ) %>% mutate(
    model = "2: Tissue + Batch",
    tier1_cluster = "C6C7"
  )
# option 3
C6C7_DS_3 <- C6C7_DS_3 %>%
  select(cluster_id, marker_id, p_val, p_adj) %>%
  dplyr::rename(
    cluster = cluster_id,
    marker = marker_id,
    p_val_DS = p_val,
    p_adj_DS = p_adj
  ) %>% mutate(
    model = "3: Mixed",
    tier1_cluster = "C6C7"
  )

# C6C7_DS_all <- bind_rows(C6C7_DS_1, C6C7_DS_2, C6C7_DS_3)
C6C7_DS_all <- bind_rows(C6C7_DS_1, C6C7_DS_2)


# COMBINE ALL
DS_FINAL <- bind_rows(C1_DS_all, C2_DS_all, C3C4_DS_all, C5C9_DS_all, C6C7_DS_all)
rownames(DS_FINAL) <- NULL

file_name <- paste("TBNK_subclust_DS.txt")
write.table(DS_FINAL, file.path(res_dir, file_name), sep = "\t", row.names = FALSE, quote = FALSE)



# ------------------------------------------------------------------------------
# Bind only mixed OPTION 3 model
#   - OPTION 3 used: diffcyt-DS-LMM, instead of diffcyt-DS-limma because
#     # designF & contrastF are different acceptable type from limma

DS_FINAL_MIXED <- bind_rows(C1_DS_3, C2_DS_3, C3C4_DS_3, C5C9_DS_3, C6C7_DS_3)
rownames(DS_FINAL_MIXED) <- NULL

file_name <- paste("TBNK_subclust_DS_Mixed.txt")
write.table(DS_FINAL_MIXED, file.path(res_dir, file_name), sep = "\t", row.names = FALSE, quote = FALSE)

