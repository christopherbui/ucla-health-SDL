sel_panel <- "TBNK"

rds_subclust_dir <- paste0("D:/CHRISTOPHER_BUI/MPE_CYTOF_RDS/", sel_panel, "/Results_Subcluster/RDS_Subcluster")

rds_files <- list.files(rds_subclust_dir)

# sce_c1 <- readRDS(file.path(rds_subclust_dir, rds_files[1]))
sce_c2 <- readRDS(file.path(rds_subclust_dir, rds_files[2]))
# sce_c3c4 <- readRDS(file.path(rds_subclust_dir, rds_files[3]))
# sce_c5c9 <- readRDS(file.path(rds_subclust_dir, rds_files[4]))
# sce_c6c7 <- readRDS(file.path(rds_subclust_dir, rds_files[5]))
# sce_c8 <- readRDS(file.path(rds_subclust_dir, rds_files[6]))

sce_main <- readRDS("D:/CHRISTOPHER_BUI/MPE_CYTOF_RDS/TBNK/sce_subset_lineage.rds")

# sce with all markers
# 1. remove outliers from meta9 C1
# 2. project sce_c1 metadata to sce_main
sce_main_c1 <- sce_main[, sce_main$meta9 == 1]

umap_matrix <- reducedDim(sce_main_c1, "UMAP")
tmp_sel <- umap_matrix[, 1] < (-10) | umap_matrix[, 2] > 10
colData(sce_main)$c1_outlier <- FALSE
colData(sce_main_c1)$c1_outlier[colData(sce_main)$meta9 == 1] <- tmp_sel

sce_main <- sce_main[, !(sce_main$meta9 == 1 & colData(sce_main)$c1_outlier)]

################################## C2 ##########################################
metak <- "meta6"

# select SCE
sce_tmp <- sce_c2
sce_tmp$meta6 <- cluster_ids(sce_tmp, "meta6")

# set PBMC as reference
colData(sce_tmp)$tissue_type <- factor(colData(sce_tmp)$tissue_type, levels = c("PBMC", "MPE"))
tmp_metadata <- ei(sce_tmp)

# set batch 2 & 3 as reference
tmp_metadata$BATCH_NEW <- as.character(tmp_metadata$BATCH)
tmp_metadata$BATCH_NEW[tmp_metadata$BATCH_NEW %in% c("2", "3")] <- "2_3"  # combine 2 and 3
tmp_metadata$BATCH_NEW <- factor(tmp_metadata$BATCH_NEW, levels = c("2_3", "1", "4", "5", "6", "7"))

# select markers
sel_markers <- rownames(sce_tmp) %in% rownames(sce_tmp)

# OPTION 1: TISSUE ONLY -----------------------------------------------------
design <- createDesignMatrix(tmp_metadata, cols_design = c("tissue_type"))
contrast <- createContrast(c(0,1))

res_DS_1 <- diffcyt(sce_tmp,
                    clustering_to_use = metak,
                    analysis_type = "DS",
                    method_DS = "diffcyt-DS-limma",
                    design = design,
                    contrast = contrast,
                    verbose = FALSE,
                    transform = FALSE)



# OPTION 2: BATCH AS FIXED -----------------------------------------------------
design <- createDesignMatrix(tmp_metadata, cols_design = c("tissue_type", "BATCH_NEW"))
colnames(design)
contrast <- createContrast(c(0, 1, rep(0, 5)))
# intercept, tissue_type, batch 1, batch4, batch5, batch6, batch7

res_DS_2 <- diffcyt(sce_tmp,
                    clustering_to_use = metak,
                    analysis_type = "DS",
                    method_DS = "diffcyt-DS-limma",
                    design = design,
                    contrast = contrast,
                    verbose = FALSE,
                    transform = FALSE)
