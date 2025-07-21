sel_panel <- "Myeloid"

rds_subclust_dir <- paste0("D:/CHRISTOPHER_BUI/MPE_CYTOF_RDS/", sel_panel, "/Results_Subcluster/RDS_Subcluster")


rds_files <- list.files(rds_subclust_dir)
# rds files
sce_c1 <- readRDS(file.path(rds_subclust_dir, "C1_sce_subclust_full_lineage.rds"))
# sce_c234 <- readRDS(file.path(rds_subclust_dir, "C2C3C4_sce_subclust_full_lineage.rds"))
# sce_c5 <- readRDS(file.path(rds_subclust_dir, "C5_sce_subclust_full_lineage.rds"))
# sce_c6 <- readRDS(file.path(rds_subclust_dir, "C6_sce_subclust_full_lineage.rds"))


output_dir <- paste0("D:/CHRISTOPHER_BUI/MPE_CYTOF_RDS/", sel_panel, "/Results_Subcluster")

#-------------------------------------------------------------------------------
# NOTE:
#   - REMOVED EPCAM & C5 FROM original TIER 1
#-------------------------------------------------------------------------------

# C1 ---------------------------------------------------------------------------
prefix <- "C1"

# DOTPLOTS
m3_dp <- dotplotTables(sce_c1, "meta3", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m4_dp <- dotplotTables(sce_c1, "meta4", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m5_dp <- dotplotTables(sce_c1, "meta5", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m7_dp <- dotplotTables(sce_c1, "meta7", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)

sel_aggregate <- c("median")
scale_option <- FALSE
label_for_dotplot <- paste(ifelse(scale_option, "scaled", "non-scaled"), sel_aggregate, "exprs", sep = " ")

# meta3
tmp_fig <- dotplotFig(m3_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta3", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)

# meta4
tmp_fig <- dotplotFig(m4_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta4", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)

# meta5
tmp_fig <- dotplotFig(m5_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta5", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)

# meta7
tmp_fig <- dotplotFig(m7_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta7", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)


# CORR
corr_m3 <- cor(m3_dp$expr_wide)
corr_m4 <- cor(m4_dp$expr_wide)
corr_m5 <- cor(m5_dp$expr_wide)
corr_m7 <- cor(m7_dp$expr_wide)

pheatmap(corr_m3,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta3 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

pheatmap(corr_m4,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta4 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

pheatmap(corr_m5,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta5 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

pheatmap(corr_m7,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta7 clusters",
         fontsize_row = 10,
         fontsize_col = 10)


# meta3 & meta4
corr_mat_m3_m4 <- cor(m3_dp$expr_wide, m4_dp$expr_wide)
pheatmap(corr_mat_m3_m4,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta4 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# meta3 & meta5
corr_mat_m3_m5 <- cor(m3_dp$expr_wide, m5_dp$expr_wide)
pheatmap(corr_mat_m3_m5,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta5 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# meta3 & meta7
corr_mat_m3_m7 <- cor(m3_dp$expr_wide, m7_dp$expr_wide)
pheatmap(corr_mat_m3_m7,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta7 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# meta4 & meta7
corr_mat_m4_m7 <- cor(m4_dp$expr_wide, m7_dp$expr_wide)
pheatmap(corr_mat_m4_m7,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta4 and meta7 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# meta5 & meta7
corr_mat_m5_m7 <- cor(m5_dp$expr_wide, m7_dp$expr_wide)
pheatmap(corr_mat_m5_m7,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta5 and meta7 clusters",
         fontsize_row = 10,
         fontsize_col = 10)


# ABUNDANCE DOTPLOT
# meta3
CATALYST::plotAbundances(sce_c1, k = "meta3", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta4
CATALYST::plotAbundances(sce_c1, k = "meta4", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta5
CATALYST::plotAbundances(sce_c1, k = "meta5", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta7
CATALYST::plotAbundances(sce_c1, k = "meta7", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")


# PBMC BATCH SAMPLE
sce_c1_pbmc <- filterSCE(sce_c1, tissue_type == "PBMC")

# meta2
clust_tbl_pbmc <- table(cluster_ids(sce_c1_pbmc, "meta2"), colData(sce_c1_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta3
clust_tbl_pbmc <- table(cluster_ids(sce_c1_pbmc, "meta3"), colData(sce_c1_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta4
clust_tbl_pbmc <- table(cluster_ids(sce_c1_pbmc, "meta4"), colData(sce_c1_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta5
clust_tbl_pbmc <- table(cluster_ids(sce_c1_pbmc, "meta5"), colData(sce_c1_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta7
clust_tbl_pbmc <- table(cluster_ids(sce_c1_pbmc, "meta7"), colData(sce_c1_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")


# C2,C3,C4 ---------------------------------------------------------------------
prefix <- "C2C3C4"

# DOTPLOTS
m3_dp <- dotplotTables(sce_c234, "meta3", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m5_dp <- dotplotTables(sce_c234, "meta5", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m8_dp <- dotplotTables(sce_c234, "meta8", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)

sel_aggregate <- c("median")
scale_option <- FALSE
label_for_dotplot <- paste(ifelse(scale_option, "scaled", "non-scaled"), sel_aggregate, "exprs", sep = " ")

# meta3
tmp_fig <- dotplotFig(m3_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta3", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)

# meta5
tmp_fig <- dotplotFig(m5_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta5", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)

# meta8
tmp_fig <- dotplotFig(m8_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta8", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)


# CORR
corr_m3 <- cor(m3_dp$expr_wide)
corr_m5 <- cor(m5_dp$expr_wide)
corr_m8 <- cor(m8_dp$expr_wide)


pheatmap(corr_m3,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta3 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

pheatmap(corr_m5,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta5 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

pheatmap(corr_m8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta8 clusters",
         fontsize_row = 10,
         fontsize_col = 10)


# meta3 & meta5
corr_mat_m3_m5 <- cor(m3_dp$expr_wide, m5_dp$expr_wide)
pheatmap(corr_mat_m3_m5,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta5 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# meta3 & meta8
corr_mat_m3_m8 <- cor(m3_dp$expr_wide, m8_dp$expr_wide)
pheatmap(corr_mat_m3_m8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta8 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# meta5 & meta8
corr_mat_m5_m8 <- cor(m5_dp$expr_wide, m8_dp$expr_wide)
pheatmap(corr_mat_m5_m8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta5 and meta8 clusters",
         fontsize_row = 10,
         fontsize_col = 10)


# ABUNDANCE DOTPLOT
# meta3
CATALYST::plotAbundances(sce_c234, k = "meta3", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta5
CATALYST::plotAbundances(sce_c234, k = "meta5", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta8
CATALYST::plotAbundances(sce_c234, k = "meta8", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")


# PBMC BATCH SAMPLE
sce_c234_pbmc <- filterSCE(sce_c234, tissue_type == "PBMC")

# meta3
clust_tbl_pbmc <- table(cluster_ids(sce_c234_pbmc, "meta3"), colData(sce_c234_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta5
clust_tbl_pbmc <- table(cluster_ids(sce_c234_pbmc, "meta5"), colData(sce_c234_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta8
clust_tbl_pbmc <- table(cluster_ids(sce_c234_pbmc, "meta8"), colData(sce_c234_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")


# C5 ---------------------------------------------------------------------------
prefix <- "C5"

# DOTPLOTS
m3_dp <- dotplotTables(sce_c5, "meta3", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m6_dp <- dotplotTables(sce_c5, "meta6", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m8_dp <- dotplotTables(sce_c5, "meta8", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)

sel_aggregate <- c("median")
scale_option <- FALSE
label_for_dotplot <- paste(ifelse(scale_option, "scaled", "non-scaled"), sel_aggregate, "exprs", sep = " ")

# meta3
tmp_fig <- dotplotFig(m3_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta3", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)

# meta6
tmp_fig <- dotplotFig(m6_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta6", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)

# meta8
tmp_fig <- dotplotFig(m8_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta8", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)


# CORR
corr_m3 <- cor(m3_dp$expr_wide)
corr_m6 <- cor(m6_dp$expr_wide)
corr_m8 <- cor(m8_dp$expr_wide)


pheatmap(corr_m3,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta3 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

pheatmap(corr_m6,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta6 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

pheatmap(corr_m8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta8 clusters",
         fontsize_row = 10,
         fontsize_col = 10)


# meta3 & meta6
corr_mat_m3_m6 <- cor(m3_dp$expr_wide, m6_dp$expr_wide)
pheatmap(corr_mat_m3_m6,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta6 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# meta3 & meta8
corr_mat_m3_m8 <- cor(m3_dp$expr_wide, m8_dp$expr_wide)
pheatmap(corr_mat_m3_m8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta8 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# meta6 & meta8
corr_mat_m6_m8 <- cor(m6_dp$expr_wide, m8_dp$expr_wide)
pheatmap(corr_mat_m6_m8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta6 and meta8 clusters",
         fontsize_row = 10,
         fontsize_col = 10)


# ABUNDANCE DOTPLOT
# meta3
CATALYST::plotAbundances(sce_c5, k = "meta3", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta5
CATALYST::plotAbundances(sce_c5, k = "meta6", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta8
CATALYST::plotAbundances(sce_c5, k = "meta8", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")


# PBMC BATCH SAMPLE
sce_c5_pbmc <- filterSCE(sce_c5, tissue_type == "PBMC")

# meta3
clust_tbl_pbmc <- table(cluster_ids(sce_c5_pbmc, "meta3"), colData(sce_c5_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta6
clust_tbl_pbmc <- table(cluster_ids(sce_c5_pbmc, "meta6"), colData(sce_c5_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta8
clust_tbl_pbmc <- table(cluster_ids(sce_c5_pbmc, "meta8"), colData(sce_c5_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")


# C6 ---------------------------------------------------------------------------
prefix <- "C6"

# DOTPLOTS
m4_dp <- dotplotTables(sce_c6, "meta4", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m6_dp <- dotplotTables(sce_c6, "meta6", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m8_dp <- dotplotTables(sce_c6, "meta8", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)

sel_aggregate <- c("median")
scale_option <- FALSE
label_for_dotplot <- paste(ifelse(scale_option, "scaled", "non-scaled"), sel_aggregate, "exprs", sep = " ")

# meta4
tmp_fig <- dotplotFig(m4_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta4", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)

# meta6
tmp_fig <- dotplotFig(m6_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta6", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)

# meta8
tmp_fig <- dotplotFig(m8_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta8", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)


# CORR
corr_m4 <- cor(m4_dp$expr_wide)
corr_m6 <- cor(m6_dp$expr_wide)
corr_m8 <- cor(m8_dp$expr_wide)


pheatmap(corr_m4,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta4 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

pheatmap(corr_m6,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta6 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

pheatmap(corr_m8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta8 clusters",
         fontsize_row = 10,
         fontsize_col = 10)


# meta4 & meta6
corr_mat_m4_m6 <- cor(m4_dp$expr_wide, m6_dp$expr_wide)
pheatmap(corr_mat_m4_m6,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta4 and meta6 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# meta4 & meta8
corr_mat_m4_m8 <- cor(m4_dp$expr_wide, m8_dp$expr_wide)
pheatmap(corr_mat_m4_m8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta4 and meta8 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# meta6 & meta8
corr_mat_m6_m8 <- cor(m6_dp$expr_wide, m8_dp$expr_wide)
pheatmap(corr_mat_m6_m8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta6 and meta8 clusters",
         fontsize_row = 10,
         fontsize_col = 10)


# ABUNDANCE DOTPLOT
# meta4
CATALYST::plotAbundances(sce_c6, k = "meta4", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta6
CATALYST::plotAbundances(sce_c6, k = "meta6", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta8
CATALYST::plotAbundances(sce_c6, k = "meta8", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")


# PBMC BATCH SAMPLE
sce_c6_pbmc <- filterSCE(sce_c6, tissue_type == "PBMC")

# meta4
clust_tbl_pbmc <- table(cluster_ids(sce_c6_pbmc, "meta4"), colData(sce_c6_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta6
clust_tbl_pbmc <- table(cluster_ids(sce_c6_pbmc, "meta6"), colData(sce_c6_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta8
clust_tbl_pbmc <- table(cluster_ids(sce_c6_pbmc, "meta8"), colData(sce_c6_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

