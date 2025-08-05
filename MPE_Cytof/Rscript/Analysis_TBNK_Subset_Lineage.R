sel_panel <- "TBNK"

rds_subclust_dir <- paste0("D:/CHRISTOPHER_BUI/MPE_CYTOF_RDS/", sel_panel, "/Results_Subcluster/RDS_Subcluster")

rds_files <- list.files(rds_subclust_dir)

sce_c1 <- readRDS(file.path(rds_subclust_dir, rds_files[1]))
# sce_c2 <- readRDS(file.path(rds_subclust_dir, rds_files[2]))
# sce_c3c4 <- readRDS(file.path(rds_subclust_dir, rds_files[3]))
# sce_c5c9 <- readRDS(file.path(rds_subclust_dir, rds_files[4]))
# sce_c6c7 <- readRDS(file.path(rds_subclust_dir, rds_files[5]))
# sce_c8 <- readRDS(file.path(rds_subclust_dir, rds_files[6]))

output_dir <- paste0("D:/CHRISTOPHER_BUI/MPE_CYTOF_RDS/", sel_panel, "/Results_Subcluster")

# C1 ---------------------------------------------------------------------------
prefix <- "C1"

# DOTPLOTS
m3_dp <- dotplotTables(sce_c1, "meta3", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m5_dp <- dotplotTables(sce_c1, "meta5", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m7_dp <- dotplotTables(sce_c1, "meta7", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m10_dp <- dotplotTables(sce_c1, "meta10", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)

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

# meta7
tmp_fig <- dotplotFig(m7_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta7", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)

# meta10
tmp_fig <- dotplotFig(m10_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta10", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)


# CORR
corr_m3 <- cor(m3_dp$expr_wide)
corr_m5 <- cor(m5_dp$expr_wide)
corr_m7 <- cor(m7_dp$expr_wide)
corr_m10 <- cor(m10_dp$expr_wide)

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

pheatmap(corr_m7,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta7 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

pheatmap(corr_m10,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta10 clusters",
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

# meta5 & meta7
corr_mat_m5_m7 <- cor(m5_dp$expr_wide, m7_dp$expr_wide)
pheatmap(corr_mat_m5_m7,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta5 and meta7 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# meta7 meta10
corr_mat_m7_m10 <- cor(m7_dp$expr_wide, m10_dp$expr_wide)
pheatmap(corr_mat_m7_m10,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta7 and meta10 clusters",
         fontsize_row = 10,
         fontsize_col = 10)


# ABUNDANCE DOTPLOT
# meta3
CATALYST::plotAbundances(sce_c1, k = "meta3", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta5
CATALYST::plotAbundances(sce_c1, k = "meta5", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta7
CATALYST::plotAbundances(sce_c1, k = "meta7", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta10
CATALYST::plotAbundances(sce_c1, k = "meta10", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")


# PBMC BATCH SAMPLE
sce_c1_pbmc <- filterSCE(sce_c1, tissue_type == "PBMC")

# meta3
clust_tbl_pbmc <- table(cluster_ids(sce_c1_pbmc, "meta3"), colData(sce_c1_pbmc)$sample_id)
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

# meta10
clust_tbl_pbmc <- table(cluster_ids(sce_c1_pbmc, "meta10"), colData(sce_c1_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")



# C2 ---------------------------------------------------------------------------
prefix <- "C2"

# DOTPLOTS
m3_dp <- dotplotTables(sce_c2, "meta3", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m6_dp <- dotplotTables(sce_c2, "meta6", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m9_dp <- dotplotTables(sce_c2, "meta9", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)

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

# meta9
tmp_fig <- dotplotFig(m9_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta9", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)


# CORR
corr_m3 <- cor(m3_dp$expr_wide)
corr_m6 <- cor(m6_dp$expr_wide)
corr_m9 <- cor(m9_dp$expr_wide)

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

pheatmap(corr_m9,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta9 clusters",
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

# meta3 & meta9
corr_mat_m3_m9 <- cor(m3_dp$expr_wide, m9_dp$expr_wide)
pheatmap(corr_mat_m3_m9,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta9 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# meta6 & meta9
corr_mat_m6_m9 <- cor(m6_dp$expr_wide, m9_dp$expr_wide)
pheatmap(corr_mat_m6_m9,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta6 and meta9 clusters",
         fontsize_row = 10,
         fontsize_col = 10)


# ABUNDANCE DOTPLOT
# meta3
CATALYST::plotAbundances(sce_c2, k = "meta3", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta6
CATALYST::plotAbundances(sce_c2, k = "meta6", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta9
CATALYST::plotAbundances(sce_c2, k = "meta9", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")


# PBMC BATCH SAMPLE
sce_c2_pbmc <- filterSCE(sce_c2, tissue_type == "PBMC")

# meta3
clust_tbl_pbmc <- table(cluster_ids(sce_c2_pbmc, "meta3"), colData(sce_c2_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta6
clust_tbl_pbmc <- table(cluster_ids(sce_c2_pbmc, "meta6"), colData(sce_c2_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta9
clust_tbl_pbmc <- table(cluster_ids(sce_c2_pbmc, "meta9"), colData(sce_c2_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")



# C3,C4 ---------------------------------------------------------------------------
prefix <- "C3C4"

# DOTPLOTS
m3_dp <- dotplotTables(sce_c3c4, "meta3", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m6_dp <- dotplotTables(sce_c3c4, "meta6", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m12_dp <- dotplotTables(sce_c3c4, "meta12", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m4_dp <- dotplotTables(sce_c3c4, "meta4", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)

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

# meta12
tmp_fig <- dotplotFig(m12_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta12", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)

# meta4
tmp_fig <- dotplotFig(m4_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta4", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)



# CORR
corr_m3 <- cor(m3_dp$expr_wide)
corr_m6 <- cor(m6_dp$expr_wide)
corr_m12 <- cor(m12_dp$expr_wide)
corr_m4 <- cor(m4_dp$expr_wide)

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

pheatmap(corr_m12,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta12 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

pheatmap(corr_m4,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta4 clusters",
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

# meta3 & meta4
corr_mat_m3_m4 <- cor(m3_dp$expr_wide, m4_dp$expr_wide)
pheatmap(corr_mat_m3_m4,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta4 clusters",
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




# ABUNDANCE DOTPLOT
# meta3
CATALYST::plotAbundances(sce_c3c4, k = "meta3", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta6
CATALYST::plotAbundances(sce_c3c4, k = "meta6", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta12
CATALYST::plotAbundances(sce_c3c4, k = "meta12", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta4
CATALYST::plotAbundances(sce_c3c4, k = "meta4", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")


# PBMC BATCH SAMPLE
sce_c3c4_pbmc <- filterSCE(sce_c3c4, tissue_type == "PBMC")

# meta3
clust_tbl_pbmc <- table(cluster_ids(sce_c3c4_pbmc, "meta3"), colData(sce_c3c4_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta6
clust_tbl_pbmc <- table(cluster_ids(sce_c3c4_pbmc, "meta6"), colData(sce_c3c4_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta12
clust_tbl_pbmc <- table(cluster_ids(sce_c3c4_pbmc, "meta12"), colData(sce_c3c4_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta4
clust_tbl_pbmc <- table(cluster_ids(sce_c3c4_pbmc, "meta4"), colData(sce_c3c4_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")



# C5,C9 ---------------------------------------------------------------------------
prefix <- "C5C9"

# DOTPLOTS
m3_dp <- dotplotTables(sce_c5c9, "meta3", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m4_dp <- dotplotTables(sce_c5c9, "meta4", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m5_dp <- dotplotTables(sce_c5c9, "meta5", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m6_dp <- dotplotTables(sce_c5c9, "meta6", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)

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

# meta6
tmp_fig <- dotplotFig(m6_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta6", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)



# CORR
corr_m3 <- cor(m3_dp$expr_wide)
corr_m4 <- cor(m4_dp$expr_wide)
corr_m5 <- cor(m5_dp$expr_wide)
corr_m6 <- cor(m6_dp$expr_wide)

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

pheatmap(corr_m6,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta6 clusters",
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

# meta5 & meta6
corr_mat_m5_m6 <- cor(m5_dp$expr_wide, m6_dp$expr_wide)
pheatmap(corr_mat_m5_m6,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta5 and meta6 clusters",
         fontsize_row = 10,
         fontsize_col = 10)


# ABUNDANCE DOTPLOT
# meta3
CATALYST::plotAbundances(sce_c5c9, k = "meta3", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta4
CATALYST::plotAbundances(sce_c5c9, k = "meta4", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta5
CATALYST::plotAbundances(sce_c5c9, k = "meta5", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta6
CATALYST::plotAbundances(sce_c5c9, k = "meta6", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")



# PBMC BATCH SAMPLE
sce_c5c9_pbmc <- filterSCE(sce_c5c9, tissue_type == "PBMC")

# meta3
clust_tbl_pbmc <- table(cluster_ids(sce_c5c9_pbmc, "meta3"), colData(sce_c5c9_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta4
clust_tbl_pbmc <- table(cluster_ids(sce_c5c9_pbmc, "meta4"), colData(sce_c5c9_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta5
clust_tbl_pbmc <- table(cluster_ids(sce_c5c9_pbmc, "meta5"), colData(sce_c5c9_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta6
clust_tbl_pbmc <- table(cluster_ids(sce_c5c9_pbmc, "meta6"), colData(sce_c5c9_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")



# C6,C7 ---------------------------------------------------------------------------
prefix <- "C6C7"

# DOTPLOTS
m3_dp <- dotplotTables(sce_c6c7, "meta3", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m5_dp <- dotplotTables(sce_c6c7, "meta5", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m6_dp <- dotplotTables(sce_c6c7, "meta6", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m7_dp <- dotplotTables(sce_c6c7, "meta7", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)

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

# meta6
tmp_fig <- dotplotFig(m6_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta6", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)

# meta7
tmp_fig <- dotplotFig(m7_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta7", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)



# CORR
corr_m3 <- cor(m3_dp$expr_wide)
corr_m5 <- cor(m5_dp$expr_wide)
corr_m6 <- cor(m6_dp$expr_wide)
corr_m7 <- cor(m7_dp$expr_wide)

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

pheatmap(corr_m6,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta6 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

pheatmap(corr_m7,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta7 clusters",
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

# meta3 & meta6
corr_mat_m3_m6 <- cor(m3_dp$expr_wide, m6_dp$expr_wide)
pheatmap(corr_mat_m3_m6,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta6 clusters",
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

# meta5 & meta6
corr_mat_m5_m6 <- cor(m5_dp$expr_wide, m6_dp$expr_wide)
pheatmap(corr_mat_m5_m6,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta5 and meta6 clusters",
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

# meta6 & meta7
corr_mat_m6_m7 <- cor(m6_dp$expr_wide, m7_dp$expr_wide)
pheatmap(corr_mat_m6_m7,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta6 and meta7 clusters",
         fontsize_row = 10,
         fontsize_col = 10)


# ABUNDANCE DOTPLOT
# meta3
CATALYST::plotAbundances(sce_c6c7, k = "meta3", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta5
CATALYST::plotAbundances(sce_c6c7, k = "meta5", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta6
CATALYST::plotAbundances(sce_c6c7, k = "meta6", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta7
CATALYST::plotAbundances(sce_c6c7, k = "meta7", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")



# PBMC BATCH SAMPLE
sce_c6c7_pbmc <- filterSCE(sce_c6c7, tissue_type == "PBMC")

# meta3
clust_tbl_pbmc <- table(cluster_ids(sce_c6c7_pbmc, "meta3"), colData(sce_c6c7_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta5
clust_tbl_pbmc <- table(cluster_ids(sce_c6c7_pbmc, "meta5"), colData(sce_c6c7_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta6
clust_tbl_pbmc <- table(cluster_ids(sce_c6c7_pbmc, "meta6"), colData(sce_c6c7_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta7
clust_tbl_pbmc <- table(cluster_ids(sce_c6c7_pbmc, "meta7"), colData(sce_c6c7_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")


# C8 ---------------------------------------------------------------------------
prefix <- "C8"

# DOTPLOTS
m4_dp <- dotplotTables(sce_c8, "meta4", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m5_dp <- dotplotTables(sce_c8, "meta5", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)

sel_aggregate <- c("median")
scale_option <- FALSE
label_for_dotplot <- paste(ifelse(scale_option, "scaled", "non-scaled"), sel_aggregate, "exprs", sep = " ")

# meta4
tmp_fig <- dotplotFig(m4_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta4", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)

# meta5
tmp_fig <- dotplotFig(m5_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta5", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)


# CORR
corr_m4 <- cor(m4_dp$expr_wide)
corr_m5 <- cor(m5_dp$expr_wide)

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


# meta3 & meta5
corr_mat_m3_m5 <- cor(m3_dp$expr_wide, m5_dp$expr_wide)
pheatmap(corr_mat_m3_m5,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta5 clusters",
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

# meta3 & meta7
corr_mat_m3_m7 <- cor(m3_dp$expr_wide, m7_dp$expr_wide)
pheatmap(corr_mat_m3_m7,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta7 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# meta5 & meta6
corr_mat_m5_m6 <- cor(m5_dp$expr_wide, m6_dp$expr_wide)
pheatmap(corr_mat_m5_m6,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta5 and meta6 clusters",
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

# meta6 & meta7
corr_mat_m6_m7 <- cor(m6_dp$expr_wide, m7_dp$expr_wide)
pheatmap(corr_mat_m6_m7,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta6 and meta7 clusters",
         fontsize_row = 10,
         fontsize_col = 10)


# ABUNDANCE DOTPLOT
# meta3
CATALYST::plotAbundances(sce_c6c7, k = "meta3", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta5
CATALYST::plotAbundances(sce_c6c7, k = "meta5", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta6
CATALYST::plotAbundances(sce_c6c7, k = "meta6", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta7
CATALYST::plotAbundances(sce_c6c7, k = "meta7", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")



# PBMC BATCH SAMPLE
sce_c6c7_pbmc <- filterSCE(sce_c6c7, tissue_type == "PBMC")

# meta3
clust_tbl_pbmc <- table(cluster_ids(sce_c6c7_pbmc, "meta3"), colData(sce_c6c7_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta5
clust_tbl_pbmc <- table(cluster_ids(sce_c6c7_pbmc, "meta5"), colData(sce_c6c7_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta6
clust_tbl_pbmc <- table(cluster_ids(sce_c6c7_pbmc, "meta6"), colData(sce_c6c7_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

# meta7
clust_tbl_pbmc <- table(cluster_ids(sce_c6c7_pbmc, "meta7"), colData(sce_c6c7_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")

