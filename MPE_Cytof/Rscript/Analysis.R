rds_subclust_dir <- c("D:/CHRISTOPHER_BUI/MPE_CYTOF_RDS/Cytokines/Results_Subcluster/RDS_Subcluster")

rds_files <- list.files(rds_subclust_dir)

sce_c4 <- readRDS(file.path(rds_subclust_dir, rds_files[1]))
sce_c5 <- readRDS(file.path(rds_subclust_dir, rds_files[2]))
sce_c6 <- readRDS(file.path(rds_subclust_dir, rds_files[3]))
sce_myl <- readRDS(file.path(rds_subclust_dir, rds_files[4]))


# C4 ---------------------------------------------------------------------------
output_dir <- "D:/CHRISTOPHER_BUI/MPE_CYTOF_RDS/Cytokines/Results_Subcluster"
prefix <- "C4"

m3_dp <- dotplotTables(sce_c4, "meta3", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m6_dp <- dotplotTables(sce_c4, "meta6", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m4_dp <- dotplotTables(sce_c4, "meta4", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m8_dp <- dotplotTables(sce_c4, "meta8", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)

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

corr_m3 <- cor(m3_dp$expr_wide)
corr_m6 <- cor(m6_dp$expr_wide)
corr_m4 <- cor(m4_dp$expr_wide)
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

pheatmap(corr_m4,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta4 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

pheatmap(corr_m8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta8 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# corr between metacluster groups
# meta3 & meta6
corr_mat_c3_c6 <- cor(m3_dp$expr_wide, m6_dp$expr_wide)
pheatmap(corr_mat_c3_c6,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta6 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# meta3 & meta8
corr_mat_c3_c8 <- cor(m3_dp$expr_wide, m8_dp$expr_wide)
pheatmap(corr_mat_c3_c8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta8 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# meta6 & meta8
corr_mat_c6_c8 <- cor(m6_dp$expr_wide, m8_dp$expr_wide)
pheatmap(corr_mat_c6_c8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta6 and meta8 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# abundance boxplot
# meta3
CATALYST::plotAbundances(sce_c4, k = "meta3", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta6
CATALYST::plotAbundances(sce_c4, k = "meta6", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta8
CATALYST::plotAbundances(sce_c4, k = "meta8", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")

# view batch effect
sce_c4_pbmc <- filterSCE(sce_c4, tissue_type == "PBMC")

# meta3
clust_tbl_pbmc <- table(cluster_ids(sce_c4_pbmc, "meta3"), colData(sce_c4_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
png(file.path(output_dir, paste0("Cytokines_", prefix, "_meta3", "_pbmc_sample_cluster_heatmap.png")), width = 800, height = 600)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")
dev.off()

# meta6
clust_tbl_pbmc <- table(cluster_ids(sce_c4_pbmc, "meta6"), colData(sce_c4_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
png(file.path(output_dir, paste0("Cytokines_", prefix, "_meta6", "_pbmc_sample_cluster_heatmap.png")), width = 800, height = 600)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")
dev.off()

# meta8
clust_tbl_pbmc <- table(cluster_ids(sce_c4_pbmc, "meta8"), colData(sce_c4_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
png(file.path(output_dir, paste0("Cytokines_", prefix, "_meta8", "_pbmc_sample_cluster_heatmap.png")), width = 800, height = 600)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")
dev.off()


# C5 ---------------------------------------------------------------------------
output_dir <- "D:/CHRISTOPHER_BUI/MPE_CYTOF_RDS/Cytokines/Results_Subcluster"
prefix <- "C5"

m3_dp <- dotplotTables(sce_c5, "meta3", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m5_dp <- dotplotTables(sce_c5, "meta5", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m7_dp <- dotplotTables(sce_c5, "meta7", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)

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

corr_m3 <- cor(m3_dp$expr_wide)
corr_m5 <- cor(m5_dp$expr_wide)
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

pheatmap(corr_m7,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta7 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

corr_mat_c3_c5 <- cor(m3_dp$expr_wide, m5_dp$expr_wide)
pheatmap(corr_mat_c3_c5,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta5 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

corr_mat_c3_c7 <- cor(m3_dp$expr_wide, m7_dp$expr_wide)
pheatmap(corr_mat_c3_c7,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta7 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

corr_mat_c5_c7 <- cor(m5_dp$expr_wide, m7_dp$expr_wide)
pheatmap(corr_mat_c5_c7,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta5 and meta7 clusters",
         fontsize_row = 10,
         fontsize_col = 10)


# abundance boxplot
# meta3
CATALYST::plotAbundances(sce_c5, k = "meta3", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta5
CATALYST::plotAbundances(sce_c5, k = "meta5", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta7
CATALYST::plotAbundances(sce_c5, k = "meta7", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")


# view batch effect
sce_c5_pbmc <- filterSCE(sce_c5, tissue_type == "PBMC")

# meta3
clust_tbl_pbmc <- table(cluster_ids(sce_c5_pbmc, "meta3"), colData(sce_c5_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
png(file.path(output_dir, paste0("Cytokines_", prefix, "_meta3", "_pbmc_batch_heatmap.png")), width = 800, height = 600)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")
dev.off()

# meta5
clust_tbl_pbmc <- table(cluster_ids(sce_c5_pbmc, "meta5"), colData(sce_c5_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
png(file.path(output_dir, paste0("Cytokines_", prefix, "_meta5", "_pbmc_batch_heatmap.png")), width = 800, height = 600)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")
dev.off()

# meta7
clust_tbl_pbmc <- table(cluster_ids(sce_c5_pbmc, "meta7"), colData(sce_c5_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
png(file.path(output_dir, paste0("Cytokines_", prefix, "_meta7", "_pbmc_batch_heatmap.png")), width = 800, height = 600)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")
dev.off()



# C6 ---------------------------------------------------------------------------
prefix <- c("C6")

m3_dp <- dotplotTables(sce_c6, "meta3", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m7_dp <- dotplotTables(sce_c6, "meta7", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)

sel_aggregate <- c("median")
scale_option <- FALSE
label_for_dotplot <- paste(ifelse(scale_option, "scaled", "non-scaled"), sel_aggregate, "exprs", sep = " ")

# meta3
tmp_fig <- dotplotFig(m3_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta3", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)

# meta7
tmp_fig <- dotplotFig(m7_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta7", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)


corr_m3 <- cor(m3_dp$expr_wide)
corr_m7 <- cor(m7_dp$expr_wide)

pheatmap(corr_m3,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta3 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

pheatmap(corr_m7,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta7 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# corr between metacluster groups
corr_mat_c3_c7 <- cor(m3_dp$expr_wide, m7_dp$expr_wide)
pheatmap(corr_mat_c3_c7,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta7 clusters",
         fontsize_row = 10,
         fontsize_col = 10)



# abundance boxplot
# meta3
CATALYST::plotAbundances(sce_c6, k = "meta3", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")

# meta7
CATALYST::plotAbundances(sce_c6, k = "meta7", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")



# view batch effect
sce_c6_pbmc <- filterSCE(sce_c6, tissue_type == "PBMC")

# meta3
clust_tbl_pbmc <- table(cluster_ids(sce_c6_pbmc, "meta3"), colData(sce_c6_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
png(file.path(output_dir, paste0("Cytokines_", prefix, "_meta3", "_pbmc_batch_heatmap.png")), width = 800, height = 600)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")
dev.off()

# meta7
clust_tbl_pbmc <- table(cluster_ids(sce_c6_pbmc, "meta7"), colData(sce_c6_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
png(file.path(output_dir, paste0("Cytokines_", prefix, "_meta7", "_pbmc_batch_heatmap.png")), width = 800, height = 600)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")
dev.off()


# Myeloid ----------------------------------------------------------------------
prefix <- c("Myeloid")

m3_dp <- dotplotTables(sce_myl, "meta3", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m8_dp <- dotplotTables(sce_myl, "meta8", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m11_dp <- dotplotTables(sce_myl, "meta11", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)

sel_aggregate <- c("median")
scale_option <- FALSE
label_for_dotplot <- paste(ifelse(scale_option, "scaled", "non-scaled"), sel_aggregate, "exprs", sep = " ")

# meta3
tmp_fig <- dotplotFig(m3_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta3", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)

# meta8
tmp_fig <- dotplotFig(m8_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta8", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)

# meta11
tmp_fig <- dotplotFig(m11_dp$expr_long, lab = label_for_dotplot)
file_name <- paste(sel_panel, prefix, "meta11", "dotplot", sel_aggregate, ifelse(scale_option, "SCALED.png", "NON_SCALED.png"), sep = "_")
ggsave(filename = file.path(output_dir, file_name), plot = tmp_fig, width = 10, height = 8)


corr_m3 <- cor(m3_dp$expr_wide)
corr_m8 <- cor(m8_dp$expr_wide)
corr_m11 <- cor(m11_dp$expr_wide)

pheatmap(corr_m3,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta3 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

pheatmap(corr_m8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta8 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

pheatmap(corr_m11,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta11 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# corr between metacluster groups
# meta3 & meta8
corr_mat_c3_c8 <- cor(m3_dp$expr_wide, m8_dp$expr_wide)
pheatmap(corr_mat_c3_c8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta8 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# meta3 & meta11
corr_mat_c3_c11 <- cor(m3_dp$expr_wide, m11_dp$expr_wide)
pheatmap(corr_mat_c3_c11,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta11 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# meta8 & meta11
corr_mat_c8_c11 <- cor(m8_dp$expr_wide, m11_dp$expr_wide)
pheatmap(corr_mat_c8_c11,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta8 and meta11 clusters",
         fontsize_row = 10,
         fontsize_col = 10)



# abundance boxplot
# meta3
CATALYST::plotAbundances(sce_myl, k = "meta3", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")

CATALYST::plotAbundances(sce_myl, k = "meta8", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")

CATALYST::plotAbundances(sce_myl, k = "meta11", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")



# view batch effect
sce_myl_pbmc <- filterSCE(sce_myl, tissue_type == "PBMC")

# meta3
clust_tbl_pbmc <- table(cluster_ids(sce_myl_pbmc, "meta3"), colData(sce_myl_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
png(file.path(output_dir, paste0("Cytokines_", prefix, "_meta3", "_pbmc_batch_heatmap.png")), width = 800, height = 600)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")
dev.off()

# meta8
clust_tbl_pbmc <- table(cluster_ids(sce_myl_pbmc, "meta8"), colData(sce_myl_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
png(file.path(output_dir, paste0("Cytokines_", prefix, "_meta8", "_pbmc_batch_heatmap.png")), width = 800, height = 600)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")
dev.off()

# meta11
clust_tbl_pbmc <- table(cluster_ids(sce_myl_pbmc, "meta11"), colData(sce_myl_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
png(file.path(output_dir, paste0("Cytokines_", prefix, "_meta11", "_pbmc_batch_heatmap.png")), width = 800, height = 600)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")
dev.off()
