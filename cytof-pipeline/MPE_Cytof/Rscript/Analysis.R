rds_subclust_dir <- c("D:/CHRISTOPHER_BUI/MPE_CYTOF_RDS/Cytokines/Results_Subcluster/RDS_Subcluster")

rds_files <- list.files(rds_subclust_dir)

sce_c4 <- readRDS(file.path(rds_subclust_dir, rds_files[1]))
sce_c5 <- readRDS(file.path(rds_subclust_dir, rds_files[2]))
sce_c6 <- readRDS(file.path(rds_subclust_dir, rds_files[3]))
sce_myl <- readRDS(file.path(rds_subclust_dir, rds_files[4]))


# C4 ---------------------------------------------------------------------------
m3_dp <- dotplotTables(sce_c4, "meta3", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m6_dp <- dotplotTables(sce_c4, "meta6", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m4_dp <- dotplotTables(sce_c4, "meta4", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)

corr_m3 <- cor(m3_dp$expr_wide)
corr_m6 <- cor(m6_dp$expr_wide)
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

pheatmap(corr_m4,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta4 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# corr between metacluster groups
corr_mat_c3_c6 <- cor(m3_dp$expr_wide, m6_dp$expr_wide)
pheatmap(corr_mat_c3_c6,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Correlation between meta3 and meta6 clusters",
         fontsize_row = 10,
         fontsize_col = 10)


# abundance boxplot
# meta3
CATALYST::plotAbundances(sce_c4, k = "meta3", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")
# meta6
CATALYST::plotAbundances(sce_c4, k = "meta6", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")

# view batch effect
sce_c4_pbmc <- filterSCE(sce_c4, tissue_type == "PBMC")

# meta3
clust_tbl_pbmc <- table(cluster_ids(sce_c4_pbmc, "meta3"), colData(sce_c4_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
png(file.path(rds_subclust_dir, paste0("meta3_", subclust_prefix, "_pbmc_sample_cluster_heatmap.png")), width = 800, height = 600)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")
dev.off()

# meta6
clust_tbl_pbmc <- table(cluster_ids(sce_c4_pbmc, "meta6"), colData(sce_c4_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
png(file.path(rds_subclust_dir, paste0("meta6_", subclust_prefix, "_pbmc_sample_cluster_heatmap.png")), width = 800, height = 600)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")
dev.off()


# view batch - PBMC + MPE
clust_tbl_pbmc <- table(cluster_ids(sce_c4, "meta3"), colData(sce_c4)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
png(file.path(rds_subclust_dir, paste0("meta3_", subclust_prefix, "_pbmc_sample_cluster_heatmap.png")), width = 800, height = 600)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")
dev.off()

clust_tbl_pbmc <- table(cluster_ids(sce_c4, "meta6"), colData(sce_c4)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
png(file.path(rds_subclust_dir, paste0("meta6_", subclust_prefix, "_pbmc_sample_cluster_heatmap.png")), width = 800, height = 600)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")
dev.off()




# C6 ---------------------------------------------------------------------------
subclust_prefix <- c("C6")

m3_dp <- dotplotTables(sce_c6, "meta3", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m4_dp <- dotplotTables(sce_c6, "meta4", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m7_dp <- dotplotTables(sce_c6, "meta7", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)


corr_m3 <- cor(m3_dp$expr_wide)
corr_m4 <- cor(m4_dp$expr_wide)
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
         main = "Correlation between meta3 and meta4 clusters",
         fontsize_row = 10,
         fontsize_col = 10)



# abundance boxplot
# meta3
CATALYST::plotAbundances(sce_c6, k = "meta3", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")



# view batch effect
sce_c6_pbmc <- filterSCE(sce_c6, tissue_type == "PBMC")

# meta3
clust_tbl_pbmc <- table(cluster_ids(sce_c6_pbmc, "meta3"), colData(sce_c6_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
png(file.path(rds_subclust_dir, paste0("meta3_", subclust_prefix, "_pbmc_sample_cluster_heatmap.png")), width = 800, height = 600)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")
dev.off()



# Myeloid ----------------------------------------------------------------------
subclust_prefix <- c("Myeloid")

m3_dp <- dotplotTables(sce_myl, "meta3", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)
m11_dp <- dotplotTables(sce_myl, "meta11", assay = "exprs", fun = "median", scale = FALSE, q = 0.01)

corr_m3 <- cor(m3_dp$expr_wide)
corr_m11 <- cor(m11_dp$expr_wide)

pheatmap(corr_m3,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta3 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

pheatmap(corr_m11,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Corr between meta11 clusters",
         fontsize_row = 10,
         fontsize_col = 10)

# abundance boxplot
# meta3
CATALYST::plotAbundances(sce_myl, k = "meta3", by = "cluster_id", group_by = "BATCH") +
  ggtitle("Abundance by Cluster & Batch")



# view batch effect
sce_myl_pbmc <- filterSCE(sce_myl, tissue_type == "PBMC")

# meta3
clust_tbl_pbmc <- table(cluster_ids(sce_myl_pbmc, "meta3"), colData(sce_myl_pbmc)$sample_id)
clust_prop_tbl_pbmc <- prop.table(clust_tbl_pbmc, 2)
png(file.path(rds_subclust_dir, paste0("meta3_", subclust_prefix, "_pbmc_sample_cluster_heatmap.png")), width = 800, height = 600)
pheatmap(clust_prop_tbl_pbmc, cluster_rows = FALSE, cluster_cols = TRUE, main = "PBMC Samples Cluster Proportions")
dev.off()
