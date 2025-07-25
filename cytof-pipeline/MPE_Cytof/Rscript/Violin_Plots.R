library(cowplot)

# LOAD DATA
DE_c5 <- read_table(r'(D:\CHRISTOPHER_BUI\MPE_CYTOF_RDS\Cytokines\Results_Subcluster\Diff_Expr\Cytokines_C5_meta5_DE_table.txt)')

# ------------------------------------------------------------------------------
# 1. Filter for significance with p_adj_DS < 0.05, in cluster 3
# 2. Get associated markers
# 3. Violin plot
#       - For each marker
#           - pbmc vs mpe
#           - pmbc vs mpe, facet by batch
# ------------------------------------------------------------------------------

# Prepare SCE
sce_main <- readRDS(r'(D:\CHRISTOPHER_BUI\MPE_CYTOF_RDS\Cytokines\sce_main_with_meta6_full_lin_clusters.rds)')
sce_c5 <- readRDS(r'(D:\CHRISTOPHER_BUI\MPE_CYTOF_RDS\Cytokines\Results_Subcluster\RDS_Subcluster\C5_sce_subclust_full_lineage.rds)')
sce_c5$meta5 <- cluster_ids(sce_c5, "meta5")

# need sce with all markers
sce_main_c5 <- filterSCE(sce_main, meta6_full_lin == 5)
sce_main_c5$meta5 <- sce_c5$meta5

# set pbmc as reference
colData(sce_main_c5)$tissue_type <- factor(colData(sce_main_c5)$tissue_type, levels = c("PBMC", "MPE"))

sce_tmp <- filterSCE(sce_main_c5, meta5 == 3)

# set batch levels in desired order
# c("1_PBMC", "1_MPE", "2_PBMC", "2_MPE", etc.)
batch_levels <- unlist(lapply(levels(sce_tmp$BATCH), function(b) c(paste0(b, "_PBMC"), paste0(b, "_MPE"))))

sce_tmp$BATCH_tissue <- factor(
  paste0(sce_tmp$BATCH, "_", sce_tmp$tissue_type),
  levels = batch_levels
)

# ------------------------------------------------------------------------------

# select markers
sel_panel <- "Cytokines"

p_thres <- 0.05

sel_DE_c5 <- DE_c5 %>%
  dplyr::filter(cluster == 3, p_adj_DS < p_thres)

sel_markers <- sel_DE_c5$marker

# ------------------------------------------------------------------------------

# tissue type facet by markers
violin_PLOT <- scater::plotExpression(
  sce_tmp,
  features = sel_markers,
  x = "tissue_type",
  colour_by = "tissue_type",
  exprs_values = "exprs",
  scales = "free",
  one_facet = TRUE,
  scattermore = TRUE,
  point_size = 1
) +
  ggtitle(paste0(sel_panel, " - Marker Exprs by Tissue Type; p < ", p_thres)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

output_dir <- r'(D:\CHRISTOPHER_BUI\MPE_CYTOF_RDS\Cytokines)'
file_name <- file.path(output_dir, "TEST_Cytokines_C5_C3_sig_markers_by_tissue.png")
ggsave(filename = file_name, plot = violin_PLOT, width = 10, height = 12, dpi = 300)

# ------------------------------------------------------------------------------

# tissue type facet by batch, for each marker
marker_violin_PLOTS <- lapply(sel_markers, function(m) {
  scater::plotExpression(
    sce_tmp,
    features = m,
    x = "BATCH_tissue",
    colour_by = "tissue_type",
    exprs_values = "exprs",
    scales = "free",
    one_facet = TRUE,
    scattermore = TRUE,
    point_size = 1
  ) +
    scale_x_discrete(labels = function(x) gsub("_.*", "", x)) +
    ggtitle(paste0(m, " - Cytokines C5-C3; p < 0.05")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
})

final_violin_plot <- cowplot::plot_grid(plotlist = marker_violin_PLOTS, ncol = 1, labels = "AUTO")

# ADJUST BELOW FOR BETTER TITLE; USE ABOVE FOR NOW
# final_plot <- plot_grid(
#   ggdraw() +
#     draw_label(paste0(sel_panel, " - Marker Exprs by Batch & Tissue Type; p < ", p_thres),
#                hjust = 0.5,
#                size = 18),
#   big_violin_plot,
#   ncol = 1,
#   rel_heights = c(0.05, 1)
# )

output_dir <- r'(D:\CHRISTOPHER_BUI\MPE_CYTOF_RDS\Cytokines)'
file_name <- file.path(output_dir, "Cytokines_C5_C3_sig_markers_by_tissue_batch.png")
ggsave(filename = file_name, plot = final_violin_plot, width = 12, height = 4 * length(sel_markers), dpi = 300)
