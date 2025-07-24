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
sel_panel <- "Cytokines"
# select markers
p_thres <- 0.05

sel_DE_c5 <- DE_c5 %>%
  dplyr::filter(cluster == 3, p_adj_DS < p_thres)

sel_markers <- sel_DE_c5$marker


# tissue type facet by markers
violin_PLOT <- scater::plotExpression(
  sce_c5c3,
  features = sel_markers,
  x = "tissue_type",
  colour_by = "tissue_type",
  exprs_values = "exprs",
  scales = "free",
  one_facet = TRUE,
  scattermore = TRUE,
  point_size = 1
) +
  ggtitle(paste0(sel_panel, " - Marker Exprs By Tissue Type; p < ", p_thres)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

output_dir <- r'(D:\CHRISTOPHER_BUI\MPE_CYTOF_RDS\Cytokines)'
file_name <- file.path(output_dir, "Cytokines_C5_C3_sig_markers_by_tissue.png")
ggsave(filename = file_name, plot = violin_PLOT, width = 10, height = 12, dpi = 300)


# tissue type facet by batch, for each marker
for (m in sel_markers) {
  message("Plotting Marker: ", m)
  
  violin_PLOT <- scater::plotExpression(
    sce_c5c3,
    features = m,
    x = "tissue_type",
    colour_by = "tissue_type",
    exprs_values = "exprs",
    scales = "free",
    one_facet = FALSE,
    scattermore = TRUE,
    point_size = 1
  ) +
    facet_wrap(~BATCH, ncol = 4) +
    ggtitle(paste0(sel_panel, " - Marker Exprs by Tissue Type & Batch; p < ", p_thres)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1))
  
  file_name <- file.path(output_dir, paste0("Cytokines_C5_C3_", m, "_by_tissue_batch.png"))
  ggsave(filename = file_name, plot = violin_PLOT, width = 12, height = 8, dpi = 300)
}
