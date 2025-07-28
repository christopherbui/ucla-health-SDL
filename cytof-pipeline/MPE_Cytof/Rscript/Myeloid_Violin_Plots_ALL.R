library(cowplot)

# DE tables for each TIER 2 cluster
DE_files <- list.files(r'(D:\CHRISTOPHER_BUI\MPE_CYTOF_RDS\Myeloid\Results_Subcluster\Diff_Expr)',
                       pattern = "\\.txt",
                       full.names = TRUE)

# RDS files for each TIER 2 cluster
rds_files <- list.files(r'(D:\CHRISTOPHER_BUI\MPE_CYTOF_RDS\Myeloid\Results_Subcluster\RDS_Subcluster)',
                        pattern = "\\.rds",
                        full.names = TRUE)

# metak for each TIER 2 cluster (used for DE)
metak <- list(
  C1 = "meta9",
  C2C3C4 = "meta6"
)

# grouping of TIER 1 clusters
clust_group <- list(
  C1 = c(1),
  C2C3C4 = c(2,3,4)
)
output_dir <- r'(D:\CHRISTOPHER_BUI\MPE_CYTOF_RDS\Myeloid\Results_Subcluster\TEST_VIOLIN)'


# ------------------------------------------------------------------------------
# process SCE
# sce_main <- readRDS(r'(D:\CHRISTOPHER_BUI\MPE_CYTOF_RDS\Myeloid\sce_main.rds)')
# set pbmc as reference
# colData(sce_main)$tissue_type <- factor(colData(sce_main)$tissue_type, levels = c("PBMC", "MPE"))
# # add column denoting PBMC, MPE, or Ref for plotting
# colData(sce_main)$tissue_type_mod <- case_when(
#   sce_main$tissue_type == "MPE" ~ "MPE",
#   sce_main$tissue_type == "PBMC" & sce_main$patient_id != "Ref" ~ "PBMC",
#   sce_main$tissue_type == "PBMC" & sce_main$patient_id == "Ref" ~ "Ref",
# )
# colData(sce_main)$tissue_type_mod <- factor(colData(sce_main)$tissue_type_mod, levels = c("PBMC", "MPE", "Ref"))

# load processed sce
sce_main <- readRDS(r'(D:\CHRISTOPHER_BUI\MPE_CYTOF_RDS\Myeloid\sce_main_NEW.rds)')

# set batch levels in desired order
# c("1_PBMC", "1_MPE", "1_Ref", 2_PBMC", "2_MPE", "2_Ref", etc.)
batch_levels <- unlist(lapply(levels(sce_main$BATCH), function(b) c(paste0(b, "_PBMC"), paste0(b, "_MPE"), paste0(b, "_Ref"))))

# ------------------------------------------------------------------------------
q_thres <- 0.05

# for each TIER 2 cluster group
for (i in seq_along(DE_files)) {
  # read in table
  DE <- read_tsv(DE_files[i])
  
  # filter for differentially expressed markers for each TIER 2 cluster
  sig_DE <- DE %>%
    dplyr::filter(p_adj_DS < q_thres)
  
  # if no significant markers
  if (nrow(sig_DE) == 0) {
    message("No significant markers: ", basename(DE_files[i]), "skipping...\n")
    next
  }
  
  # list of significant markers for every TIER 2 cluster
  sig_marker_list <- split(sig_DE$marker, sig_DE$cluster)
  
  message("STARTING: ", basename(DE_files[i]))
  message("Processing sce_main & sce_sub...\n")
  
  # read TIER 2 rds
  sce_sub <- readRDS(rds_files[i])
  
  tier2_name <- gsub("_sce_subclust.*", "", basename(rds_files[i])) # C1, C2C3C4
  sce_sub[[metak[[tier2_name]]]] <- cluster_ids(sce_sub, metak[[tier2_name]])
  
  sce_main_sub <- filterSCE(sce_main, meta6_full_lin %in% clust_group[[tier2_name]])
  sce_main_sub[[metak[[tier2_name]]]] <- cluster_ids(sce_sub, metak[[tier2_name]])
  
  
  # for each TIER 2 cluster that has significant marker, do violin plot
  for (c in names(sig_marker_list)) {
    sel_markers <- sig_marker_list[[c]]
    message("DOING CLUSTER: ", c)
    message(length(sel_markers), " significant markers")
    message("Preparing sce_tmp_ref & sce_tmp...\n")
    
    # contains reference
    sce_tmp_ref <- filterSCE(sce_main_sub, sce_main_sub[[metak[[tier2_name]]]] == as.numeric(c))
    # no reference
    sce_tmp <- filterSCE(sce_tmp_ref, patient_id != "Ref")
    
    sce_tmp_ref$BATCH_tissue <- factor(
      paste0(sce_tmp_ref$BATCH, "_", sce_tmp_ref$tissue_type_mod),
      levels = batch_levels
    )
    sce_tmp$BATCH_tissue <- factor(
      paste0(sce_tmp$BATCH, "_", sce_tmp$tissue_type_mod),
      levels = batch_levels
    )
    # sce_tmp has no ref cells
    sce_tmp$BATCH_tissue <- droplevels(sce_tmp$BATCH_tissue)
    
    
    # PBMC vs MPE --------------------------------------------------------------
    message("PLOTTING: PBMC vs MPE - ", tier2_name, "-C", c)
    
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
      ggtitle(paste0("Myeloid", " - Marker Exprs by Tissue Type; q < ", q_thres)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 0, hjust = 1))
    
    plot_height <- 1.6 * length(sel_markers)
    
    file_name <- file.path(output_dir, paste0("Myeloid_", tier2_name, "_C", c, "_tissue.png"))
    ggsave(filename = file_name, plot = violin_PLOT, width = 10, height = plot_height, limitsize = FALSE, dpi = 300)
    
    
    # PBMC vs MPE w/BATCH ------------------------------------------------------
    message("PLOTTING: PBMC vs MPE w/BATCH - ", tier2_name, "-C", c)
    
    violin_PLOTS_BATCH <- lapply(sel_markers, function(m) {
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
        ggtitle(paste0(m, " - Myeloid ", tier2_name, "-C", c, "; q < ", q_thres)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 0, hjust = 1))
    })
    final_violin_PLOT <- cowplot::plot_grid(plotlist = violin_PLOTS_BATCH, ncol = 1, labels = "AUTO")
    
    file_name <- file.path(output_dir, paste0("Myeloid", tier2_name, "_C", c, "_tissue_batch.png"))
    ggsave(filename = file_name, plot = final_violin_PLOT, width = 12, height = 4 * length(sel_markers), limitsize = FALSE, dpi = 300)
    
    
    # PBMC vs MPE vs Ref -------------------------------------------------------
    message("PLOTTING: PBMC vs MPE vs Ref - ", tier2_name, "-C", c)
    
    violin_PLOT <- scater::plotExpression(
      sce_tmp_ref,
      features = sel_markers,
      x = "tissue_type_mod",
      colour_by = "tissue_type_mod",
      exprs_values = "exprs",
      scales = "free",
      one_facet = TRUE,
      scattermore = TRUE,
      point_size = 1
    ) +
      ggtitle(paste0("Myeloid", " - Marker Exprs by Tissue Type; q < ", q_thres)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 0, hjust = 1))
    
    file_name <- file.path(output_dir, paste0("Myeloid_", tier2_name, "_C", c, "_tissue_ref.png"))
    ggsave(filename = file_name, plot = violin_PLOT, width = 10, height = plot_height, limitsize = FALSE, dpi = 300)
    
    
    # PBMC vs MPE vs Ref w/BATCH -----------------------------------------------
    message("PLOTTING: PBMC vs MPE vs Ref w/BATCH - ", tier2_name, "-C", c)
    
    violin_PLOTS_BATCH <- lapply(sel_markers, function(m) {
      scater::plotExpression(
        sce_tmp_ref,
        features = m,
        x = "BATCH_tissue",
        colour_by = "tissue_type_mod",
        exprs_values = "exprs",
        scales = "free",
        one_facet = TRUE,
        scattermore = TRUE,
        point_size = 1
      ) +
        scale_x_discrete(labels = function(x) gsub("_.*", "", x)) +
        ggtitle(paste0(m, " - Myeloid ", tier2_name, "-C", c, "; q < ", q_thres)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 0, hjust = 1))
    })
    final_violin_PLOT <- cowplot::plot_grid(plotlist = violin_PLOTS_BATCH, ncol = 1, labels = "AUTO")
    
    file_name <- file.path(output_dir, paste0("Myeloid_", tier2_name, "_C", c, "_tissue_batch_ref.png"))
    ggsave(filename = file_name, plot = final_violin_PLOT, width = 12, height = 4 * length(sel_markers), limitsize = FALSE, dpi = 300)
  }
}
