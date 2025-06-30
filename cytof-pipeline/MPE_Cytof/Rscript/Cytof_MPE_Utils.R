# ------------------------------------------------------------------------------
# UTILITY FUNCTIONS FOR MPE CYTOF
# ------------------------------------------------------------------------------


#' Get lineage markers
#'
#' @param sce SingleCellExperiment object
#'
#' @returns Vector of lineage markers
#' 
get_lineage_markers <- function(sce) {
  # lineage markers are of class 'type' in sce rowdata
  type_markers <- rownames(sce)[rowData(sce)$marker_class == "type"]
  
  return(type_markers)
}




#' Calculate Marker Expression Quantiles
#'
#' @param sce SingleCellExperiment object
#'
#' @returns Dataframe of marker quantiles
#' 
calculate_marker_quantiles <- function(sce) {
  
  # expression matrix
  es <- assay(sce, "exprs")
  
  marker_names <- rowData(sce)$Symbol
  
  marker_stats <- apply(es, 1, function(x) {
    c(
      mean = mean(x),
      `25%` = unname(quantile(x, 0.25)),
      `50%` = unname(quantile(x, 0.50)),
      `75%` = unname(quantile(x, 0.75)),
      `90%` = unname(quantile(x, 0.90)),
      `95%` = unname(quantile(x, 0.95)),
      max = max(x)
    )
  })
  
  marker_stats_df <- as.data.frame(t(marker_stats))
  
  rownames(marker_stats_df) <- marker_names
  
  return(marker_stats_df)
}




#' Get dataframe with PC and variance explained
#'
#' @param sce SingleCellExperiment object
#'
#' @returns Dataframe with PC number and variance explained
#' 
get_pca_variance_explained <- function(sce) {
  pca_matrix <- reducedDim(sce, "PCA")
  pca_variance <- apply(pca_matrix, 2, var)
  pca_variance_explained <- pca_variance / sum(pca_variance)
  
  elbow_df <- data.frame(
    PC = seq_along(pca_variance_explained),
    VarianceExplained = pca_variance_explained
  )
  
  return(elbow_df)
}




#' Plot elbow curve
#'
#' @param var_percent vector of PC variance explained, of PC percent of variance explained
#'
#' @returns Eblow plot (ggplot object)
#' 
plot_elbow_plot <- function(var_percent) {
  elbow_df <- data.frame(
    PC = seq_along(var_percent),
    VarianceExplained = var_percent
  )
  
  elbow_PLOT <- ggplot(elbow_df, aes(x = PC, y = VarianceExplained)) +
    geom_point(color = "#4287f5", size = 2) +
    geom_line(color = "#4287f5", linewidth = 0.9) +
    geom_text(aes(label = round(VarianceExplained, 2)),
              hjust = -0.5, vjust = -0.5, size = 3.5) +
    scale_x_continuous(breaks = seq_along(var_percent)) +
    labs(title = paste0("PCA Elbow Plot", " - ", selPanel),
         x = "Principal Component",
         y = "Percent of Variance Explained") +
    theme_light(base_size = 14)
  
  return(elbow_PLOT)
}



scran_analysis <- function(sce, cluster_name, test_type = "wilcox", average = "median") {
  clusters <- sort(unique(colData(sce)[[cluster_name]]))
  
  res_by_cluster <- list()
  
  for (c in clusters) {
    message("Processing Cluster: ", c)
    
    # subset by cluster
    sce_tmp <- filterSCE(sce, !!sym(cluster_name) == c)
    
    tmp_de_pv <- scran::findMarkers(sce_tmp,
                                    groups = colData(sce_tmp)$tissue_type,
                                    assay.type = "exprs",
                                    test.type = test_type,
                                    pval.type = "all",
                                    min.prop = 0.25,
                                    direction = "up")
    
    sum_out <- scran::summaryMarkerStats(sce_tmp,
                                         groups = colData(sce_tmp)$tissue_type,
                                         assay.type = "exprs",
                                         average = average)
    
    cluster_tissue_list <- list()
    
    for (tissue in c("PBMC", "MPE")) {
      if (tissue %in% names(tmp_de_pv)) {
        # Extract and tag marker names
        de_df <- as.data.frame(tmp_de_pv[[tissue]])
        de_df$marker <- rownames(de_df)
        
        sum_df <- as.data.frame(sum_out[[tissue]])
        sum_df$marker <- rownames(sum_df)
        
        # inner join on marker
        merged <- inner_join(sum_df, de_df, by = "marker")
        
        # add cluster id & tissue type info
        merged$cluster <- c
        merged$tissue_type <- tissue
        
        cluster_tissue_list[[tissue]] <- merged
      }
    }
    
    # concatenate PBMC and MPE results for this cluster
    cluster_all_tissue <- bind_rows(cluster_tissue_list)
    
    # save to master list
    res_by_cluster[[as.character(c)]] <- cluster_all_tissue
  }
  
  # concatenate all info
  all_info <- bind_rows(res_by_cluster)
  
  # reorder columns
  all_info <- all_info[, c("marker", "cluster", "tissue_type", setdiff(colnames(all_info), c("marker", "cluster", "tissue_type")))]
  
  return(all_info)
}

violin_plot <- function()











