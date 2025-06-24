# Run in R v4.3.3
library(stringr)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(pals)
library(scales)
library(ggplot2)
library(pheatmap)
library(cowplot)
library(reshape2)
library(ComplexHeatmap)

library(flowCore)
library(flowAI)
library(flowCut)
library(cytutils)   #clauclting aof/change fcs name

library(CATALYST)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(diffcyt)
library(scater)
library(FlowSOM)
library(bluster)
library(class)    #req

dotplotTables <- function(sce,
                       cluster_name,
                       assay = "exprs",
                       fun = "median",
                       scale = TRUE,
                       q = 0.01){

  # this allows us to select both raw (from FSOM cluster()) & processed clusters
  isIDExist <- which(colnames(colData(sce))==cluster_name)
  if (length(isIDExist)==1){
     es_cluster_id <- colData(sce)[,cluster_name]
  } else {
     es_cluster_id <- cluster_ids(sce, cluster_name)
  }
  
  es <- assay(sce, assay) # expression matrix
  if (length(grep("Symbol",colnames(rowData(sce))))>0){
     rownames(es) <- rowData(sce)$Symbol
     marker_symbols <- rowData(sce)$Symbol
  }
  th <- rowMedians(es)  # median threshold
  
  cs <- seq_len(ncol(es))  # indices vector with length = number of cells
  cs <- split(cs, es_cluster_id)   # split indices into groups based on cid_sel
  
  # return fraction of cells > th per marker
  fq <- sapply(cs, function(i) {
      rowMeans(es[, i, drop = FALSE] > th)
  })
  
  # compute mean/median expression by cluster
  # ms is markers x clusters
  lab <- paste(fun, assay)
  if (fun =="median"){
      ms <- sapply(cs, function(i) rowMedians(es[, i, drop = FALSE]))
  }else{
      ms <- sapply(cs, function(i) rowMeans(es[, i, drop = FALSE]))
  }
  rownames(ms) <- rownames(es)   

  # if scale=TRUE
  if(scale==TRUE){
     lab <- paste("scaled", lab)
     ms <- CATALYST:::.scale_exprs(ms, q = q)    
  } 
  
  # long format for ggplot
  df <- melt(ms)            # melt(ms): Var1 = marker name, Var2 = cluster id, value = fun() value of marker across clusters
  df$fq <- melt(fq)$value   # melt(fq): Var1 = marker name, Var2 = cluster id, value = proportion > th
  
  colnames(df)[colnames(df) == "Var1"] <- "Symbol"
  colnames(df)[colnames(df) == "Var2"] <- "cluster"
  
  # need to turn cluster column into factor for ggplot
  df$cluster <- factor(df$cluster)

  list(expr_wide = ms, expr_long = df)
}

dotplotFig <- function(df,
                   pal = hcl.colors(11, "viridis"),
                   lab = c("scaled median exprs")){

  marker_symbols <- df$Symbol[!duplicated(df$Symbol)]
  max_expr <- floor(max(df$value))
  
  dot_PLOT <- ggplot(df, aes(Symbol, cluster, col = value, size = fq, label = sprintf("%.2f", value))) +
    geom_point() +
    scale_x_discrete("marker", limits = marker_symbols, expand = c(0, 1)) +
    scale_y_discrete("cluster_id", limits = unique(df$cluster), expand = c(0, 0.5)) +
    scale_color_gradientn(lab, breaks = seq(0, max_expr+0.5, max_expr/4), colors = pal) +
    scale_size_continuous(
      range = c(0, 5),
      labels = formatC(seq(0, 1, 0.25), 2, format = "f")
    ) +
    guides(
      color = guide_colorbar(order = 1),
      size = guide_legend("frac. cells with expr >\n global marker median")
    ) +
    coord_equal() +
    theme_linedraw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8, margin = margin(r = 4))
    )
  
  return(dot_PLOT)
}


workFolder <- paste("D:/","MPE_Cytof/",sep="")
setwd(workFolder)

# -------------------
inputFolder <- c("V:/cdbui/MPE_Cytof/CYTOF_data/Analysis/TBNK/RDS")

tmp_fin <- c("sce_main.rds")
fin <- file.path(inputFolder,tmp_fin)        
sce_all <- readRDS(file=fin)

tmp_fin <- c("sce_full_lineage.rds")
fin <- file.path(inputFolder,tmp_fin)        
sub_sce_1 <- readRDS(file=fin)

tmp_fin <- c("sce_full_lineage_AFTER_FSOM.rds")
fin <- file.path(inputFolder,"After_FSOM",tmp_fin)        
sub_sce_1_after <- readRDS(file=fin)


tmp_fin <- c("sce_subset_lineage.rds")
fin <- file.path(inputFolder,tmp_fin)        
sub_sce_2 <- readRDS(file=fin)


tmp_fin <- c("sce_subset_lineage_AFTER_FSOM.rds")
fin <- file.path(inputFolder,"After_FSOM",tmp_fin)        
sub_sce_2_after <- readRDS(file=fin)

rowMeans(assay(sub_sce_2,"exprs"))
rowMeans(assay(sub_sce_2_after,"exprs"))
rowMeans(assay(sce_all,"exprs"))

# check expression data
es_2 <- assay(sub_sce_2_after,assay = "exprs")
es_1 <- assay(sce_all,assay = "exprs")
idx_1v2 <- left_join(data.frame(id=rownames(es_2)),
              data.frame(id=rownames(es_1),c(1:nrow(es_1))))
es_1s <- es_1[idx_1v2[,2],]

tmp <- cor(t(es_1s),t(es_2))
pheatmap(tmp)

# redo FSOM
maxk <- 20  #***specify number of consensus clusters, default = 20
sel_markers <- rownames(sub_sce_1)
sub_sce_1_LT <- cluster(sub_sce_1,
               features = sel_markers,
               xdim = 10,
               ydim = 10,
               maxK = maxk,
               seed = 1234)

maxk <- 20  #***specify number of consensus clusters, default = 20
sel_markers <- rownames(sub_sce_2)
set.seed(100)
sub_sce_2_LT <- cluster(sub_sce_2,
               features = sel_markers,
               xdim = 10,
               ydim = 10,
               maxK = maxk)
              

## regenerate the susbet and rerun FSOM
tmp_sel_idx <- is.element(rownames(sce_all),
                     rownames(sub_sce_2))
sel_markers <- rownames(sce_all)[tmp_sel_idx]
sub_sce_2LT <- sce_all[sel_markers,]
set.seed(1234)
sub_sce_2_LT2 <- CATALYST::cluster(sub_sce_2LT,
               features = rownames(sub_sce_2LT),
               xdim = 10,
               ydim = 10,
               maxK = maxk)
               ##seed = 1234)


# check cell distributions
selMetaK <- c("meta8")   #****
tmp_org <- cluster_ids(sub_sce_2_after, selMetaK)  #****
tmp_mod <- cluster_ids(sub_sce_2_LT2, selMetaK)     #*****
tmp_1v0 <- prop.table(table(tmp_org,tmp_mod),1)
pheatmap(tmp_1v0,cluster_row=FALSE,cluster_col=FALSE)

# check expression data of clusters
selMetaK <- c("meta20")   #****
scale_option <- FALSE
sel_aggregateFxn <- c("median")
lab4dotplot <- paste(ifelse(scale_option,"scaled","non-scaled"),
                     sel_aggregateFxn,"exprs", sep=" ")
tmp_ftab_0 <- dotplotTables(sub_sce_2_after,cluster_name=selMetaK,
                   assay = "exprs",
                   fun = sel_aggregateFxn,
                   scale = scale_option,
                   q = 0.01)
tmp_fig0 <- dotplotFig(tmp_ftab_0$expr_long,lab=lab4dotplot)

tmp_ftab_1 <- dotplotTables(sub_sce_2_LT2,cluster_name=selMetaK,
                   assay = "exprs",
                   fun = sel_aggregateFxn,
                   scale = scale_option,
                   q = 0.01)
tmp_ftab_1$expr_wide
tmp_fig1 <- dotplotFig(tmp_ftab_1$expr_long,lab=lab4dotplot)



