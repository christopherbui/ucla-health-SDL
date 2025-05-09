.dotplot <- function(x, k, 
                     assay = "exprs", fun = "median", 
                     scale = TRUE, q = 0.01, 
                     pal = hcl.colors(11, "viridis")) {

x <- sce_ref  
  
k <- "meta20"
assay <- "exprs"
fun <- "mean"
scale <- TRUE
q <- 0.01
pal <- hcl.colors(11, "viridis")
  
x$cluster_id <- cluster_ids(x, k)
es <- assay(x, assay)
th <- rowMedians(es)

cs <- seq_len(ncol(x))
cs <- split(cs, x$cluster_id)
fq <- sapply(cs, function(i)
  rowMeans(es[, i, drop = FALSE] > th))

# compute median expression by cluster
lab <- paste(fun, assay)
ms <- CATALYST:::.agg(x, by = "cluster_id", assay = assay, fun = fun)
if (scale) {
  lab <- paste("scaled", lab)
  ms <- CATALYST:::.scale_exprs(ms, q = q)
}

# do hierarchical clustering on rows & columns
cluster_order <- function(x) order.dendrogram(as.dendrogram(hclust(dist(x))))
ro <- colnames(ms)[cluster_order(t(ms))]
co <- rownames(ms)[cluster_order(ms)]

df <- cbind(melt(ms), fq = melt(fq)$value)

df_wide <- dcast(df, Var1 ~ Var2, value.var = "value")

# write.table(df, file = paste0("dotplot_", fun, "_expression_matrix.txt"), sep = "\t", quote = FALSE, col.names = NA)
write.table(df_wide, file = paste0("dotplot", fun, "_expression_matrix_wide.txt"), sep = "\t", quote = FALSE, col.names = NA)

ggplot(df, aes(Var1, Var2, col = value, size = fq, label = sprintf("%.2f", value))) + 
  geom_point() +
  geom_text(color = "black", size = 2.5, vjust = 0.5) +  # show mean/median value
  scale_x_discrete("marker", limits = co, expand = c(0, 0.5)) +
  scale_y_discrete("cluster_id", limits = ro, expand = c(0, 0.5)) +
  scale_color_gradientn(lab, breaks = seq(0, 1, 0.5), colors = pal) +
  scale_size_continuous(range = c(0, 5), 
                        labels = formatC(seq(0, 1, 0.25), 2, format = "f")) +
  guides(
    color = guide_colorbar(order = 1),
    size = guide_legend("% cells with expr. above\n global marker median")) +
  coord_equal() + theme_linedraw() + theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1))
}