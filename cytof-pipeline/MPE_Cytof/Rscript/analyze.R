# load raw & processed dotplot matrices
raw_dp <- read.delim(file.path(res_dir, "Dotplot/PC_6/meta20_dotplot_mean_expression_matrix_wide.txt"), header = TRUE, check.names = FALSE)
proc_dp <- read.delim(file.path(res_dir, "Dotplot/PC_6/meta20_processed_dotplot_mean_expression_matrix_wide.txt"), header = TRUE, check.names = FALSE)

# order by markers for consistency
common_markers <- intersect(raw_dp$Symbol, proc_dp$Symbol)
raw_dp <- raw_dp[order(raw_dp$Symbol), ]
proc_dp <- proc_dp[order(proc_dp$Symbol), ]

# only get numeric columns for correlation calculation
raw_dp <- raw_dp[, !(names(raw_dp) %in% "Symbol")][, -1]
rownames(raw_dp) <- NULL

proc_dp <- proc_dp[, !(names(proc_dp) %in% "Symbol")][, -1]
rownames(proc_dp) <- NULL

c1_raw <- raw_dp[, 1]
c5_proc <- proc_dp[, 5]

cor(c1_raw, c5_proc)
