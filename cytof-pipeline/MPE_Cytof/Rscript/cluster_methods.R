# ------------------------------------------------------------------------------
# Import packages
# ------------------------------------------------------------------------------
library(CATALYST)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scran)


# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------
# set working directory
workFolder <- paste0("U:/cdbui/", "MPE_Cytof")
setwd(workFolder)

# get sample info
fin_info <- file.path("Ranalysis","mpe_cytof_sampleInfo_022625.txt")
allSampleInfo <- read.delim(fin_info, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
allSampleInfo <- allSampleInfo[, -1]

# select panel
selPanel <- c("TBNK")
# selPanel <- c("Myeloid")
# selPanel <- c("Cytokines")

# get panel info
fin_panel <- paste(selPanel, "_markers_022625.txt", sep="")
panel_info <- read.delim(file.path("Ranalysis", fin_panel), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# get channel info for fsom
lineage_idx <- which(panel_info$marker_class == "lineage")
lineage_markers <- panel_info$fcs_desc[lineage_idx]
function_idx <- which(panel_info$marker_class == "function")
function_markers <- panel_info$fcs_desc[function_idx]

# set panel input for CATALYST, requires 3 columns
tmp_panel <- panel_info[, c("fcs_colname", "antigen", "marker_class")]
tmpi <- which(grepl("function", panel_info$marker_class) == TRUE)
tmp_panel$marker_class[tmpi] <- c("state")
tmpi <- which(panel_info$marker_class == "lineage")
tmp_panel$marker_class[tmpi] <- c("type")


# set input directory
norm_dir <- file.path(workFolder, "CYTOF_data", "CytoNormed", selPanel, "updateFileName")

# set parent output directory
analysis_dir <- file.path(workFolder, "CYTOF_data", "Analysis", selPanel)
if (!dir.exists(analysis_dir)) dir.create(analysis_dir, recursive = TRUE)

# set results directory
res_dir <- file.path(analysis_dir, "Results")
if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)

# get fcs files
fcs_files <- list.files(norm_dir,
                        pattern = ".fcs$",
                        full.names = TRUE)

# get fcs corename
tmp_core <- sapply(basename(fcs_files), getCoreID, simplify = TRUE)

# match filenames and metadata
# prepare metadata for CATALYST
sampleInfo <- dplyr::filter(allSampleInfo, panel_id == selPanel)
md_info <- left_join(data.frame(file = basename(fcs_files), sample_id = tmp_core),
                     sampleInfo[, -1],
                     by = join_by(sample_id == corename))

md_info$BATCH <- factor(md_info$BATCH)



# ------------------------------------------------------------------------------
# Load SCE
# ------------------------------------------------------------------------------

# full lineage sce
file_path <- "U:/cdbui/MPE_Cytof/CYTOF_data/Analysis/TBNK/RDS/After_FSOM/sce_full_lineage_AFTER_FSOM.rds"
sce_tmp <- readRDS(file_path)




# ------------------------------------------------------------------------------
# runPCA
# NOTE:
#   - runPCA on sce if it does not have reduced dimension data
# ------------------------------------------------------------------------------

pc_prop <- 0.75 #*** proportion of total features to use
num_pcs <- 10

sce_tmp <- runPCA(sce_tmp, exprs_values = "exprs", ncomponents = num_pcs)


  

# ------------------------------------------------------------------------------
# Graph Cluster
# ------------------------------------------------------------------------------  

nn.clusters <- clusterCells(sce_tmp, use.dimred = "PCA")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  