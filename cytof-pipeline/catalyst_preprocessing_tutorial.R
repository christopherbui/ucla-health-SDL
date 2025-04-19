# load required packages
library(CATALYST)
library(cowplot)
library(flowCore)
library(ggplot2)
library(SingleCellExperiment)

# load "raw_data" dataset from CATALYST
data("raw_data")

# INSPECT ----------------------------------------
class(raw_data) # flowSet containing 2 flowFrames
sampleNames(raw_data)

# exprs() shows the expression matrix for flowFrame object
# raw_data_1.fcs with 2500 cells & 67 features
ff_1 <- raw_data[[1]]
ff_2 <- raw_data[[2]]

colnames(ff_1)
flowCore::parameters(ff_1)$desc
pData(parameters(ff_1))
keyword(ff_1)

write.FCS(ff_1, filename = file.path("cat_v1.fcs"))
write.FCS(ff_2, filename = file.path("cat_v2.fcs"))

# view column names, and shape
colnames(ff_1)
dim(ff_1)
# ------------------------------------------------


# construct SingelCellExperiment object
sce <- prepData(raw_data)
dim(sce)

# view number of events per sample/flowFrame
table(sce$sample_id)

# view non-isotope features
names(int_colData(sce))

# NORMALIZATION ----------------------------------------
sce <- prepData(raw_data)

# apply normalization; keep raw data
res <- normCytof(sce, beads = "dvs", k = 50, assays = c("counts", "exprs"), overwrite = FALSE)

# check proportion of beads & removed events
n <- ncol(sce); ns <- c(ncol(res$beads), ncol(res$removed))
data.frame(
    check.names = FALSE,
    "#" = c(ns[1], ns[2]),
    "%" = 100 * c(ns[1]/n, ns[2]/n),
    row.names = c("beads", "removed")
)

# extract data excluding beads & doublets, including normalized data
sce <- res$data
assayNames(sce)

# plot DNA intensity vs bead signal
res$scatter
# plot smoothed bead intensities
res$lines
# ------------------------------------------------



# DEBARCODING ----------------------------------------
