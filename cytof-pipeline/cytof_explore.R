library(flowCore)


data_dir <- file.path(path.expand("~"), "cytof-pipeline/Catalyst/PBMC8_fcs_files")

fcs_file_name <- "PBMC8_30min_patient1_BCR-XL.fcs"
fcs_file_path <- file.path(data_dir, fcs_file_name)

fcs_ref <- "PBMC8_30min_patient1_Reference.fcs"
fcs_ref_path <- file.path(data_dir, fcs_ref)

# load sample fcs
fcs <- read.FCS(fcs_file_path, transformation = FALSE)

# load sample fcs ref
fcs_ref <- read.FCS(fcs_ref_path, transformation = FALSE)

param_info <- pData(parameters(fcs))
param_info

sce <- prepData(fcs)
rownames(sce)
colnames(sce)
colnames(fcs)
pData(parameters(fcs))

