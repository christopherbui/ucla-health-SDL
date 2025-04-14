source("functions.R")

library(CATALYST)
library(flowCore)
library(ggplot2)

dir <- file.path(getwd(), "data")
raw_data_dir <- file.path(dir, "RawFiles")
signal_clean_dir <- file.path(data_dir, "SignalClean")

files <- list.files(signal_clean_dir,
                    pattern = ".fcs$",
                    full.names = TRUE)

# files

ff1_path <- file.path(raw_data_dir, "cat_sample1.fcs")
ff2_path <- file.path(raw_data_dir, "cat_sample2.fcs")

# read data
ff_1 <- read.FCS(ff1_path, transformation = FALSE)
ff_2 <- read.FCS(ff2_path, transformation = FALSE)

exp <- exprs(ff_2)
dim(exp)

exp_u <- unique(exp)
dim(exp_u)


# remove duplicates
for (file in files) {
  # read fcs file
  ff <- read.FCS(file, transformation = FALSE)
  
  # extract expression matrix
  expr <- exprs(ff)
  expr_dedup <- unique(expr)
}