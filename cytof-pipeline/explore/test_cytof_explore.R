library(flowCore)

# fixes DNA channel naming
standardize_channel_names <- function(ff) {
  chs <- colnames(ff)
  # Extract mass from parentheses, e.g., Ir191 => Ir191Di
  for (i in seq_along(chs)) {
    matches <- regmatches(chs[i], regexec("\\(([^)]+)\\)", chs[i]))
    if (length(matches[[1]]) > 1) {
      metal <- matches[[1]][2]
      chs[i] <- paste0(metal, "Di")
    }
  }
  colnames(ff) <- chs
  return(ff)
}

# load fcs sample
fcs_path1 = "/Users/cdbui/Documents/cytof-pipeline/Catalyst/PBMC8_fcs_files/PBMC8_30min_patient1_BCR-XL.fcs"
fcs_path2 = "/Users/cdbui/Documents/cytof-pipeline/Catalyst/PBMC8_fcs_files/PBMC8_30min_patient2_BCR-XL.fcs"

ff_1 <- read.FCS(fcs_path1, transformation = FALSE)
ff_2 <- read.FCS(fcs_path2, transformation = FALSE)

ff_fix1 <- standardize_channel_names(ff_1)
ff_fix2 <- standardize_channel_names(ff_2)


colnames(ff_1)
colnames(ff_fix1)


sce <- prepData(ff_1)
int_colData(sce)

colnames(sce)
rowData(sce)

exprs_data <- exprs(ff_fix1)
colnames(exprs_data)

channels(sce)
assayNames(sce)

# write fcs files
# write.FCS(ff_fix1,
#           filename = file.path(getwd(), "ff_fix1.fcs"))
# write.FCS(ff_fix2,
#           filename = file.path(getwd(), "ff_fix2.fcs"))

# -------------------------------------
# COMPARE
# -------------------------------------
data("raw_data")
ff_c <- raw_data[[1]]
pData(parameters(ff_c))
colnames(ff_c)

sce <- prepData(ff_c)
ff_norm <- normCytof(sce,
                     beads = "dvs",
                     # beads = c(139, 141, 146, 159, 165, 169, 175),
                     k = 50,
                     assays = c("counts", "exprs"),
                     overwrite = FALSE,
                     remove_beads = TRUE)
t <- sce2fcs(ff_norm$data)
colnames(t)







# check norm files
fpath <- "C:/Users/cdbui/Desktop/cat_data/BeadNorm/cat_sample1_norm.fcs"
ff <- read.FCS(fpath, transformation = FALSE)
colnames(ff)
ff$data
