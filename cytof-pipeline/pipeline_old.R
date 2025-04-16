#source("installation.R")
source("functions.R")

library(CATALYST)
library(flowCore)
library(ggplot2)


# create data folder where all analysis will be stored
if(!dir.exists("data")) dir.create("data")

# set data dir as main directory
dir <- file.path(getwd(), "data")


# ------------------------------------------------------------------------------
# Bead normalization -----------------------------------------------------------
#-------------------------------------------------------------------------------

# set input directory to raw fcs files (pathway to files that are to be normalized)
raw_data_dir <- file.path(dir, "RawFiles")

# set & create directory for normalized fcs files
bead_norm_dir <- file.path(dir, "BeadNorm")
if(!dir.exists(bead_norm_dir)) dir.create(bead_norm_dir)

# define full pathway to files that are to be normalized
files <- list.files(file.path(raw_data_dir),
                    pattern = ".fcs$",
                    full.names = TRUE)
files





#-------------------------------------------------------------------------------

ff <- read.FCS(file.path(raw_data_dir, "PBMC8_30min_patient4_BCR-XL_mod.fcs"))
chs <- colnames(ff)
beads <- c(139, 141, 146, 159, 165, 169, 175)
#CATALYST:::.get_bead_cols(chs, beads)
chs


# fix DNA Channel names for normalization (using CATALYST sample fcs file)
ff <- read.FCS(file.path(raw_data_dir, "PBMC8_30min_patient4_BCR-XL.fcs"), transformation = FALSE)
colnames(ff)[grep("Ir191", colnames(ff))] <- "Ir191Di"
colnames(ff)[grep("Ir193", colnames(ff))] <- "Ir193Di"
write.FCS(ff, filename = file.path(raw_data_dir, "PBMC8_30min_patient4_BCR-XL_mod.fcs"))


rename_bead_channels <- function(ff, bead_masses = c(139, 141, 146, 159, 165, 169, 175)) {
  chs <- colnames(ff)
  for (mass in bead_masses) {
    chs <- gsub(sprintf(".*?(%s).*", mass), sprintf("%sDi", mass), chs)
  }
  colnames(ff) <- chs
  return(ff)
}
# 
ff <- rename_bead_channels(ff)
grep("Ir191|Ir193", colnames(ff), value = TRUE)  # Should return both DNA channels
write.FCS(ff, filename = file.path(raw_data_dir, "PBMC8_30min_patient4_BCR-XL_mod.fcs"))
pData(parameters(ff))

expr <- exprs(ff)
dna_chs <- grep("Ir19[13]", colnames(ff), value = TRUE)
bead_chs <- grep("^(139|141|146|159|165|169|175)Di", colnames(ff), value = TRUE)

print(dna_chs)
print(bead_chs)

df <- data.frame(
  DNA = expr[, dna_chs[1]],
  Bead = expr[, bead_chs[1]]
)

ggplot(df, aes(x = DNA, y = Bead)) +
  geom_point(alpha = 0.3, color = "blue") +
  labs(title = "Bead vs DNA intensity",
       x = dna_chs[1], y = bead_chs[1]) +
  theme_minimal()





ff_x <- read.FCS(file.path(raw_data_dir, "fcs_x.fcs"))
pData(parameters(ff_x))
colnames(ff_x)[grep("Ir191", colnames(ff))] <- "Ir191Di"
colnames(ff_x)[grep("Ir193", colnames(ff))] <- "Ir193Di"
write.FCS(ff_x, file.path(raw_data_dir, "fcs_x_mod.fcs"))

#-------------------------------------------------------------------------------




# create baseline file to which all the files will be normalized
ref_sample <- baseline_file(fcs_files = files,
                            # beads = c(139, 141, 146, 159, 165, 169, 175),
                            beads = "dvs",
                            out_dir = bead_norm_dir)

# normalize file by file in the loop, saving new file with each loop execution
for (file in files) {
  # read fcs file
  ff <- flowCore::read.FCS(file,
                           transformation = FALSE,
                           truncate_max_range = FALSE)
  print("A")
  
  # bead normalize the files
  ff_norm <- bead_normalize(flow_frame = ff,
                            # beads = c(139, 141, 146, 159, 165, 169, 175),
                            beads = "dvs",
                            norm_to_ref = ref_sample,
                            out_dir = bead_norm_dir,
                            to_plot = TRUE,
                            k = 80,
                            markers_to_keep = c("CD", "HLA", "IgD", "TCR", "Ir", 
                                                "Viability", "IL", "IFNa",
                                                "TNF", "TGF", "MIP", "MCP", "Granz"))
  print("B")
  
  # save normalized fcs files
  flowCore::write.FCS(ff_norm,
                      filename = file.path(bead_norm_dir,
                                           gsub(".fcs", "_beadNorm.fcs", basename(file), ignore.case = TRUE)))
  
  print("C")
  
}





























