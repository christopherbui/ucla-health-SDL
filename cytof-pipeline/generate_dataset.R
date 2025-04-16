library(CATALYST)
library(flowCore)

# set output directory
output_dir <- file.path(getwd(), "generated_data")
if(!dir.exists(output_dir)) dir.create(output_dir)

# get raw data
data("raw_data")

ids <- sampleNames(raw_data)

for (id in ids) {
  ff <- raw_data[[id]] # subset flowFrame
  fn <- sprintf("cat_sample_%s.fcs", id) # construct file name
  
  print(fn)
  
  fn <- file.path(output_dir, fn) # construct output path
  
  write.FCS(ff, fn) # write to fcs
}
