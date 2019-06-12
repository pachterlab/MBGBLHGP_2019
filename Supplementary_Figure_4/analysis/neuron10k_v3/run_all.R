# Run all analyses for this dataset, in the correct order
library(rmarkdown)

file_ordered <- c("standard.Rmd", "cluster_markers.Rmd", "cell_type_annotation.Rmd",
                  "cluster.Rmd", "correlation.Rmd", "DE_methods.Rmd", "paralogs.Rmd",
                  "pseudotime.Rmd")

for (f in file_ordered) {
  render(f, knit_root_dir = "/home/single_cell_analysis/brain_storm")
}

# supp_figs.R and the notebook for species mixing are independent from these.
