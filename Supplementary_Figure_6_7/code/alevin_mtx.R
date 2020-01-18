# Convert all alevin output into sparse matrix
library(purrr)
###########
# Get dataset names
alevin_dirs <- list.dirs("../salmon_out/", recursive = FALSE)
name_use <- str_replace(alevin_dirs, "../salmon_out//salmon_", "") %>% 
  str_replace("_out", "")
# Not redoing the ones I've done
alevin_use <- alevin_dirs[-c(4:8,15)]
name_use_sub <- name_use[-c(4:8,15)]
for (i in seq_along(alevin_use)) {
  cat("Processing dataset", name_use_sub[i], "\n")
  readAlevin_sparse(paste0(alevin_use[i], "/alevin/quants_mat.gz"),
                    save_mtx = TRUE, 
                    save_path = paste0("./output/alevin_mtx/", name_use_sub[i]))
}
