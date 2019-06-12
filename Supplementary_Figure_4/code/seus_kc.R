# Get half of seus and use my computer to do Leiden clustering
library(Seurat)
library(purrr)
seus <- readRDS("./output/neuron10k_v3/seurat.rds")
seus_use <- seus[c("kallisto", "cellranger")]
seus_use <- map(seus_use, SetAssayData, slot = "scale.data", new.data = matrix(nrow = 0, ncol = 0))
saveRDS(seus_use, file = "./output/neuron10k_v3/seus_kc.rds")
# Then I'll modify the seus_kc.rds on my computer for Leiden clustering
# The code below was run on my Mac
library(reticulate)
use_python("/anaconda3/bin/python")
library(Seurat)
library(purrr)
seus <- readRDS("./seus_kc.rds")
seus <- map(seus, FindClusters, algorithm = 4, resolution = 1)
saveRDS(seus, "./seus_kc_leiden.rds")
