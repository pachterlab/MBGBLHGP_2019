# Make the gsea bar plot for all datasets
source("./code/gsea_bar_full.R")

dataset_meta <- tribble(~name_regex,	~ species,
                        "EMTAB7320_v2",	"Mus musculus",
                        "heart_?1k_v2",	"Mus musculus",
                        "heart_?1k_v3",	"Mus musculus",
                        "hgmm_?10k_v3",	c("Homo sapiens", "Mus musculus"),
                        "hgmm_?1k_v2",	c("Homo sapiens", "Mus musculus"),
                        "hgmm_?1k_v3",	c("Homo sapiens", "Mus musculus"),
                        "neuron_?10k_v3",	"Mus musculus",
                        "pbmc_?10k_v3",	"Homo sapiens",
                        "pbmc_?1k_v3",	"Homo sapiens",
                        "SRR6956073_v2", "Danio rerio",
                        "SRR6998058_v2", "Mus musculus",
                        "SRR7299563_v2", "Rattus norvegicus",
                        "SRR7692543_v2", "Caenorhabditis elegans",
                        "SRR8206317_v2", "Mus musculus",
                        "SRR8257100_v2", "Arabidopsis thaliana",
                        "SRR8327928_v2", "Homo sapiens",
                        "SRR8513910_v2", "Drosophila melanogaster",
                        "SRR8524760_v2", "Homo sapiens",
                        "SRR8599150_v2", "Mus musculus",
                        "SRR8611943_v2", "Macaca fascicularis",
                        "SRR8639063_v2", "Mus musculus")
dataset_meta <- dataset_meta %>% 
  mutate(name_use = str_remove(name_regex, "_\\?"))

# Now begins the run
dataset_meta <- dataset_meta %>% 
  mutate(plot = pmap(list(species, name_regex, name_use),
                     ~ try(gsea_bar_full(..1, ..2, ..3))))
