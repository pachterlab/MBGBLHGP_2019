# Make the gsea bar plot for all datasets
source("./code/gsea_bar_full.R")

dataset_meta <- tribble(~name_regex,	~ species, ~mapping,
                        "EMTAB7320_v2",	"Mus musculus", "org.Mm.eg.db",
                        "heart_?1k_v2",	"Mus musculus", "org.Mm.eg.db",
                        "heart_?1k_v3",	"Mus musculus", "org.Mm.eg.db",
                        "hgmm_?10k_v3",	c("Homo sapiens", "Mus musculus"), c("org.Hs.eg.db", "org.Mm.eg.db"),
                        "hgmm_?1k_v2",	c("Homo sapiens", "Mus musculus"), c("org.Hs.eg.db", "org.Mm.eg.db"),
                        "hgmm_?1k_v3",	c("Homo sapiens", "Mus musculus"), c("org.Hs.eg.db", "org.Mm.eg.db"),
                        "neuron_?10k_v3",	"Mus musculus", "org.Mm.eg.db",
                        "pbmc_?10k_v3",	"Homo sapiens", "org.Hs.eg.db",
                        "pbmc_?1k_v3",	"Homo sapiens", "org.Hs.eg.db",
                        "SRR6956073_v2", "Danio rerio", "org.Dr.eg.db",
                        "SRR6998058_v2", "Mus musculus", "org.Mm.eg.db",
                        "SRR7299563_v2", "Rattus norvegicus", "org.Rn.eg.db",
                        #"SRR7692543_v2", "Macaca fascicularis",
                        "SRR8206317_v2", "Mus musculus", "org.Mm.eg.db",
                        "SRR8257100_v2", "Arabidopsis thaliana", "org.At.tair.db",
                        "SRR8327928_v2", "Homo sapiens", "org.Hs.eg.db",
                        "SRR8513910_v2", "Drosophila melanogaster", "org.Dm.eg.db",
                        "SRR8524760_v2", "Homo sapiens", "org.Hs.eg.db",
                        "SRR8599150_v2", "Mus musculus", "org.Mm.eg.db",
                        "SRR8611943_v2", "Caenorhabditis elegans", "org.Ce.eg.db",
                        "SRR8639063_v2", "Mus musculus", "org.Mm.eg.db")
dataset_meta <- dataset_meta %>% 
  mutate(name_use = str_remove(name_regex, "_\\?"))

# Now begins the run
dataset_meta <- dataset_meta %>% 
  mutate(plot = pmap(list(species, name_regex, mapping, name_use),
                     ~ gsea_bar_full(..1, ..2, ..3, ..4)))
