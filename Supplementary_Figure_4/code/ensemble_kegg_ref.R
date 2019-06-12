library(biomaRt)
library(tidyverse)
library(gage)

species <- c("Mus musculus", "Homo sapiens", "Rattus norvegicus", 
             "Danio rerio", "Drosophila melanogaster",
             "Macaca fascicularis", "Caenorhabditis elegans", "Arabidopsis thaliana")

for (s in species) {
  # Ensembl
  if (s != "Arabidopsis thaliana") {
    ds <- BUSpaRse:::species2dataset(s)
    v <- 94
    host <- "www.ensembl.org"
    m <- "ensembl"
    para_ds <- str_replace(ds, "gene_ensembl", "paralog_ensembl_gene")
  } else {
    ds <- "athaliana_eg_gene"
    v <- NULL
    host <- "plants.ensembl.org"
    m <- "plants_mart"
    para_ds <- str_replace(ds, "eg_gene", "paralog_ensembl_gene")
  }
  mart <- useEnsembl(m, ds, version = v, host = host)
  gns <- getBM(c("external_gene_name", "ensembl_gene_id", "description", "entrezgene"), 
               mart = mart)
  names(gns)[1:2] <- c("gene_name", "gene")
  write_csv(gns, paste0("./data/Ensembl/gns_", str_replace(s, " ", "_"), ".csv"))
  paralog <- getBM(c("ensembl_gene_id", para_ds), mart = mart)
  write_csv(paralog, paste0("./data/Ensembl/paralog_", str_replace(s, " ", "_"), ".csv"))
   Kegg
  kegg <- kegg.gsets(s, id.type = "entrez")
  # A data frame version of the kegg pathways
  kegg_df <- tibble(gene_set = names(kegg$kg.sets),
                    entrezgene = unname(kegg$kg.sets)) %>% 
    unnest()
  saveRDS(kegg, paste0("./data/kegg/kegg_", str_replace(s, " ", "_"), ".rds"))
  write_csv(kegg_df, paste0("./data/kegg/kegg_df_", str_replace(s, " ", "_"), ".csv"))
}
