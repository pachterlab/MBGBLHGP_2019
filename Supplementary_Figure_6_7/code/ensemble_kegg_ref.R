library(biomaRt)
library(tidyverse)
library(gage)
library(AnnotationDbi)
library(AnnotationHub)
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
  fields <- c("external_gene_name", "ensembl_gene_id", "description", "entrezgene")
  if (s == "Arabidopsis thaliana") {
    fields[4] <- "entrezgene_id"
  }
  gns <- getBM(fields, mart = mart)
  names(gns)[1:2] <- c("gene_name", "ensembl")
  names(gns)[4] <- "entrezgene"
  write_csv(gns, paste0("./data/Ensembl/gns_", str_replace(s, " ", "_"), ".csv"))
  #paralog <- getBM(c("ensembl_gene_id", para_ds), mart = mart)
  #write_csv(paralog, paste0("./data/Ensembl/paralog_", str_replace(s, " ", "_"), ".csv"))
}

# Get GO annotations as a data frame
pkgs <- c(paste("org", c("Mm", "Hs", "Rn", "Dr", "Dm", "Ce"), "eg.db", sep = "."),
          "org.At.tair.db") # No macaque since Bioconductor doesn't have the data
species2 <- species[-6]
for (i in seq_along(pkgs)) {
  p <- pkgs[i]
  library(p, character.only = TRUE)
  db <- get(p)
  if (species2[i] == "Arabidopsis thaliana") {
    kt <- "TAIR"
  } else {
    kt <- "ENSEMBL"
  }
  go_df <- select(db, keys = keys(db, keytype = kt), keytype = kt, 
                  columns = c("ENTREZID", "GO", "ONTOLOGY", "EVIDENCE"))
  names(go_df)[1:3] <- c("ensembl", "entrezgene", "GO.ID")
  saveRDS(go_df, paste0("./data/go/go_", str_replace(species2[i], " ", "_"), ".rds"))
}

for (s in species) {
  # Kegg
  kegg <- kegg.gsets(s, id.type = "entrez")
  # A data frame version of the kegg pathways
  kegg_df <- tibble(gene_set = names(kegg$kg.sets),
                    entrezgene = unname(kegg$kg.sets)) %>% 
    unnest()
  saveRDS(kegg, paste0("./data/kegg/kegg_", str_replace(s, " ", "_"), ".rds"))
  write_csv(kegg_df, paste0("./data/kegg/kegg_df_", str_replace(s, " ", "_"), ".csv"))
}
