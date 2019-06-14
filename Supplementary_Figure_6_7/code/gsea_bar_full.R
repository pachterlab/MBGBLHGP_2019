library(BUSpaRse)
library(tidyverse)
library(DropletUtils)
library(Matrix)
library(Seurat)
code_files <- paste0("./code/", c("find_clust_mirror", "gsea", "plot_comparison",
                                  "read_count_output", "utils"), ".R")
walk(code_files, source)

#' Plot the GSEA bar plot directly from matrix
#' 
#' This function directly reads the gene count matrix for kallisto and CellRanger, 
#' filters the barcodes and genes, normalizes it with Seurat, finds genes 
#' differentially "expressed" between kallisto and CellRanger, and does GSEA on
#' "marker genes" with adjusted p-values less than 0.05. Then it makes the bar plot
#' showing the number of marker genes in each enriched 
#' 
#' @param species Latin name of the species of interest.
#' @param name_regex Regex for the name of the dataset, such as heart1k_v2 and
#' SRR6956073_v2.
#' @param name_use Name to use, not the regex.
#' @param min_cells Minimum number of cells, used for filtering barcodes.
#' @param verbose Whether to display progress.
#' @return Nothing in the R session. The GSEA method DE bar plot is saved, as
#' well as the data frame used to produce it.
gsea_bar_full <- function(species, name_regex, 
                          name_use = name_regex, 
                          min_cells = 1000, verbose = TRUE) {
  gns <- map_dfr(species, ~ read_csv(paste0("./data/Ensembl/gns_", 
                                        str_replace(.x, " ", "_"), ".csv")))
  kegg_kgsets <- map(species, ~ readRDS(paste0("./data/kegg/kegg_", 
                                        str_replace(.x, " ", "_"), 
                                        ".rds"))$kg.sets) %>% 
    purrr::reduce(c)
  kegg_df <- map_dfr(species, ~ read_csv(paste0("./data/kegg/kegg_df_", 
                                            str_replace(.x, " ", "_"), ".csv")))
  
  # Matrices----------------
  mats <- list()
  # Read in kallisto matrix
  mats$kallisto <- read_count_output(get_dir(name_regex, 
                                             "/home/single_cell_analysis/kallisto_out_single"),
                                     "genes", tcc = FALSE)
  # Filter barcodes
  bc_rank <- barcodeRanks(mats$kallisto, lower = min_cells)
  mats$kallisto <- mats$kallisto[, Matrix::colSums(mats$kallisto) > bc_rank$inflection]
  # Read in filtered CellRanger matrix
  mats$cellranger <- read_cellranger(paste0(get_dir(name_regex, "../cellranger_out"),
                                            "/outs/filtered_feature_bc_matrix"))
  colnames(mats$cellranger) <- colnames(mats$cellranger) %>% str_remove("-1")
  # Only keep genes that are detected
  mats <- map(mats, ~ .x[Matrix::rowSums(.x) > 0, ])
  
  # Only keep barcodes and genes present in both
  bcs_use <- intersect(colnames(mats$kallisto), colnames(mats$cellranger))
  genes_use <- intersect(rownames(mats$kallisto), rownames(mats$cellranger))
  mats <- map(mats, ~ .x[genes_use, bcs_use])
  
  # Determine gene ID type
  if (all(str_detect(genes_use, "(ENS[A-Z]*G\\d+)|(^AT.G\\d+)|(WBGene\\d+)|(FBgn\\d+)"))) {
    col_use <- "gene"
  } else if (all(str_detect(genes_use, "^\\d+$"))) {
    col_use <- "entrezgene"
  } else {
    col_use <- "gene_name"
  }
  
  # clean Ensembl gene symbols for mixed species datasets
  if (col_use == "gene" && !str_detect(genes_use, "^ENS[A-Z]*G\\d+$")) {
    genes_use <- str_remove(genes_use, "^hg19_")
    genes_use <- str_remove(genes_use, "^mm10_")
    for (i in seq_along(mats)) {
      rownames(mats[[i]]) <- genes_use
    }
  }
  # The universe for GSEA
  entrez_inter <- gns$entrezgene[gns[[col_use]] %in% genes_use]
  entrez_inter <- entrez_inter[!is.na(entrez_inter)]
  # Seurat-----------------
  seus <- imap(mats,
               ~ CreateSeuratObject(.x, project = .y) %>% 
                 NormalizeData(verbose = verbose))
  
  method_markers <- methodDE(seus$kallisto, seus$cellranger, test.use = "LR", 
                             latent.vars = "nCount_RNA", verbose = verbose) %>% 
    filter(p_val_adj < 0.05)
  
  # GSEA------------------
  method_markers <- method_markers %>% 
    left_join(gns, by = c("gene" = col_use))
  method_gsea <- egsea_ora(method_markers$gene, species, gns, entrez_inter, kegg_kgsets,
                           logFC = method_markers$avg_logFC) %>% 
    egsea_results_df()
  out <- method_markers %>% 
    left_join(kegg_df, by = "entrezgene") %>% 
    inner_join(method_gsea, by = "gene_set") %>% 
    filter(p.adj < 0.05) %>% 
    mutate(change = cut(avg_logFC, -20:20))
  if (length(species) == 1) {
    out <- out %>% 
      mutate(gene_set = str_remove(gene_set, species_gs(species)))
  }
  # Save the data frame
  write_csv(out, paste0("./output/gsea_bar/", name_use, ".csv"))
  if (nrow(out) == 0) {
    p <- ggplot(out) +
      annotate("text", x = 1, y = 1, label = "No significant gene set") +
      theme_void()
  } else if (length(unique(out$gene_set)) == 1) {
    # Plot histogram if there's one gene set
    p <- ggplot(out, aes(avg_logFC)) +
      geom_histogram(bins = 10) +
      theme_bw() +
      labs(x = paste(unique(out$gene_set), "avg logFC"))
  } else {
    p <- ggplot(out, aes(fct_reorder(gene_set, gene_set, length, .desc = TRUE), 
                         fill = fct_rev(change))) +
      geom_bar() +
      scale_fill_viridis_d(option = "E", name = "logFC", direction = -1) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      labs(y = "Number of genes", x = "gene set")
  }
  p <- p + ggtitle(paste(name_use, paste(species, collapse = ", ")))
  ggsave(paste0("./output/gsea_bar/", name_use, ".png"), 
         p, device = "png", width = 9, height = 6)
}
