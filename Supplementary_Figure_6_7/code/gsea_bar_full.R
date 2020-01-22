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
gsea_bar_full <- function(species, name_regex, mapping,
                          name_use = name_regex, 
                          min_cells = 1000,
                          verbose = TRUE) {
  gns <- map_dfr(species, ~ read_csv(paste0("./data/Ensembl/gns_", 
                                        str_replace(.x, " ", "_"), ".csv")))
  go_df <- map_dfr(species, ~ readRDS(paste0("./data/go/go_", 
                                             str_replace(.x, " ", "_"), ".rds")))
  
  # Matrices----------------
  mats <- list()
  # Read in kallisto matrix
  mats$kallisto <- read_count_output(get_dir(name_regex, 
                                             "/home/single_cell_analysis/kallisto_out_single"),
                                     "genes", tcc = FALSE)
  # Filter barcodes
  bc_rank <- barcodeRanks(mats$kallisto, lower = min_cells)
  mats$kallisto <- mats$kallisto[, Matrix::colSums(mats$kallisto) > 
                                   metadata(bc_rank)$inflection]
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
    col_use <- "ensembl"
  } else if (all(str_detect(genes_use, "^\\d+$"))) {
    col_use <- "entrezgene"
  } else {
    col_use <- "gene_name"
  }
  
  # clean Ensembl gene symbols for mixed species datasets
  if (col_use == "ensembl" && !str_detect(genes_use, "^ENS[A-Z]*G\\d+$")) {
    genes_use <- str_remove(genes_use, "^hg19_")
    genes_use <- str_remove(genes_use, "^mm10_")
    for (i in seq_along(mats)) {
      rownames(mats[[i]]) <- genes_use
    }
  }
  
  # Seurat-----------------
  seus <- imap(mats,
               ~ CreateSeuratObject(.x, project = .y) %>% 
                 NormalizeData(verbose = verbose))
  
  method_markers <- methodDE(seus$kallisto, seus$cellranger, test.use = "LR", 
                             latent.vars = "nCount_RNA", verbose = verbose)
  
  # GSEA------------------
  method_markers <- method_markers %>% 
    left_join(gns, by = c("gene" = col_use))
  if (col_use != "ensembl") {
    names(method_markers)[names(method_markers) == "gene"] <- "tmp"
    names(method_markers)[names(method_markers) == "ensembl"] <- "gene"
    names(method_markers)[names(method_markers) == "tmp"] <- col_use
  }
  if (species == "Arabidopsis thaliana") {
    id_use <- "entrezgene"
    # To deal with some peculiarity of org.At.tair.db, which thinks ensembl ID is entrez
    names(gns)[names(gns) == "entrezgene"] <- "tmp"
    names(gns)[names(gns) == "ensembl"] <- "entrezgene"
    names(gns)[names(gns) == "tmp"] <- "other"
  } else {
    id_use <- "ensembl"
  }
  if (length(mapping) > 1L) { # mixed species, mouse and human
    gs <- rownames(GetAssayData(seus$kallisto))
    universes <- list(human = gs[str_detect(gs, "ENSG")],
                      mouse = gs[str_detect(gs, "ENMUSG")])
    method_gsea <- map2_dfr(mapping, universes, function(m, u) {
      res <- topgo_df(method_markers$gene, method_markers$p_val_adj < 0.05, gns,
                      universe = u, mapping = m, ID = id_use,
                      n_bonferroni = 69,
                      statistic = "fisher") # 69 ontologies tested across 20 datasets in total
      res$mapping <- m
      res
    })
  } else {
    method_gsea <- topgo_df(method_markers$gene, method_markers$p_val_adj < 0.05, gns,
                            universe = rownames(GetAssayData(seus$kallisto)), 
                            mapping = mapping, ID = id_use,
                            n_bonferroni = 69, statistic = "fisher")
  }
  
  method_gsea$dataset <- name_use
  out <- method_markers %>% 
    left_join(go_df, by = c("gene" = "ensembl")) %>% 
    inner_join(method_gsea, by = "GO.ID") %>% 
    dplyr::filter(p_log > uniform_log) %>% 
    mutate(change = cut(avg_logFC, -20:20),
           dataset = name_use)
  # Save the data frame for bar plot
  write_csv(out, paste0("./output/gsea_bar/", name_use, ".csv"))
  write_csv(method_markers, paste0("./output/method_de/", name_use, ".csv"))
  if (nrow(out) == 0) {
    p <- ggplot(out) +
      annotate("text", x = 1, y = 1, label = "No significant gene set") +
      theme_void()
  } else if (length(unique(out$GO.ID)) == 1) {
    # Plot histogram if there's one gene set
    p <- ggplot(out, aes(avg_logFC)) +
      geom_histogram(bins = 10) +
      theme_bw() +
      labs(x = paste(unique(out$Term), "avg logFC"))
  } else {
    max_abs <- abs(out$avg_logFC[which.max(abs(out$avg_logFC))])
    p <- ggplot(out, aes(fct_reorder(Term, Term, length, .desc = TRUE), 
                         fill = fct_rev(change))) +
      geom_bar() +
      scale_fill_viridis_d(option = "E", name = "logFC", direction = -1,
                           begin = 1 - abs(min(out$avg_logFC))/max_abs,
                           end = abs(max(out$avg_logFC))/max_abs) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      labs(y = "Number of genes", x = "gene set")
  }
  p <- p + ggtitle(paste(name_use, paste(species, collapse = ", ")))
  ggsave(paste0("./output/gsea_bar/", name_use, ".png"), 
         p, device = "png", width = 9, height = 6)
  
  # Save the QQ plot data frame
  write_csv(method_gsea, paste0("./output/gsea_qq/", name_use, ".csv"))
  # qq plot
  p <- plot_qq(method_gsea) +
    ggtitle(paste(name_use, paste(species, collapse = ", ")))
  if (length(mapping) > 1L) {
    p <- p +
      facet_wrap(~ mapping)
  }
  ggsave(paste0("./output/gsea_qq/", name_use, ".png"), 
         p, device = "png", width = 9, height = 6)
}
