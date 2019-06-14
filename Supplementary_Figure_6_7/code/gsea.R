library(EGSEA)
source("./code/utils.R")
#' Avoid boilerplate in GSEA
#' 
#' This function converts the Ensembl gene IDs into Entrez IDs and then do GSEA
#' with kegg pathways.
#' 
#' @inheritParams EGSEA::buildKEGGIdx
#' @inheritParams EGSEA::egsea.ora
#' @param genes Ensembl IDs of genes for which GSEA is done. Should be a character vector.
#' @param gns A data frame relating Ensembl IDs to Entrez IDs. Ensembl ID should
#' be in a column called "gene", and the Entrez IDs in a column "entrezgene".
#' @return Same thing EGSEA::buildKEGGIdx returns.
egsea_ora <- function(genes, species, gns, universe, gsets, logFC = NULL) {
  if (!is.null(logFC) && (length(genes) != length(logFC))) {
    stop("genes and logFC must have the same length.")
  }
  # In case there's a problem with the annotation
  genes_use <- genes %in% gns$gene
  genes <- genes[genes_use]
  logFC <- logFC[genes_use]
  entrez <- gns$entrezgene[gns$gene %in% genes]
  entrez_use <- !is.na(entrez)
  entrez <- entrez[entrez_use]
  logFC <- logFC[entrez_use]
  gs_annot <- try(buildCustomIdx(entrez, gsets = gsets, species = species, label = "kegg"))
  if (class(gs_annot) == "try-error") return("No geneset")
  egsea.ora(entrez, universe = universe, gs.annots = gs_annot, logFC = logFC,
            report = FALSE, verbose = FALSE)
}

#' Get kegg results data frame from EGSEA results
#' 
#' @param gsa An object of class EGSEAResults, with kegg.
#' @return A data frame.
egsea_results_df <- function(gsa) {
  if (!is(gsa, "EGSEAResults")) {
    df <- data.frame(p.value = NA_real_,
                     p.adj = NA_real_,
                     gene_set = NA_character_)
  } else {
    df <- gsa@results$kegg$test.results$ExperimentalContrast
    df$gene_set <- rownames(df)
    df <- df[,c(1,2,7)]
  }
  df
}

#' Run GSEA for each cluster in each method benchmarked
#' 
#' This function first splits DE output by method for form a list of data frames,
#' and then splits each of these data frames in the list by cluster, and then
#' does GSEA separately for each cluster.
#' 
#' @param markers A data frame output of Seurat's FindAllMarkers, with only the
#' genes to use GSEA. In addition it should have a column called method.
#' @param \dots Arguments passed to egsea_ora (my implementation in this file),
#' except gene and logFC.
#' @return A data frame with columns p.value, p.adj, gene_set, cluster, and method.
cluster_wise_gsea <- function(markers, ...) {
  markers_method <- markers %>% 
    dplyr::select(avg_logFC, cluster, gene, method) %>% 
    group_by(method) %>% 
    group_split()
  markers_cluster <- map(markers_method,
                         ~ .x %>% group_by(cluster) %>% group_split())
  markers_gsea <- map(markers_cluster, function(m) {
    map(m, function(.x) {
      df <- egsea_ora(.x$gene, logFC = .x$avg_logFC, ...) %>% 
        egsea_results_df() 
      df$cluster <- unique(as.character(.x$cluster))
      df$method <- unique(as.character(.x$method))
      df
    }) %>% bind_rows()
  }) %>% bind_rows()
  return(markers_gsea)
}

#' Prepare data frame for some GSEA plots
#' 
#' The GSEA bubble plot and bar plot for logFC require the same set of data frame
#' operations before making the plot, so I'm writing a function for that.
#' 
#' @param markers A data frame output of Seurat's FindAllMarkers, with only the
#' genes used for GSEA. In addition it should have a column called method.
#' @param kegg_df Data frame with two columns: gene_set for kegg gene sets, and
#' entrezgene for Entrez gene IDs.
#' @param gsea_res Data frame output from egsea_results_df, with GSEA done to
#' each cluster, rbinded by method. Must have columns gene_set, method, and cluster.
#' @param species Latin name of species.
#' @return The formatted data frame.
df4gsea_plot <- function(markers, kegg_df, gsea_res, species) {
  if ("cluster" %in% names(markers)) {
    cols_use <- c("gene_set", "method", "cluster")
  } else {
    cols_use <- c("gene_set", "method")
  }
  markers %>% 
    ungroup() %>% 
    left_join(kegg_df, by = "entrezgene") %>% 
    inner_join(gsea_res, by = cols_use) %>% 
    filter(p.adj < 0.05) %>% 
    mutate(gene_set = str_remove(gene_set, species_gs(species)))
}

#' Bubble plot for cluster GSEA results
#' 
#' This function makes a bubble plot with clusters in the x axis and gene sets 
#' in the y axis. Different methods benchmarked will be plotted as different
#' facets. The colors stand for -log10 of adjusted p-value of gene set enrichment,
#' and the size of the points stands for number of marker genes in this gene set.
#' 
#' @inheritParams df4gsea_plot
#' @param ncol Numbers of columns for the facetted plot.
#' @return A ggplot2 object.
gsea_bubble <- function(markers, kegg_df, gsea_res, species, ncol) {
  n_clusts <- max(as.integer(as.character(markers$cluster)))
  markers %>% 
    df4gsea_plot(kegg_df, gsea_res, species) %>% 
    mutate(cluster = factor(cluster, levels = as.character(1:(n_clusts)))) %>% 
    group_by(method, cluster, gene_set, p.adj) %>% 
    dplyr::count() %>% 
    ggplot(aes(cluster, gene_set, color = -log10(p.adj), size = n)) +
    geom_point(alpha = 0.5) +
    scale_color_viridis_c(option = "E") +
    scale_size_continuous(name = "Number of\ngenes") +
    facet_wrap(~ method, ncol = ncol) +
    labs(y = "Gene set")
}

#' Bar chart showing logFC of DE genes in enriched gene sets
#' 
#' This function plots number of DE genes in each interval of logFC in each
#' enriched gene set, facetted by the method benchmarked.
#' 
#' @inheritParams gsea_bubble
#' @param strip_position Where the facet strip should be placed, passed to facet_wrap, 
#' defaults to "right".
#' @param text_angle Angle to tilt x axis labels, useful for long gene set names.
#' @param hjust Horizontal justification of x axis label.
#' @param vjust Vertical justification of y axis label.
#' @param FC_low Lower limit of log fold change. If you see NAs in logFC in the
#' plot, consider adjusting FC_low and FC_high.
#' @param FC_high Upper limit of log fold change.
#' @return A ggplot2 object.
gsea_logFC_bar <- function(markers, kegg_df, gsea_res, species, ncol, 
                           strip_position = "right", text_angle = 45,
                           hjust = 1, vjust = 1,
                           FC_low = -5, FC_high = 5) {
  markers %>% 
    df4gsea_plot(kegg_df, gsea_res, species) %>% 
    mutate(change = cut(avg_logFC, FC_low:FC_high)) %>% 
    ggplot(aes(fct_reorder(gene_set, gene_set, length, .desc = TRUE), fill = fct_rev(change))) +
    geom_bar() +
    scale_fill_viridis_d(option = "E", name = "logFC", direction = -1) +
    facet_wrap(~ method, ncol = ncol, strip.position = strip_position) +
    theme(axis.text.x = element_text(angle = text_angle, hjust = hjust, vjust = vjust)) +
    labs(y = "Number of genes", x = "gene set")
}
