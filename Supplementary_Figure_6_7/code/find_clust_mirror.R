library(dplyr)
library(purrr)
#' Find proportion of shared markers bewteen two clusters
#' 
#' @param markers Data frame output from FindAllMarkers in Seurat.
#' @param clust1 String, name of one of the two clusters.
#' @param clust2 String, name of the other cluster.
#' @return Logical, whether the proportion of shared markers is greater than the
#' threshold.
if_shared <- function(markers, clust1, clust2, threshold = 0.5) {
  df <- markers %>% 
    filter(p_val_adj < 0.05, cluster %in% c(clust1, clust2)) %>% 
    group_by(gene) %>% 
    dplyr::count() %>% 
    group_by(n > 1) %>% 
    dplyr::count()
  if (nrow(df) == 1) {
    if (df[["n > 1"]]) return(TRUE) else return(FALSE)
  } else {
    return((df$n[df[["n > 1"]]] / sum(df$n)) > threshold)
  }
}

#' Find mirroring clusters across methods
#' 
#' This function finds clusters mirrored across kallisto bus and CellRanger
#' outputs by looking at marker genes shared.
#' 
#' @param markers Data frame output from FindAllMarkers in Seurat, with marker
#' genes for individual clusters rather than methods in general. Here the
#' cluster names are assumed to be named with a number (larger number for smaller
#' cluster as in Seurat's cluster naming convention) followed by the method 
#' (kallisto_bus or CellRanger), separated by one underscore. 
#' @param markers_method Data frame output from FindAllMarkers in Seurat, with
#' markers for kallisto bus vs CellRanger in general, not for individual clusters.
#' This argument is used to remove genes that mark methods in general.
#' @param threshold At least this proportion of markers must be shared in order
#' for two clusters to be considered shared. Defaults to 0.5.
#' @param max_ind_diff The cluster index from Seurat is ordered in such a way
#' that smaller index means larger clusters. For instance, cluster 0 is the largest,
#' followed by 1, 2, and etc. Mirrored clusters should have similar size, so to
#' avoid unnecessary work, search is restricted to indices within a certain
#' range from the current index of interest. Defaults to 3.
#' @return A data frame with two columns. One column for clusters for kallisto bus,
#' and the other column for corresponding columns for CellRanger.
#' 
find_clust_mirror <- function(markers, markers_method, threshold = 0.5,
                              max_ind_diff = 3L) {
  # Remove genes that mark methods in general
  markers <- markers %>% 
    mutate_if(is.factor, as.character) %>% 
    anti_join(markers_method, by = "gene")
  clust_kb <- unique(markers$cluster[str_detect(markers$cluster, "kallisto_bus$")])
  clust_cr <- unique(markers$cluster[str_detect(markers$cluster, "CellRanger$")])
  out <- cross_df(list(kallisto_bus = clust_kb,
                       CellRanger = clust_cr)) %>% 
    separate(kallisto_bus, into = c("ind_kb", "kb"), sep = "_", remove = FALSE,
             extra = "drop") %>% 
    separate(CellRanger, into = c("ind_cr", "cr"), sep = "_", remove = FALSE) %>% 
    dplyr::select(-kb, -cr) %>% 
    mutate(ind_kb = as.integer(ind_kb),
           ind_cr = as.integer(ind_cr)) %>% 
    filter(abs(ind_kb - ind_cr) <= max_ind_diff)
  out <- out %>% 
    mutate(is_mirror = map2_lgl(kallisto_bus, CellRanger,
                                if_shared, markers = markers, 
                                threshold = threshold)) %>% 
    filter(is_mirror) %>% 
    dplyr::select(-ind_kb, -ind_cr, -is_mirror) %>% 
    distinct()
  return(out)
}

#' Reformat output of FindMarkers
#' 
#' I hate rownames, and no wonder setting rownames is deprecated in Tidyverse.
#' But probably the Seurat team disagrees, so the output of FindMarkers has
#' genes as rownames rather than a separate column. Here I reformat that output
#' so genes are a column on its own.
#' 
#' @param markers Output from FindMarkers in Seurat.
#' @return A data frame with genes as a column and without rownames.
#' 
clean_markers <- function(markers) {
  markers$gene <- rownames(markers)
  rownames(markers) <- NULL
  markers
}

#' Add mean log normalized expression of genes to markers
#' 
#' For MA plots and future reference.
#' 
#' @param markers Output from FindAllMarkers in Seurat. Must have gene column.
#' @param seu Seurat object where markers were found.
#' @return A data frame, like \code{markers}, with mean log count for each gene
#' of interest.
add_marker_log_mean <- function(markers, seu) {
  gene_log_mean <- tibble(gene = unique(markers$gene),
                          log_norm_mean = Matrix::rowMeans(seu@data[gene, ]))
  markers %>% 
    left_join(gene_log_mean, by = "gene")
}
