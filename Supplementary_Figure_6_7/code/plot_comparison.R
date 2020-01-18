library(zeallot)
library(gridExtra)
library(ggpubr)
library(Seurat)
library(ggrepel)
library(OmicsMarkeR)
library(tidyverse)
#' Knee plot
#' 
#' This function makes the knee plot with total counts in the y axis and barcode
#' rank in the x axis along with annotations showing knee point and inflection
#' point. This function is specific to this project. It can be used to compare
#' kallisto bus and CellRanger results.
#' 
#' @param bc_ranks A named list of barcode ranks, from the \code{barcodeRanks}
#' function.
#' @return A ggplot2 object.
knee_plot <- function(bc_ranks) {
  knee_plt <- tibble(rank = map(bc_ranks, ~ .x[["rank"]]),
                     total = map(bc_ranks, ~ .x[["total"]]),
                     method = names(bc_ranks)) %>% 
    unnest() %>% 
    distinct() %>% 
    filter(total > 0)
  annot <- tibble(inflection = map_dbl(bc_ranks, ~ metadata(.x)[["inflection"]]),
                  method = names(bc_ranks))
  p <- ggplot(knee_plt, aes(rank, total, color = method)) +
    geom_line() +
    geom_hline(aes(yintercept = inflection, color = method), 
               data = annot, linetype = 2) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Rank", y = "Total reads")
  return(p)
}

#' Connect points with the same barcode across methods
#' 
#' This function is used for comparing kallisto bus and CellRanger, to see how
#' the change in method shifts cells in gene expression space.
#' 
#' @param seu Seurat object with dimension reductions done. There must be a
#' column in the metadata called "method" indicating whether a cells is from
#' kallisto bus or CellRanger.
#' @param dr String for dimension reduction to use. Should be name of an dr 
#' object in the dr slot of the Seurat object, such as "pca".
#' @param pt_size Point size.
#' @param alpha Transparency of points.
#' @param seg_size Line width of the segments.
#' @param seg_alpha Transparency of segments.
#' @param seg_color Color of segments.
#' @param seg_type Line type of segments.
#' @param dim.1 Numeric, the first dimension to plot.
#' @param dim.2 Numeric, the second dimension to plot.
#' @param use_species Logical, whether the dataset has mixed species. If TRUE,
#' then there must be a column called "species" in the metadata.
#' @return A ggplot2 object.
connect_plot <- function(seu, dr, pt_size = 1, pt_alpha = 0.5, 
                         seg_size = 0.3, seg_alpha = 0.1, seg_color = "blue",
                         seg_type = 1, dim1 = 1, dim2 = 2, 
                         use_species = FALSE) {
  pca_plt <- tibble(x = seu@dr[[dr]]@cell.embeddings[, dim1],
                    y = seu@dr[[dr]]@cell.embeddings[, dim2],
                    method = seu@meta.data$method,
                    barcode = colnames(seu@data)) %>% 
    mutate(barcode = str_remove(barcode, "(_kb)|(_cr)"))
  if (use_species) {
    pca_plt$species <- seu@meta.data$species
    pca_plt2 <- pca_plt %>% 
      select(-species)
  } else pca_plt2 <- pca_plt
  c(kb_pca, cr_pca) %<-% (pca_plt2 %>% 
                            split(pca_plt2$method))
  pca_seg <- kb_pca %>% 
    full_join(cr_pca, by = "barcode")
  p <- ggplot(pca_plt)
  if (use_species) {
    p <- p +
      geom_point(aes(x, y, color = method, shape = species), 
                 alpha = alpha, size = pt_size) +
      scale_shape_manual(values = c(16, 5, 8))
  } else {
    p <- p +
      geom_point(aes(x, y, color = method), alpha = pt_alpha, size = pt_size)
  }
  dr1 <- ifelse(dr == "pca", "PC", dr)
  p <- p +
    geom_segment(data = pca_seg, 
                 aes(x = x.x, y = y.x, xend = x.y, yend = y.y),
                 size = seg_size, alpha = seg_alpha, color = seg_color,
                 linetype = seg_type) +
    labs(x = paste0(dr1, dim1), y = paste0(dr1, dim2))
  return(p)
}

#' Plot expression of top n genes in cluster
#' 
#' @param markers Data frame returned by FindAllMarkers in Seurat.
#' @param clust_use String, which cluster's top marker genes to be plotted.
#' @param seu Seurat object from which the clustering and differential expression
#' were done.
#' @param n How many top genes to plot. Genes are sorted by fold change in expression
#' of the cluster of interest compared to everything else. If negative, will
#' plot the most depleted n genes.
#' @param \dots Arguments to be passed to FeaturePlot in Seurat.
#' @return No return value, only a graphical output.
plot_top_n <- function(markers, clust_use, seu, n = 4, ...) {
  if (n > 0) f <- function(x) x else f <- function(x) -x
  df <- markers %>% 
    filter(cluster == clust_use, p_val_adj < 0.05) %>% 
    top_n(n, avg_logFC) %>% 
    arrange(desc(f(avg_logFC)))
  FeaturePlot(seu, features.plot = df$gene, ...)
}

#' Plot nGene vs nUMI
#' 
#' Another QC plot for total counts. Axes are in log scale.
#' 
#' @param seus A list of Seurat object, with columns nGene and nUMI in metadata.
#' @param pt_size Point size.
#' @param pt_alpha Point transparency.
#' @param compare_methods Logical, whether comparing kallisto_bus and CellRanger.
#' If TRUE, then there should be a column called method in metadata.
#' @return A ggplot2 object.
gene_umi_plot <- function(seus, pt_size = 0.5, pt_alpha = 0.5) {
  metas <- map(seus, ~ .x@meta.data) %>% 
    bind_rows()
  ggplot(metas, aes(nUMI, nGene)) +
    geom_point(size = pt_size, alpha = pt_alpha) +
    scale_x_log10() +
    scale_y_log10()
}

#' MA plot
#' 
#' Plot mean of normalized counts in x axis and log fold change between two
#' conditions in the y axis, highlighting genes with reasonably large log fold 
#' change.
#' 
#' @param markers Output of FindMarkers (after making gene a
#' column from \code{clean_markers}). Must have a column log_norm_mean for mean 
#' count in log normalized data for the genes of interest (see 
#' \code{add_marker_log_mean}).
#' @param gene_labels Character, column name in \code{markers} with gene names
#' for labeling. Ignored if \code{do_label = FALSE}.
#' @param logFC_threshold Positive number, absolute value of LogFC threshold 
#' for coloring.
#' @param do_label Logical, whether to label DE genes.
#' @param pt_size Point size.
#' @param pt_alpha Point transparency.
#' @param color_ns Color for genes that don't have enough logFC.
#' @param color_sig Color for genes with enough logFC
#' @param n_label Number of genes to label. If NULL, then all genes, unless 
#' there're too many genes to label. At most top 20 genes will be labeled if NULL.
#' @return A ggplot2 object.
plotMA <- function(markers, gene_labels = "gene", logFC_threshold = 1,
                   do_label = TRUE, pt_size = 0.5, pt_alpha = 0.5,
                   color_ns = "gray30", color_sig = "blue", n_label = NULL) {
  markers <- markers %>% 
    mutate(is_sig = abs(avg_logFC) > logFC_threshold)
  if (!any(markers$is_sig)) {
    p <- ggplot(markers, aes(log_norm_mean, avg_logFC)) +
      geom_point(size = pt_size, alpha = pt_slpha, color = color_ns)
  } else {
    if (max(markers$avg_logFC) < logFC_threshold) {
      thresh_use <- -logFC_threshold
    } else if (min(markers$avg_FC) > -logFC_threshold) {
      thresh_use <- logFC_threshold
    } else {
      thresh_use <- c(logFC_threshold, -logFC_threshold)
    }
    p <- ggplot(markers, aes(log_norm_mean, avg_logFC, color = is_sig)) +
      geom_point(size = pt_size, alpha = pt_slpha) +
      geom_hline(yintercept = thresh_use, linetype = 2, color = "gray60") +
      scale_color_manual(c(color_ns, color_sig)) +
      theme(legend.position = "none")
    if (do_label) {
      n_sig <- sum(markers$is_sig)
      if (is.null(n_label)) {
        if (n_sig > 20) {
          n_label <- 20
        } else {
          n_label <- n_sig
        }
      } else if (n_label > n_sig) {
        n_label <- n_sig
      }
      ord <- order(abs(markers$avg_logFC), decreasing = TRUE)[1:n_label]
      gl <- markers[ord, c("avg_logFC", "log_norm_mean", gene_labels)]
      p <- p +
        geom_text_repel(aes_string("log_norm_mean", "avg_logFC", 
                                   label = gene_labels),
                        data = gl)
    }
  }
  p <- p +
    labs(x = "Mean of log normalized counts", y = "log2 fold change")
  return(p)
}

#' Plot number of markers per cluster
#' 
#' Bar plot of number of markers per cluster that has at least a logFC. Also includes
#' a text number on top of the bar.
#' 
#' @param markers Output of FindAllMarkers.
#' @param logFC_theshold Only marker genes with at least this logFC will be
#' considered in the plot.
#' @param nudge_y Y offset of text.
#' @param positive_only Logical, whether to only look at genes enriched in
#' a cluster.
#' @return A ggplot2 object.
plot_n_markers <- function(markers, ncol = 2, nrow = ceiling(length(seus) / ncol),
                           logFC_threshold = 0.5, nudge_y = 10,
                           positive_only = FALSE) {
  if (positive_only) f <- function(x) x else f <- abs
  n_markers <- markers %>% 
    filter(p_val_adj < 0.05, f(avg_logFC) > logFC_threshold) %>% 
    group_by(method, cluster) %>% 
    dplyr::count()
  ggplot(n_markers, aes(cluster, n)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = n), nudge_y = nudge_y) +
    facet_wrap(~ method, ncol = ncol, nrow = nrow)
}

#' Compute matrix of pairwise Jaccard indices
#' 
#' When comparing more than 2 sets, each pair has a Jaccard index, and we need
#' a matrix of Jaccard indices.
#' 
#' @param set_list Named list of sets.
#' @param plot Logical, whether a plot of the jaccard indices should be shown.
#' @return A data frame with each unique pair of different sets along with their
#' jaccard index. The data frame is returned invisibly.
get_jaccard_matrix <- function(set_list, plot = TRUE) {
  out <- cross_df(list(set1 = seq_along(set_list),
                       set2 = seq_along(set_list)),
                  .filter = `>=`) %>% 
    mutate(jaccard = map2_dbl(set1, set2, ~ jaccard(set_list[[.x]], set_list[[.y]])),
           set1 = names(set_list)[set1],
           set2 = names(set_list)[set2])
  if (plot) {
    p <- ggplot(out, aes(set1, set2, fill = jaccard)) +
      geom_tile() +
      geom_text(aes(label = round(jaccard, digits = 2)), color = "gray") +
      coord_equal() +
      scale_x_discrete(expand = c(0,0)) +
      scale_y_discrete(expand = c(0,0))
    print(p)
  }
  invisible(out)
}

#' Plot dimension reductions for a list of Seurat objects
#' 
#' To avoid boiler plate.
#' 
#' @param seus A list of Seurat object all with the dimension reduction of interest.
#' @param reduction The dimension reduction to use, such as pca or tsne.
#' @param ncol Number of columns in the final plot.
#' @param nrow Number of rows in the final plot.
#' @param \dots Other arguments to pass to Seurat::DimPlot.
#' @return An object of class ggarrange.
DimPlot_list <- function(seus, reduction = "pca", dims = c(1,2), ncol = 2, 
                         nrow = ceiling(length(seus) / ncol), arrange_grob = FALSE, ...) {
  if (reduction == "pca") {
    axis_labels <- paste("PC", dims)
  } else if (reduction == "tsne") {
    axis_labels <- paste("t-SNE", dims)
  } else {
    axis_labels <- paste(toupper(reduction), dims)
  }
  plots <- imap(seus, 
                function(.x, .y) {
                  DimPlot(.x, reduction = reduction, dims = dims, ...) +
                    labs(x = axis_labels[1], y = axis_labels[2], title = .y) +
                    theme_bw() +
                    theme(legend.position = "none")
                })
  ggarrange(plotlist = plots, ncol = ncol, nrow = nrow)
}

#' Compare row or column sums between intersecting and non-intersecting
#' 
#' This is to see whether the non-intersecting genes or barcodes have different
#' total counts from intersecting genes or barcodes.
#' 
#' @param inter Character vector of intersecting genes or barcodes.
#' @param diff Character vector of non-intersecting genes or barcodes. This means
#' those not in the intersection of the 4 methods.
#' @param sums A named list of column or row sums from matrices. I'm not calculating
#' this inside this function as this can be reused.
#' @param dim Which dimension. If 1, then row, and if 2, then column.
#' @return A ggplot2 object.
#' 
intersect_count_plot <- function(inter, diff, sums, ncol = 2) {
  tot_inter <- map(sums, ~ .x[inter]) %>% 
    as_tibble() %>% 
    mutate(category = "intersecting")
  tot_diff <- map(sums, ~ .x[diff]) %>% 
    as_tibble() %>% 
    mutate(category = "non-intersecting")
  tot <- bind_rows(tot_inter, tot_diff) %>% 
    gather(key = "method", value = "value", -category) %>% 
    filter(!is.na(value)) %>% 
    distinct()
  ggplot(tot, aes(category, value)) +
    geom_violin() +
    geom_jitter(alpha = 0.1, size = 0.1) +
    facet_wrap(~ method, ncol = ncol) +
    scale_y_log10()
}

#' Differential expression between two methods
#' 
#' This function does DE for the same cells and genes but from different methods,
#' with all clusters together.
#' 
#' @param seu1 Seurat object.
#' @param seu2 Seurat object.
#' @param \dots Other arguments passed to \code{Seurat::FindMarkers}.
#' @return A data frame like the output of \code{Seurat::FindMarkers}.
methodDE <- function(seu1, seu2, ...) {
  # Make cell names unique
  seu1 <- RenameCells(seu1, add.cell.id = "a")
  seu2 <- RenameCells(seu2, add.cell.id = "b")
  seu <- merge(seu1, seu2)
  Idents(seu) <- "orig.ident"
  markers <- FindMarkers(seu, ident.1 = unique(as.character(seu1$orig.ident)), ...) %>% 
    clean_markers()
}

#' Differential expression between two methods (cluster-wise)
#' 
#' This function does DE for the same cells and genes but from different methods
#' cluster-wise. The clustering from seu1 will be used, and the comparison will be done on the
#' same barcodes in seu1 and seu2.
#' 
#' @param seu1 Seurat object. Clusters from this object will be used.
#' @param seu2 Seurat object.
#' @param \dots Other arguments passed to \code{Seurat::FindMarkers}.
#' @return A data frame like the output of \code{Seurat::FindMarkers}, with an
#' extra column "cluster" indicating which cluster in seu1.
methodDE_cluster <- function(seu1, seu2, ...) {
  clusts <- unique(as.character(seu1$RNA_snn_res.1))
  Idents(seu1) <- "RNA_snn_res.1"
  dfs <- map(clusts, function(x) {
    seu1_sub <- subset(seu1, idents = x)
    seu2_sub <- subset(seu2, cells = Cells(seu1_sub))
    df <- methodDE(seu1_sub, seu2_sub, ...)
    df$cluster <- x
    df
  }) %>% bind_rows()
  dfs
}
