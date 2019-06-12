library(dplyr)
#' Assign species to cells
#' 
#' This function assigns species to cells in mixed species experiments based on
#' the proportion of UMIs in the cell that belongs to a species. Right now only
#' for human and mice.
#' 
#' @param res_mat Gene count matrix with genes in rows and cells in columns with 
#' Ensembl gene IDs in row names.
#' @param proportion A cell must have at least this proportion of UMIs belonging
#' to a species in order to be assigned to that species.
#' @return A data frame with proportion and number of human and mouse UMIs as
#' well as species assignment.
assign_species <- function(res_mat, proportion = 0.9) {
  gene_species <- ifelse(str_detect(rownames(res_mat), "ENSMUSG"), "mouse", "human")
  mouse_inds <- gene_species == "mouse"
  human_inds <- gene_species == "human"
  # mark cells as mouse or human
  cell_species <- tibble(barcode = colnames(res_mat),
                         n_mouse_umi = Matrix::colSums(res_mat[mouse_inds,]),
                         n_human_umi = Matrix::colSums(res_mat[human_inds,]),
                         tot_umi = Matrix::colSums(res_mat),
                         prop_mouse = n_mouse_umi / tot_umi,
                         prop_human = n_human_umi / tot_umi)
  # Classify species based on proportion of UMI, with cutoff of 90%
  cell_species <- cell_species %>% 
    mutate(species = case_when(
      prop_mouse > 0.9 ~ "mouse",
      prop_human > 0.9 ~ "human",
      TRUE ~ "mixed"
    ))
  return(cell_species)
}
