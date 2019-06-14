library(matrixStats)
#' Convert Latin species name into regex to clean gene set names
#' 
#' Kegg gene set names have something like "hsa00010" in front of them to 
#' indicate species, but that doesn't look good on plots. That's why I wrote this
#' function, to remove those prefixes and make gene set names look nicer on plots.
#' 
#' @param species Latin name of species. If length is more than 1, then only
#' the first element will be used.
#' @return A regex that can be used in str_remove to remove those prefixes.
species_gs <- function(species) {
  if (length(species) > 1) {
    message("Only first element of species is used.")
  }
  species <-  str_split(tolower(species), pattern = " ")[[1]]
  paste0("^", str_sub(species[1], 1, 1), str_sub(species[2], 1, 2), "\\d+ ")
}

#' Converrt SingleR output into data frame
#' 
#' This function grabs the most important information from SingleR output and
#' puts it into a data frame. As SingleR is based on correlation with bulk RNA-seq
#' data of pure cell types, this is more than just a benchmark on cell type
#' annotation.
#' 
#' @param singler SingleR output.
#' @param cluster A factor or character vector for cluster of each barcode.
#' Should be in the same order as the barcodes.
#' @return A data frame with 4 columns: label for cell type assigned, barcode,
#' cluster for the cluster of the cell from Seurat, and score for the maximum
#' score for one cell among all cell types.
singler2df <- function(singler, cluster) {
  tibble(label = singler$labels[,1],
         barcode = singler$cell.names,
         cluster = cluster,
         max_score = rowMaxs(singler$scores))
}

importCDS <- function (otherCDS, seurat_scale=FALSE, import_all = FALSE) 
{
  if (class(otherCDS)[1] == "Seurat") {
    requireNamespace("Seurat")
    if (!seurat_scale) {
      data <- otherCDS@assays$RNA@counts
    } else {
      data <- otherCDS@assays$RNA@scale.data
    }
    if (class(data) == "data.frame") {
      data <- as(as.matrix(data), "sparseMatrix")
    }
    pd <- tryCatch({
      pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
      pd
    }, error = function(e) {
      pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd <- new("AnnotatedDataFrame", data = pData)
      message("This Seurat object doesn't provide any meta data")
      pd
    })
    if (length(setdiff(colnames(data), rownames(pd))) > 0) {
      data <- data[, rownames(pd)]
    }
    fData <- data.frame(gene_short_name = row.names(data), 
                        row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    #lowerDetectionLimit <- otherCDS@is.expr
    if (all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
      expr <- "negbinomial.size"
    }
    else if (any(data < 0)) {
      expressionFamily <- uninormal()
      expr <- "unimormal"
    }
    else {
      expressionFamily <- tobit()
      expr <- "tobit"
    }
    print(paste0("expressionFamily ",expr))
    # valid_data <- data[, row.names(pd)]
    monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
                                  #lowerDetectionLimit = lowerDetectionLimit,
                                  expressionFamily = expressionFamily)
    if (import_all) {
      if ("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
      }
      else {
        mist_list <- otherCDS
      }
    }
    else {
      mist_list <- list()
    }
    if ("var.genes" %in% slotNames(otherCDS)) {
      var.genes <- setOrderingFilter(monocle_cds, otherCDS@var.genes)
    }
    monocle_cds@auxClusteringData$seurat <- mist_list
  }
  else if (class(otherCDS)[1] == "SCESet") {
    requireNamespace("scater")
    message("Converting the exprs data in log scale back to original scale ...")
    data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset
    fd <- otherCDS@featureData
    pd <- otherCDS@phenoData
    experimentData = otherCDS@experimentData
    if ("is.expr" %in% slotNames(otherCDS)) 
      lowerDetectionLimit <- otherCDS@is.expr
    else lowerDetectionLimit <- 1
    if (all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    }
    else if (any(data < 0)) {
      expressionFamily <- uninormal()
    }
    else {
      expressionFamily <- tobit()
    }
    if (import_all) {
      mist_list <- otherCDS
    }
    else {
      mist_list <- list()
    }
    monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
                                  lowerDetectionLimit = lowerDetectionLimit, expressionFamily = expressionFamily)
    monocle_cds@auxOrderingData$scran <- mist_list
  }
  else {
    stop("the object type you want to export to is not supported yet")
  }
  return(monocle_cds)
}
