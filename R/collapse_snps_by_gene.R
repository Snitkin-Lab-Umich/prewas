#' Get names of genes from a binary allele matrix
#'
#' @param bin_mat Matrix.
#'
#' @return gene_names: Character. Vector of characters. Length = nrow(bin_mat).
#'
get_gene_names <- function(bin_mat){
  check_is_this_class(bin_mat, "matrix")
  check_if_binary_matrix(bin_mat)

  variant_names <- row.names(bin_mat)
  gene_names <- gsub(".*[|]", "", variant_names)
  return(gene_names)
}

#' Collapse SNPs into the gene(s) which they are from
#'
#' @param bin_mat Matrix.
#' @param gene_vec Character. Vector of gene names.
#'
#' @return gene_mat: Marix.
#'
collapse_snps_into_genes <- function(bin_mat, gene_vec){
  check_is_this_class(gene_vec, "character")
  check_if_binary_matrix(bin_mat)

  if (length(gene_vec) != nrow(bin_mat)) {
    stop("gene_vec should have same length as rows of bin_mat")
  }
  unique_gene_names <- unique(gene_vec)
  num_unique_genes <- length(unique_gene_names)

  gene_mat <- matrix(NA, nrow = num_unique_genes, ncol = ncol(bin_mat))
  row.names(gene_mat) <- unique_gene_names
  colnames(gene_mat) <- colnames(bin_mat)

  for (i in 1:num_unique_genes) {
    current_gene <- unique_gene_names[i]
    temp_snp_mat <- bin_mat[gene_vec == current_gene, , drop = FALSE]
    gene_mat[i, ] <- colSums(temp_snp_mat)
    gene_mat[i, ] <-  as.numeric(as.logical(gene_mat[i, ]))
  }

  return(gene_mat)
}
