get_gene_names <- function(bin_mat){
  variant_names <- row.names(bin_mat)
  gene_names <- gsub(".*[|]", "", variant_names)
  return(gene_names)
}

collapse_snps_into_genes <- function(bin_mat, gene_vec){
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


