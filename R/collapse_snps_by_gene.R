#' Get names of genes from a binary allele matrix
#'
#' @param bin_mat Matrix.
#'
#' @return gene_names: Character. Vector of characters. Length = nrow(bin_mat).
#' @noRd
get_gene_names <- function(bin_mat){
  check_if_binary_matrix(bin_mat)
  variant_names <- row.names(bin_mat)
  gene_names <- gsub(".*[|]", "", variant_names)
  return(gene_names)
}

#' Get names of alleles from a binary allele matrix
#'
#' @param bin_mat Matrix.
#'
#' @return allele_names: Character. Vector of characters. Length = nrow(bin_mat).
#' @noRd
get_allele_names <- function(bin_mat){
  check_if_binary_matrix(bin_mat)
  variant_names <- row.names(bin_mat)
  allele_names <- gsub('\\|.*|\\..*','', variant_names)
  return(allele_names)
}

#' Collapse SNPs into the gene(s) which they are from
#'
#' @param bin_mat Matrix.
#' @param gene_vec Character. Vector of gene names.
#'
#' @return gene_mat: Matrix.
#' @noRd
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

#' Aggregate SNPs by gene and snpeff impact
#'
#' @param num_unique_genes Numeric. number of unique gene names
#' @param unique_gene_names Character vector of unique gene names
#' @param gene_vec passed into collapse_snps_into_genes_by_impact
#' @param bin_mat passed into collapse_snps_into_genes_by_impact
#' @param pred_impact passed into collapse_snps_into_genes_by_impact
#' @param impact character indicating snpeff impact of HIGH, MODERATE, LOW,
#' MODIFIER, or ALL
#'
#' @return list of gene_mats aggregated by impact
#' @noRd
get_gene_mat_by_impact <- function(num_unique_genes, unique_gene_names, gene_vec, bin_mat, pred_impact, impact) {
  gene_mat <- matrix(NA, nrow = num_unique_genes, ncol = ncol(bin_mat))

  if (length(impact) > 1) {
    row.names(gene_mat) <- paste(unique_gene_names, paste(impact,collapse = '-'), sep = '|')
  }else {
    row.names(gene_mat) <- paste(unique_gene_names, impact, sep = '|')
    }
  colnames(gene_mat) <- colnames(bin_mat)

  if (length(impact) == 1) {
    if (impact == 'ALL') {
      impact <- c('MODERATE', 'MODIFIER', 'HIGH', 'LOW')
    }
  }


  for (i in 1:num_unique_genes) {
    current_gene <- unique_gene_names[i]
    temp_snp_mat <- bin_mat[gene_vec == current_gene & pred_impact %in% impact, , drop = FALSE]
    gene_mat[i, ] <- colSums(temp_snp_mat)
    gene_mat[i, ] <-  as.numeric(as.logical(gene_mat[i, ]))
  }

  gene_mat = gene_mat[rowSums(gene_mat) > 0, ,drop = FALSE]
  return(gene_mat)
}

#' Collapse SNPs into the gene(s) which they are from and by snpeff impcat
#'
#' @param bin_mat Matrix.
#' @param gene_vec Character. Vector of gene names.
#' @param predicted_impact Character. Vector of predicted functional impacts.
#' @param snpeff_grouping Character. Vector or single string of impacts of interest.
#' @return a list of gene_mats, collapsed by gene and snpeff impact
#' @noRd
collapse_snps_into_genes_by_impact <- function(bin_mat, gene_vec, predicted_impact, snpeff_grouping){
  check_is_this_class(gene_vec, "character")
  check_if_binary_matrix(bin_mat)

  if (length(gene_vec) != nrow(bin_mat)) {
    stop("gene_vec should have same length as rows of bin_mat")
  }

  unique_gene_names <- unique(gene_vec)
  num_unique_genes <- length(unique_gene_names)

  gene_mat_modifier <- get_gene_mat_by_impact(num_unique_genes, unique_gene_names,
                                              gene_vec, bin_mat, predicted_impact,
                                              'MODIFIER')
  gene_mat_high <- get_gene_mat_by_impact(num_unique_genes, unique_gene_names,
                                          gene_vec, bin_mat, predicted_impact,
                                          'HIGH')
  gene_mat_moderate <- get_gene_mat_by_impact(num_unique_genes, unique_gene_names,
                                              gene_vec, bin_mat, predicted_impact,
                                              'MODERATE')
  gene_mat_low <- get_gene_mat_by_impact(num_unique_genes, unique_gene_names,
                                         gene_vec, bin_mat, predicted_impact, 'LOW')

  gene_mat_all <- get_gene_mat_by_impact(num_unique_genes, unique_gene_names,
                                         gene_vec, bin_mat, predicted_impact,
                                        'ALL')
  if (!is.null(snpeff_grouping)) {
    gene_mat_custom <- get_gene_mat_by_impact(num_unique_genes, unique_gene_names,
                                         gene_vec, bin_mat, predicted_impact,
                                         snpeff_grouping)
  }else{
    gene_mat_custom <- NULL
  }

  return(list(gene_mat_all = gene_mat_all,
              gene_mat_modifier = gene_mat_modifier,
              gene_mat_high = gene_mat_high,
              gene_mat_moderate = gene_mat_moderate,
              gene_mat_low = gene_mat_low,
              gene_mat_custom = gene_mat_custom))

}


