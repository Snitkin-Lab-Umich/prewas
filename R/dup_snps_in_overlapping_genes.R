#' Duplicate SNPs found in overlapping gene pairs.
#'
#' @param bin_mat Matrix. Binary matrix.
#' @param gff_mat Matrix. GFF information.
#' @noRd
#' @return bin_mat_dup: Matrix.
dup_snps_in_overlapping_genes <- function(bin_mat, gff_mat){
  check_is_this_class(bin_mat, "matrix")
  check_is_this_class(gff_mat, "matrix")
  check_if_binary_matrix(bin_mat)
  if (ncol(gff_mat) != 9) {
    stop("GFF matrix must have exactly 9 columns")
  }

  # position of SNP
  pos <- gsub("[.].*$", "", row.names(bin_mat))

  # what genes are present at the pos?
  genes_at_pos <- lapply(pos, function(p) {
    gff_mat[gff_mat[, 4] <= p & gff_mat[, 5] >= p, 9]
  })

  # number of genes per position
  num_of_genes_at_pos <- sapply(genes_at_pos, length)

  # indices of bin_mat (number of rows)
  indices <- 1:nrow(bin_mat)

  # row names: pos.num|gene_name
  # there will be a .num. if its a multiallelic site
  pos_dup_gene <- paste(rep(row.names(bin_mat), num_of_genes_at_pos),
                        unlist(genes_at_pos),
                        sep = "|")

  # duplicate overlapping genes
  bin_mat_dup <- bin_mat[rep(indices, num_of_genes_at_pos), ]
  row.names(bin_mat_dup) <- pos_dup_gene

  return(bin_mat_dup)
}

#' Duplicate SNPs found in overlapping gene pairs when there are snpeff
#' annotations in the vcf file.
#'
#' @param bin_mat Matrix. Binary matrix.
#' @param predicted_impact vector of impacts, separated by pipes if pos in
#' overlapping gene
#' @param gene vector of genes, separated by pipes if pos in overlapping gene
#' @noRd
#' @return bin_mat_dup and predicted_impact_split
#'
dup_snps_in_overlapping_genes_snpeff <- function(bin_mat, predicted_impact, gene) {
  check_is_this_class(bin_mat, "matrix")
  check_if_binary_matrix(bin_mat)

  # number of genes per position
  num_ofgenes_at_pos <- sapply(gene, function(g) {
    length(unlist(strsplit(g, '[|]')))
  })
  print(num_ofgenes_at_pos)

  # get gene names
  genes_at_pos <- sapply(gene, function(g) {
    unlist(strsplit(g, '[|]'))
  })

  # snpeff impacts
  predicted_impact_split <- unname(unlist(sapply(predicted_impact, function(g) {
    unlist(strsplit(g, '[|]'))
  })))

  # indices of bin_mat (number of rows)
  indices <- 1:nrow(bin_mat)

  # row names: pos.num|gene_name
  # there will be a .num. if its a multiallelic site
  pos_dup_gene <- paste(rep(row.names(bin_mat), num_ofgenes_at_pos),
                        unlist(genes_at_pos),
                        sep = "|")
  print(pos_dup_gene)

  # duplicate overlapping genes
  bin_mat_dup <- bin_mat[rep(indices, num_ofgenes_at_pos), ]
  row.names(bin_mat_dup) <- pos_dup_gene

  return(list(bin_mat_dup = bin_mat_dup,
              predicted_impact_split = predicted_impact_split))

}
