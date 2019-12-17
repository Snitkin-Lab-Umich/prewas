#' Duplicate SNPs found in overlapping gene pairs.
#'
#' @param bin_mat Matrix. Binary matrix.
#' @param gff_mat Matrix. GFF information.
#'
#' @return bin_mat_dup: Matrix.
#' @export
#'
dup_snps_in_overlapping_genes <- function(bin_mat, gff_mat){
  check_is_this_class(bin_mat, "matrix")
  check_is_this_class(gff_mat, "matrix")

  # position of SNP
  pos <- gsub('[.].*$', '', row.names(bin_mat))

  # what genes are present at the pos?
  genes_at_pos <- lapply(pos, function(p){
    gff_mat[gff_mat[, 4] <= p & gff_mat[, 5] >= p, 9]
  })

  # number of genes per pos
  num_of_genes_at_pos <- sapply(genes_at_pos, length)

  # indices of bin_mat (number of rows)
  indices <- 1:nrow(bin_mat)

  # row names: pos.num.gene_name
  # there will be a .num. if its a multiallelic site
  pos_dup_gene <- paste(rep(row.names(bin_mat), num_of_genes_at_pos),
                        unlist(genes_at_pos),
                        sep = '|')

  # duplicate overlapping genes
  bin_mat_dup = bin_mat[rep(indices, num_of_genes_at_pos), ]
  row.names(bin_mat_dup) = pos_dup_gene

  return(bin_mat_dup)
}
