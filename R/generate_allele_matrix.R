replace_non_ATGC_with_N <- function(mat){
  mat <- apply(mat, 2, toupper)
  mat[is.na(mat)] <- "N"
  mat[!(mat %in% c("A", "T", "G", "C"))] <- "N"
  return(mat)
}

identify_variant_sites <- function(mat){
  rows_to_keep <- apply(mat, 1, function(row){
    sum(unique(row) %in% c("A", "T", "G", "C"))
  })
  rows_to_keep <- rows_to_keep > 1
  rows_to_keep <- unname(rows_to_keep)
  return(rows_to_keep)
}

remove_invariant_sites <- function(mat, rows_to_keep){
  mat <- mat[rows_to_keep, , drop = FALSE]
  return(mat)
}

keep_only_variant_sites <- function(dna_mat){
  dna_mat <- replace_non_ATGC_with_N(dna_mat)
  variant_site_log <- identify_variant_sites(dna_mat)
  variant_only_dna_mat <- remove_invariant_sites(dna_mat, variant_site_log)
  return(variant_only_dna_mat)
}
