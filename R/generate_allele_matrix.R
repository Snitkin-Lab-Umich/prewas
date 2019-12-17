#' In allele matrix replace non-standard alleles with N
#'
#' @param mat Matrix.
#'
#' @return mat: Matrix.
#' @export
#'
replace_non_ATGC_with_N <- function(mat){
  check_is_this_class(mat, "mat")
  mat <- apply(mat, 2, toupper)
  mat[is.na(mat)] <- "N"
  mat[!(mat %in% c("A", "T", "G", "C"))] <- "N"
  return(mat)
}

#' Find all rows in allele matrix with SNPs
#'
#' @param mat Matrix.
#'
#' @return rows_to_keep: vector.
#' @export
#'
identify_variant_sites <- function(mat){
  check_is_this_class(mat, "mat")
  rows_to_keep <- apply(mat, 1, function(row) {
    sum(unique(row) %in% c("A", "T", "G", "C"))
  })
  rows_to_keep <- rows_to_keep > 1
  rows_to_keep <- unname(rows_to_keep)
  return(rows_to_keep)
}

#' Remove any invariant sites from allele matrix.
#'
#' @param mat Matrix.
#' @param rows_to_keep Vector. Row numbers to drop.
#'
#' @return mat: Matrix.
#' @export
#'
remove_invariant_sites <- function(mat, rows_to_keep){
  check_is_this_class(mat, "mat")
  check_is_this_class(rows_to_keep, "vector")
  mat <- mat[rows_to_keep, , drop = FALSE]
  return(mat)
}

#' Given a DNA matrix keep only those rows with SNPs.
#'
#' @param dna_mat Matrix.
#'
#' @return variant_only_dna_mat: Matrix.
#' @export
#'
keep_only_variant_sites <- function(dna_mat){
  check_is_this_class(dna_mat, "mat")
  dna_mat <- replace_non_ATGC_with_N(dna_mat)
  variant_site_log <- identify_variant_sites(dna_mat)
  variant_only_dna_mat <- remove_invariant_sites(dna_mat, variant_site_log)
  return(variant_only_dna_mat)
}
