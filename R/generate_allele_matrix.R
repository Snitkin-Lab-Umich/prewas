#' In allele matrix replace non-standard alleles with N
#'
#' @param mat Matrix.
#'
#' @return mat: Matrix.
#' @export
#'
replace_non_ATGC_with_N <- function(mat){
  check_is_this_class(mat, "matrix")
  mat <- apply(mat, 2, toupper)
  mat[is.na(mat)] <- "N"
  mat[!(mat %in% c("A", "T", "G", "C"))] <- "N"
  return(mat)
}

#' Find all rows in allele matrix with SNPs.
#'
#' Give an error if there are no variant sites in the dataset.
#'
#' @param mat Matrix.
#'
#' @return rows_to_keep: Logical vector.
#' @export
#'
identify_variant_sites <- function(mat){
  check_is_this_class(mat, "matrix")
  rows_to_keep <- apply(mat, 1, function(row) {
    sum(unique(row) %in% c("A", "T", "G", "C"))
  })
  rows_to_keep <- rows_to_keep > 1
  rows_to_keep <- unname(rows_to_keep)

  if (sum(rows_to_keep) == 0) {
    stop("There are no valid variant sites in this dataset")
  }
  return(rows_to_keep)
}

#' Remove any invariant sites from allele matrix.
#'
#' @param mat Matrix.
#' @param rows_to_keep Logical. Vector of TRUE/FALSE. FALSE corresponds to rows
#'   that should be dropped (invariant loci). TRUE corresponds to rows to keep
#'   (variant loci).
#'
#' @return mat: Matrix.
#' @export
#'
remove_invariant_sites <- function(mat, rows_to_keep){
  check_is_this_class(mat, "matrix")
  check_is_this_class(rows_to_keep, "logical")

  if (length(rows_to_keep) != nrow(mat)) {
    stop("Logical must have length == nrow(mat)")
  }

  if (sum(rows_to_keep) == 0) {
    stop("There are no valid variant sites in this dataset")
  }

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
  check_is_this_class(dna_mat, "matrix")
  dna_mat <- replace_non_ATGC_with_N(dna_mat)
  variant_site_log <- identify_variant_sites(dna_mat)
  variant_only_dna_mat <- remove_invariant_sites(dna_mat, variant_site_log)
  return(variant_only_dna_mat)
}
