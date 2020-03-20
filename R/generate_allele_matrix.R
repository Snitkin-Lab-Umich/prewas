#' In allele matrix replace non-standard alleles with N
#'
#' @param mat Matrix.
#'
#' @return mat: Matrix.
#' @noRd
replace_non_ATGC_with_N <- function(mat){
  check_is_this_class(mat, "matrix")

  # Make all nucleotides upper case
  mat <- apply(mat, 2, toupper)

  # Replace NAs with "N"
  mat[is.na(mat)] <- "N"


  # For single nucleotide variants, replace non-A/T/G/C with N (ignore long
  # strings like insertions)
  mat[nchar(mat) == 1 & !(mat %in% c("A", "T", "G", "C"))] <- "N"
  return(mat)
}

#' Find all rows in allele matrix with SNPs.
#'
#' Give an error if there are no variant sites in the dataset.
#'
#' @param mat Matrix.
#'
#' @return rows_to_keep: Logical vector.
#' @noRd
identify_variant_sites <- function(mat){
  check_is_this_class(mat, "matrix")
  rows_to_keep <- apply(mat, 1, function(row) {
    sum(unique(row) != "N")
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
#' @noRd
#' TODO add descriptions of o_ref, o_alt, snpeff to BOTH params and return
remove_invariant_sites <- function(mat, o_ref, o_alt, snpeff, rows_to_keep){
  check_is_this_class(mat, "matrix")
  check_is_this_class(rows_to_keep, "logical")
  if (length(rows_to_keep) != nrow(mat)) {
    stop("Logical must have length == nrow(mat)")
  }
  if (sum(rows_to_keep) == 0) {
    stop("There are no valid variant sites in this dataset")
  }
  if (length(o_ref) != length(o_alt)) {
    stop("o_ref and o_alt need to have same length")
  }
  if (length(o_ref) != nrow(mat)) {
    stop("o_ref & o_alt need to have an entry for every mat row")
  }
  if (!is.null(snpeff)) {
    check_is_this_class(snpeff, "list")
    check_is_this_class(snpeff[[1]], "character")
    if (length(snpeff) != length(o_ref)) {
      stop("If not null snpeff needs to have an entry for every row in mat")
    }
  }

  mat <- mat[rows_to_keep, , drop = FALSE]
  o_ref <- o_ref[rows_to_keep]
  o_alt <- o_alt[rows_to_keep]
  snpeff <- snpeff[rows_to_keep]

  return(list('mat' = mat,
              'o_ref' = o_ref,
              'o_alt' = o_alt,
              'snpeff' = snpeff))
}

#' Given a DNA matrix keep only those rows with SNPs.
#'
#' @param dna_mat Matrix.
#'
#' @return variant_only_dna_mat: Matrix.
#' @noRd
#' TODO add descriptions of o_ref, o_alt, snpeff to params
#' TODO add descriptions of o_ref_var_pos, o_alt_var_pos, snpeff_var_pos to return
keep_only_variant_sites <- function(dna_mat, o_ref, o_alt, snpeff){
  check_is_this_class(dna_mat, "matrix")
  if (length(o_ref) != length(o_alt)) {
    stop("Reference and alternative allele vectors must be same length")
  }
  if (nrow(dna_mat) != length(o_ref)) {
    stop("Length of reference and alternative allele vectors must equal nrow(dna_mat)")
  }
  if (!is.null(snpeff)) {
    check_is_this_class(snpeff, "list")
    check_is_this_class(snpeff[[1]], "character")
    if (length(snpeff) != length(o_ref)) {
      stop("If not null snpeff needs to have an entry for every row in mat")
    }
  }

  dna_mat <- replace_non_ATGC_with_N(dna_mat)
  variant_site_log <- identify_variant_sites(dna_mat)
  invariant_sites_removed <- remove_invariant_sites(dna_mat, o_ref, o_alt, snpeff, variant_site_log)
  variant_only_dna_mat <- invariant_sites_removed$mat
  o_ref_var_pos <- invariant_sites_removed$o_ref
  o_alt_var_pos <- invariant_sites_removed$o_alt
  snpeff_var_pos <- invariant_sites_removed$snpeff

  return(list('variant_only_dna_mat' = variant_only_dna_mat,
              'o_ref_var_pos' = o_ref_var_pos,
              'o_alt_var_pos' = o_alt_var_pos,
              'snpeff_var_pos' = snpeff_var_pos))
}
