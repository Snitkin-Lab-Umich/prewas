#' Single line to multiline representation of multiallelic sites
#'
#' @description Changes the single line representation of multiallelic sites to
#'   a multiline representation. That is, rows that have N alleles are
#'   replicated N-1 times. The rows are not yet converted to binary.
#'
#' @param mat Matrix.
#' @param ar_results Data.frame.
#' @param o_ref Character vector. Original reference alleles. Length = Number of
#'   genotypes.
#' @param o_alt Character vector. Original alternative alleles. Length = Number
#'   of genotypes.
#' @param snpeff NULL or list of character vectors. If list, then length(list) =
#'   Number of genotypes. Each list entry can have one or more SnpEff
#'   annotations.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{mat_split}{Matrix. Alleles. Rows refer to genomic loci. Lines
#'      generated from splitting multiline to single line have .1 or .2 appended
#'      to the rowname. Columns are samples.}
#'     \item{ar_results_split}{Data.frame. First column is either the ancestral
#'      allele or the major allele. If ancestral reconstruction was performed
#'      there is a second column, probability. Rows refer to genomic loci. Lines
#'      generated from splitting multiline to single line have .1 or .2 appended
#'      to the rowname.}
#'     \item{split_rows_flag}{Numeric. The count of times a number appears in
#'     the split_rows_flag is one less than the number of alleles at that
#'     locus.}
#'     \item{o_ref}{Character vector. Original reference alleles. Length =
#'     Number of genotypes after splitting.}
#'     \item{o_alt}{Character vector. Original alternative alleles. Length =
#'      Number of genotype after splittings.}
#'     \item{snpeff}{NULL or list of character vectors. If list, then
#'     length(list) = Number of genotypes after splitting. Each list entry can
#'     have one or more SnpEff annotations.}
#'   }
#' @noRd
split_multi_to_biallelic_snps <- function(mat, ar_results, o_ref, o_alt, snpeff){
  check_is_this_class(mat, "matrix")
  check_is_this_class(ar_results, "data.frame")
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

  # GET NUMBER OF UNIQUE ALLELES AT A POSITION
  num_alleles <- apply(mat, 1, function(row) {
    length(unique(row))
  })

  # GET ROW INDICES
  row_indices <- 1:nrow(mat)

  # ROWS THAT HAVE N ALLELES ARE REPLICATED N-1 TIMES
  split_rows_flag <- rep(row_indices, (num_alleles - 1))

  # SPLIT MAT AND AR_RESULTS AND SNPEFF
  mat_split <- mat[split_rows_flag, ]
  ar_results_split <- ar_results[split_rows_flag, , drop = FALSE]
  rownames(mat_split) <- rownames(ar_results_split)
  o_ref_split <- o_ref[split_rows_flag]

  # assign allele to multiallelic site
  o_alt_split <- unlist(sapply(unique(split_rows_flag), function(i) {
    alleles = rep(o_alt[i], sum(split_rows_flag == i))

    sapply(1:length(alleles), function(j) {
      unlist(strsplit(alleles[j], ','))[j]
      })

    }))

  snpeff_split <- snpeff[split_rows_flag]

  return(list(mat_split = mat_split,
              ar_results_split = ar_results_split,
              o_ref_split = o_ref_split,
              o_alt_split = o_alt_split,
              snpeff_split = snpeff_split,
              split_rows_flag = split_rows_flag))
}
