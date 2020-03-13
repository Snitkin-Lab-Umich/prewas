#' Single line to multiline representation of multiallelic sites
#'
#' @description Changes the single line representation of multiallelic sites to
#'   a multiline representation. That is, rows that have N alleles are
#'   replicated N-1 times. The rows are not yet converted to binary.
#'
#' @param mat Matrix.
#' @param ar_results Data.frame.
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
#'     the split_rows_flag is one less than the number of alleles at that locus.
#'     }
#'   }
#' @noRd
split_multi_to_biallelic_snps <- function(mat, ar_results, snpeff){
  check_is_this_class(mat, "matrix")
  check_is_this_class(ar_results, "data.frame")

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
  snpeff_split <- snpeff[split_rows_flag]

  return(list(mat_split = mat_split,
              ar_results_split = ar_results_split,
              snpeff_split = snpeff_split,
              split_rows_flag = split_rows_flag))
}
