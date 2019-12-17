#' Single line to multiline representation of multiallelic sites
#'
#' @description Changes the single line representation of multiallelic sites to
#'   a multiline representation. That is, rows that have N alleles are
#'   replicated N-1 times.
#'
#' @param mat Matrix.
#' @param ar_results Data.frame.
#'
#' @return A list with three elements:
#'   mat_split: Matrix.
#'   ar_results_split: Data.frame.
#'   split_rows_flag: Numeric.
#' @export
#'
split_multi_to_biallelic_snps = function(mat, ar_results){
  check_is_this_class(mat, "matrix")
  check_is_this_class(ar_results, "data.frame")

  # GET NUMBER OF UNIQUE ALLELES (A, C, T, or G) AT A POSITION
  num_alleles <- apply(mat, 1, function(row) {
    sum(c('A', 'C', 'T', 'G') %in% row)
  })

  # GET ROW INDICES
  row_indices <- 1:nrow(mat)

  # ROWS THAT HAVE N ALLELES ARE REPLICATED N-1 TIMES
  split_rows_flag <- rep(row_indices, (num_alleles-1))

  # SPLIT MAT AND AR_RESULTS
  mat_split <- mat[split_rows_flag, ]
  ar_results_split <- ar_results[split_rows_flag, , drop = FALSE]
  rownames(mat_split) <- rownames(ar_results_split)

  return(list(mat_split = mat_split,
              ar_results_split = ar_results_split,
              split_rows_flag = split_rows_flag))

}
