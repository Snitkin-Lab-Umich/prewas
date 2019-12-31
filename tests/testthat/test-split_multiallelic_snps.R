# split_multi_to_biallelic_snps -----------------------------------------------#
test_that("Check that split_multi_to_biallelic_snps() gives expected output when input with only biallelic sites", {

  # Generate biallelic matrix
  test_mat <- matrix(c("A", "T", "G", "C"), nrow = 4, ncol = 5)
  test_mat[1, 1] <- "G"
  test_mat[2, 1] <- "G"
  test_mat[3, 1] <- "A"
  test_mat[4, 1] <- "G"
  colnames(test_mat) <- c("t1", "t2", "t3", "t4", "t5")
  row.names(test_mat) <- c(1, 11, 21, 31)

  test_ar <- matrix(c("A", "T", "G", "C", .9, .8, .7, .4), nrow = 4, ncol = 2)
  colnames(test_ar) <- c("ancestral_allele", "probability")
  row.names(test_ar) <- row.names(test_mat)
  test_ar <- as.data.frame(test_ar)
  test_ar$probability <- as.numeric(as.character(test_ar$probability))

  test_split <- split_multi_to_biallelic_snps(test_mat, test_ar)
  expect_equal(length(test_split), 3)
  expect_equal(class(test_split$mat_split), "matrix")
  expect_equal(class(test_split$ar_results_split), "data.frame")
  expect_equal(class(test_split$split_rows_flag), "integer")
  expect_equal(length(test_split$split_rows_flag), nrow(test_split$mat_split))
  expect_equal(test_mat, test_split$mat_split)
})

test_that("Check that split_multi_to_biallelic_snps() gives expected output when input with a quadallelic site", {
  # Generate matrix with one multiallelic line
  test_mat <- matrix(c("A", "T", "G", "C"), nrow = 4, ncol = 5)
  test_mat[1, 1] <- "G"
  test_mat[1, 2] <- "C"
  test_mat[1, 3] <- "T"
  test_mat[2, 1] <- "G"
  test_mat[3, 1] <- "A"
  test_mat[4, 1] <- "G"
  colnames(test_mat) <- c("t1", "t2", "t3", "t4", "t5")
  row.names(test_mat) <- c(1, 11, 21, 31)

  test_ar <- matrix(c("A", "T", "G", "C", .9, .8, .7, .4), nrow = 4, ncol = 2)
  colnames(test_ar) <- c("ancestral_allele", "probability")
  row.names(test_ar) <- row.names(test_mat)
  test_ar <- as.data.frame(test_ar)
  test_ar$probability <- as.numeric(as.character(test_ar$probability))

  test_split <- split_multi_to_biallelic_snps(test_mat, test_ar)
  expect_equal(length(test_split), 3)
  expect_equal(class(test_split$mat_split), "matrix")
  expect_equal(class(test_split$ar_results_split), "data.frame")
  expect_equal(class(test_split$split_rows_flag), "integer")
  expect_equal(length(test_split$split_rows_flag), nrow(test_split$mat_split))
  expect_equal(sum(grepl(1, test_split$split_rows_flag)), 3)
  expect_equal(nrow(test_split$mat_split), nrow(test_mat) + 2)
  expect_equal(test_split$mat_split[1, ], test_split$mat_split[2, ])
  expect_equal(test_split$mat_split[2, ], test_split$mat_split[3, ])
})

test_that("Check that split_multi_to_biallelic_snps() gives error given invalid inputs", {
  expect_error(split_multi_to_biallelic_snps("foo", "bar"))
  expect_error(split_multi_to_biallelic_snps(as.data.frame(matrix(0, 1, 1)), "bar"))
  expect_error(split_multi_to_biallelic_snps(matrix(0, 1, 1), "bar"))
})

