# replace_non_ATGC_with_N -----------------------------------------------------#
test_that("Check that replace_non_ATGC_with_N() gives error when given invalid input", {
  # integer
  expect_error(replace_non_ATGC_with_N(1))

  # character
  expect_error(replace_non_ATGC_with_N("foo"))

  # factor
  expect_error(replace_non_ATGC_with_N(as.factor("foo")))
})

test_that("Check that replace_non_ATGC_with_N() returns matrix when given valid input", {
  allele_mat <- matrix(c("A", "T", "*", NA), nrow = 3, ncol = 4)

  replace_results <- replace_non_ATGC_with_N(allele_mat)
  # Returns a matrix
  expect_equal(class(replace_results), "matrix")

  # Replaces all NAs
  expect_equal(0, sum(is.na(replace_results)))

  # Replaces all starts
  expect_equal(0, sum(grepl(pattern = "[*]", x = replace_results)))

  # Replaces the correct number of NAs and *s
  expect_equal(6, sum(grepl(pattern = "N", x = replace_results)))
})

# identify_variant_sites ------------------------------------------------------#
test_that("Check that identify_variant_sites() returns correct rows when given valid input", {
  allele_mat <- matrix(c("A", "A", "A", "C", "C", "T"), nrow = 3, ncol = 4)
 # TODO write this when I remember the correct orientation of the allele mat
})

