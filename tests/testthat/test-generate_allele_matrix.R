# replace_non_ATGC_with_N -----------------------------------------------------#
test_that("Check that replace_non_ATGC_with_N() gives error when given valid input", {
  # integer
  expect_error(replace_non_ATGC_with_N(1))
  # character
  expect_silent(replace_non_ATGC_with_N("foo"))
  # factor
  expect_silent(replace_non_ATGC_with_N(as.factor("foo")))
})


