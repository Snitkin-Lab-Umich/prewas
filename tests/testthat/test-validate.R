# check_is_number -------------------------------------------------------------#
test_that("Check that check_is_number() doesn't give error when given numeric input", {
  # integer
  expect_silent(check_is_number(1))
  # double
  expect_silent(check_is_number(1.5))
  # infinity
  expect_silent(check_is_number(Inf))
})

test_that("Check that check_is_number() gives error when given non-numeric input", {
  # character string
  expect_error(check_is_number("foo"))
  # factor
  expect_error(check_is_number(as.factor("foo")))
  # matrix
  expect_error(check_is_number(matrix(0, 1, 1)))
  # dataframe
  expect_error(check_is_number(as.data.frame(matrix(0, 1, 1))))
})

# is_file ---------------------------------------------------------------------#
test_that("Check that is_file() returns TRUE for files that exist in this directory", {
  # current test file
  expect_true(is_file("test-validate.R"))
})

test_that("Check that is_file() returns FALSE for files that don't exist in this directory", {
  # File that doesn't exist
  expect_false(is_file("fake_file_name.R"))
  # Empty string name
  expect_false(is_file(""))
})

test_that("Check that is_file() returns FALSE for non-character string inputs", {
  # numeric
  expect_false(is_file(1))
  # factor
  expect_false(is_file(as.factor("foo")))
  # matrix
  expect_false(is_file(matrix(0, 1, 1)))
  # dataframe
  expect_false(is_file(as.data.frame(matrix(0, 1, 1))))
})

# is_this_class ---------------------------------------------------------------#
test_that("Check that is_this_class() returns TRUE when inputs match", {
  # numeric
  expect_true(is_this_class(5, "numeric"))
  # character
  expect_true(is_this_class("5", "character"))
  # matrix
  expect_true(is_this_class(matrix(0, 1, 1), "matrix"))
  # factor
  expect_true(is_this_class(as.factor("foo"), "factor"))
})

test_that("Check that is_this_class() returns FALSE when inputs mismatch", {
  # character / numeric
  expect_true(is_this_class("5", "numeric"))
  # numeric / character
  expect_true(is_this_class(5, "character"))
  # matrix / dataframe
  expect_true(is_this_class(matrix(0, 1, 1), "dataframe"))
  # factor / character
  expect_true(is_this_class(as.factor("foo"), "character"))
})

# check_is_this_class ---------------------------------------------------------#
test_that("Check that check_is_this_class() doesn't give error when given matching inputs", {
  # numeric
  expect_silent(check_is_this_class(5, "numeric"))
  # character
  expect_silent(check_is_this_class("5", "character"))
  # matrix
  expect_silent(check_is_this_class(matrix(0, 1, 1), "matrix"))
  # factor
  expect_silent(check_is_this_class(as.factor("foo"), "factor"))
})

test_that("Check that check_is_this_class() gives error when given mismatching inputs", {
  # numeric / character
  expect_false(check_is_this_class(5, "character"))
  # character / numeric
  expect_false(check_is_this_class("5", "numeric"))
  # matrix / dataframe
  expect_false(check_is_this_class(matrix(0, 1, 1), "dataframe"))
  # factor / character
  expect_false(check_is_this_class(as.factor("foo"), "character"))
})

# check_is_tree ---------------------------------------------------------------#
test_that("Check that check_is_tree() doesn't give error when given tree input", {
  test_tree <- ape::rcoal(n = 5)
  expect_silent(check_is_tree(test_tree))
})

test_that("Check that check_is_tree() gives error when given non-tree input", {
  # numeric
  expect_error(check_is_tree(5))
  # character
  expect_error(check_is_tree("foo"))
  # matrix
  expect_error(check_is_tree(matrix(0, 1, 1)))
  # factor
  expect_error(check_is_tree(as.factor("foo")))
})


# check_tree_is_rooted --------------------------------------------------------#
test_that("Check that check_tree_is_rooted() doesn't give error when given a rooted tree input", {
  test_tree <- ape::rcoal(n = 5)
  test_tree <- phytools::midpoint.root(test_tree)
  expect_silent(check_tree_is_rooted(test_tree))
})

test_that("Check that check_tree_is_rooted() gives error when given an unrooted tree input", {
  test_tree <- ape::rcoal(n = 5)
  test_tree <- ape::unroot(test_tree)
  expect_error(check_tree_is_rooted(test_tree))
})

test_that("Check that check_tree_is_rooted() gives error when given a non-tree input", {
  # numeric
  expect_error(check_tree_is_rooted(5))
  # character
  expect_error(check_tree_is_rooted("foo"))
  # matrix
  expect_error(check_tree_is_rooted(matrix(0, 1, 1)))
  # factor
  expect_error(check_tree_is_rooted(as.factor("foo")))
})

# read_gff --------------------------------------------------------------------#
test_that("Check that read_gff() correctly parses input gff file", {
  # TODO this test needs to wait until the test data are finalized (correctly load as package data)
  # TODO as in test_gff <- read_gff(prewas::gff_input)
  # TODO then once test_gff loaded check that gff file has expected dimenions / cds / etc...
})

test_that("Check that read_gff() gives error if given a non-gff file", {
  # TODO this test needs to wait until the test data are finalized (correctly load as package data)
  # TODO as in test_gff <- read_gff(prewas::gff_input)
  # TODO then once test_gff loaded check that gff file has expected dimenions / cds / etc...
})


