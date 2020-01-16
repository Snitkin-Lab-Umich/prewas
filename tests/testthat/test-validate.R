# check_is_number -------------------------------------------------------------#
test_that("check_is_number() doesn't give error when given numeric input", {
  # integer
  expect_silent(check_is_number(1))
  # double
  expect_silent(check_is_number(1.5))
  # infinity
  expect_silent(check_is_number(Inf))
})

test_that("check_is_number() gives error when given non-numeric input", {
  # character string
  expect_error(check_is_number("foo"))
  # factor
  expect_error(check_is_number(as.factor("foo")))
  # matrix
  expect_error(check_is_number(matrix(0, 1, 1)))
  # data.frame
  expect_error(check_is_number(as.data.frame(matrix(0, 1, 1))))
})

# is_file ---------------------------------------------------------------------#
test_that("is_file() returns TRUE for files that exist in this directory", {
  # current test file
  expect_true(is_file("test-validate.R"))
})

test_that("is_file() returns FALSE for files that don't exist in this directory", {
  # File that doesn't exist
  expect_false(is_file("fake_file_name.R"))
  # Empty string name
  expect_false(is_file(""))
})

test_that("is_file() returns FALSE for non-character string inputs", {
  # numeric
  expect_false(is_file(1))
  # factor
  expect_false(is_file(as.factor("foo")))
  # matrix
  expect_false(is_file(matrix(0, 1, 1)))
  # data.frame
  expect_false(is_file(as.data.frame(matrix(0, 1, 1))))
})

# is_this_class ---------------------------------------------------------------#
test_that("is_this_class() returns TRUE when inputs match", {
  # numeric
  expect_true(is_this_class(5, "numeric"))
  # character
  expect_true(is_this_class("5", "character"))
  # matrix
  expect_true(is_this_class(matrix(0, 1, 1), "matrix"))
  # factor
  expect_true(is_this_class(as.factor("foo"), "factor"))
})

test_that("is_this_class() returns FALSE when inputs mismatch", {
  expect_false(is_this_class("5", "numeric"))
  expect_false(is_this_class(5, "character"))
  expect_false(is_this_class(matrix(0, 1, 1), "data.frame"))
  expect_false(is_this_class(as.factor("foo"), "character"))
})

# check_is_this_class ---------------------------------------------------------#
test_that("check_is_this_class() doesn't give error when given matching inputs", {
  expect_silent(check_is_this_class(5, "numeric"))
  expect_silent(check_is_this_class("5", "character"))
  expect_silent(check_is_this_class(matrix(0, 1, 1), "matrix"))
  expect_silent(check_is_this_class(as.factor("foo"), "factor"))
})

test_that("check_is_this_class() gives error when given mismatching inputs", {
  expect_error(check_is_this_class(5, "character"))
  expect_error(check_is_this_class("5", "numeric"))
  expect_error(check_is_this_class(matrix(0, 1, 1), "data.frame"))
  expect_error(check_is_this_class(as.factor("foo"), "character"))
})

# check_is_tree ---------------------------------------------------------------#
test_that("check_is_tree() doesn't give error when given tree input", {
  test_tree <- ape::rcoal(n = 5)
  expect_silent(check_is_tree(test_tree))
})

test_that("check_is_tree() gives error when given non-tree input", {
  expect_error(check_is_tree(5))
  expect_error(check_is_tree("foo"))
  expect_error(check_is_tree(matrix(0, 1, 1)))
  expect_error(check_is_tree(as.factor("foo")))
})


# check_tree_is_rooted --------------------------------------------------------#
test_that("check_tree_is_rooted() doesn't give error when given a rooted tree input", {
  test_tree <- ape::rcoal(n = 5)
  test_tree <- phangorn::midpoint(test_tree)
  expect_silent(check_tree_is_rooted(test_tree))
})

test_that("check_tree_is_rooted() gives error when given an unrooted tree input", {
  test_tree <- ape::rcoal(n = 5)
  test_tree <- ape::unroot(test_tree)
  expect_error(check_tree_is_rooted(test_tree))
})

test_that("check_tree_is_rooted() gives error when given a non-tree input", {
  expect_error(check_tree_is_rooted(5))
  expect_error(check_tree_is_rooted("foo"))
  expect_error(check_tree_is_rooted(matrix(0, 1, 1)))
  expect_error(check_tree_is_rooted(as.factor("foo")))
})

# read_gff --------------------------------------------------------------------#
test_that("read_gff() correctly parses input gff file", {
  test_gff <- read_gff(prewas::gff)
  expect_equal(ncol(test_gff), 9)
  expect_true(methods::is(test_gff[, 1], "character"))
  expect_true(methods::is(test_gff, "matrix"))
})

test_that("read_gff() gives error if given invalid input", {
  # Bad file
  expect_error(read_gff("tests/testthat/test-validate.R"))

  # Not file nor matrix
  expect_error(read_gff(prewas::tree))

  # Matrix of wrong dimensions
  expect_error(read_gff(matrix(0, nrow = 5, ncol = 10)))
})

# subset_gff ------------------------------------------------------------------#
test_that("subset_gff() correctly subsets gff to only CDS regions", {
  subsetted_temp_gff <- subset_gff(prewas::gff)
  seq_types <- unique(subsetted_temp_gff[, 3])
  num_seq_types <- length(seq_types)

  expect_equal(num_seq_types, 1)
  expect_true(seq_types == "CDS")
  expect_equal(ncol(subsetted_temp_gff), ncol(prewas::gff))
})

test_that("subset_gff() gives error if given non-matrix input", {
  expect_error(subset_gff(5))
  expect_error(subset_gff("foo"))
})

test_that("subset_gff() gives error if given no CDS regions in gff", {
  temp_gff <- prewas::gff
  temp_gff[, 3] <- "foo"
  # No CDS regions
  expect_error(subset_gff(temp_gff))
})

# clean_up_cds_name_from_gff --------------------------------------------------#
test_that("clean_up_cds_name_from_gff() doesn't change GFF size", {
  cleaned_gff <- clean_up_cds_name_from_gff(prewas::gff)
  expect_equal(ncol(cleaned_gff), ncol(prewas::gff))
})

test_that("clean_up_cds_name_from_gff() gives error for non-GFF input", {
  # Matrix with incorrect dimensions
  expect_error(clean_up_cds_name_from_gff(matrix(0, 1, 1)))
  # Wrong input types
  expect_error(clean_up_cds_name_from_gff(1))
  expect_error(clean_up_cds_name_from_gff("foo"))
})

# load_vcf_file ---------------------------------------------------------------#
test_that("load_vcf_file() works when given vcfR object", {
  vcf_output <- load_vcf_file(prewas::vcf)
  expect_true(methods::is(vcf_output, "matrix"))
  expect_equal(ncol(vcf_output), 14)
  expect_true(methods::is(vcf_output[1, ], "character"))
})

test_that("load_vcf_file() gives error when given non-VCF file input", {
  expect_error(load_vcf_file("data/test-validate.R"))
})

test_that("load_vcf_file() gives error when given non-file, non-vcfR object input", {
  expect_error(load_vcf_file(10))
  expect_error(load_vcf_file(""))
  expect_error(load_vcf_file(prewas::gff))
  expect_error(load_vcf_file(data.frame(as.matrix(0, nrow = 10, ncol = 10))))
})

# check_setequal_tree_mat -----------------------------------------------------#
test_that("check_setequal_tree_mat gives no results when tree$tip.label equals colnames in VCF matrix", {
  vcf_output <- load_vcf_file(prewas::vcf)
  vcf_colnames <- colnames(vcf_output)
  tree_tip_labels <- prewas::tree$tip.label
  expect_silent(check_setequal_tree_mat(tree_tip_labels, vcf_colnames))
})

test_that("check_setequal_tree_mat gives warning when tree$tip.label different than colnames in VCF matrix", {
  vcf_output <- load_vcf_file(prewas::vcf)
  vcf_colnames <- colnames(vcf_output)
  tree_tip_labels <- prewas::tree$tip.label

  expect_error(check_setequal_tree_mat(tree, vcf_colnames))
  expect_error(check_setequal_tree_mat("", vcf_colnames))
  expect_error(check_setequal_tree_mat(tree_tip_labels, ""))
  expect_error(check_setequal_tree_mat(tree_tip_labels, c(1:10)))
})

# check_if_binary_matrix ------------------------------------------------------#
test_that("check_if_binary_matrix() returns nothing when given binary matrix", {
  temp_bin_mat <- matrix(c(0, 1), ncol = 10, nrow = 10)
  expect_silent(check_if_binary_matrix(temp_bin_mat))

  temp_bin_mat <- matrix(1, ncol = 10, nrow = 10)
  expect_silent(check_if_binary_matrix(temp_bin_mat))

  temp_bin_mat <- matrix(0, ncol = 10, nrow = 10)
  expect_silent(check_if_binary_matrix(temp_bin_mat))
})

test_that("check_if_binary_matrix() returns error when given non-binary matrix", {
  temp_bin_mat <- matrix(c(0, 5), ncol = 10, nrow = 10)
  expect_error(check_if_binary_matrix(temp_bin_mat))

  temp_bin_mat <- matrix(LETTERS[1:25], ncol = 5, nrow = 5)
  expect_error(check_if_binary_matrix(temp_bin_mat))

  temp_bin_mat <- matrix(NA, ncol = 10, nrow = 10)
  expect_error(check_if_binary_matrix(temp_bin_mat))
})
