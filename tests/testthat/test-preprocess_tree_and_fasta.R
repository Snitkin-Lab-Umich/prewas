# preprocess_tree_and_fasta ---------------------------------------------------#
# read_in_tree ----------------------------------------------------------------#
test_that("Check that read_in_tree() works for valid tree input", {
  test_tree <- read_in_tree(prewas::tree)
  expect_equal(ape::Ntip(test_tree), ape::Ntip(prewas::tree))
})

test_that("Check that read_in_tree() gives error for invalid inputs", {
  expect_error(read_in_tree(10))
  expect_error(read_in_tree("tests/testthat/test-preprocess_tree_and_fasta.R"))
  expect_error(read_in_tree(""))
})

# build_tree ------------------------------------------------------------------#
test_that("Check that build_tree() works on valid vcf input", {
  test_dna <- load_vcf_file(prewas::vcf)
  test_tree <- build_tree(test_dna)
  expect_equal(class(test_tree), "phylo")
  expect_equal(ape::Ntip(test_tree), 14)
})

test_that("Check that build_tree() gives error when given invalid input", {
  expect_error(build_tree(prewas::tree))
  expect_error(build_tree(10))
  expect_error(build_tree(""))
})

# root_tree -------------------------------------------------------------------#

# subset_tree_and_matrix ------------------------------------------------------#
