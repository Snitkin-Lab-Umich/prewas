# preprocess_tree_and_vcf -----------------------------------------------------#
# read_in_tree ----------------------------------------------------------------#
test_that("read_in_tree() works for valid tree input", {
  test_tree <- read_in_tree(prewas::tree)
  expect_equal(ape::Ntip(test_tree), ape::Ntip(prewas::tree))
})

test_that("read_in_tree() gives error for invalid inputs", {
  expect_error(read_in_tree(10))
  expect_error(read_in_tree("tests/testthat/test-preprocess_tree_and_fasta.R"))
  expect_error(read_in_tree(""))
})

# build_tree ------------------------------------------------------------------#
test_that("build_tree() works on valid vcf input", {
  test_dna <- load_vcf_file(prewas::vcf)
  test_tree <- build_tree(test_dna)
  expect_equal(class(test_tree), "phylo")
  expect_equal(ape::Ntip(test_tree), 14)
})

test_that("build_tree() gives error when given invalid input", {
  expect_error(build_tree(prewas::tree))
  expect_error(build_tree(10))
  expect_error(build_tree(""))
})

# root_tree -------------------------------------------------------------------#
test_that("root_tree roots tree and drops outgroup when given valid inputs", {

  # Given an unrooted tree with an outgroup
  rooted_tree <- root_tree(prewas::tree, "t1")
  expect_true(ape::is.rooted(rooted_tree))
  expect_true(sum(grepl("t1$", rooted_tree$tip.label)) == 0)
  expect_equal(ape::Ntip(rooted_tree), ape::Ntip(prewas::tree) - 1)

  # Given a rooted tree and an outgroup
  tree <- ape::root(prewas::tree, "t1")
  tree <- ape::drop.tip(tree, "t1")
  rerooted_tree <- root_tree(tree, "t2")
  expect_true(ape::is.rooted(rerooted_tree))
  expect_true(sum(grepl("t2$", rerooted_tree$tip.label)) == 0)
  expect_equal(ape::Ntip(rerooted_tree), ape::Ntip(prewas::tree) - 2)

  # Given an unrooted tree and no outgroup
  rooted_tree <- root_tree(prewas::tree)
  expect_true(ape::is.rooted(rooted_tree))
  expect_equal(ape::Ntip(rooted_tree), ape::Ntip(prewas::tree))
})

test_that("root_tree() gives error for invalid inputs", {
  # No tree input
  expect_error(root_tree("foo", "bar"))
  expect_error(root_tree("tests/testthat/test-preprocess_tree_and_fasta.R"))

  # Wrong outgroup
  expect_error(root_tree(prewas::tree, "foo"))
  expect_error(root_tree(prewas::tree, 7))
})


# subset_tree_and_matrix ------------------------------------------------------#
test_that("subset_tree_and_matrix() works as expected on valid inputs", {
  # Tree tips and colnames match
  dna <- load_vcf_file(prewas::vcf)
  tree <- prewas::tree
  test_results <- subset_tree_and_matrix(tree, dna)
  expect_equal(length(test_results), 2)
  expect_identical(c("tree", "mat"), names(test_results))
  expect_equal(ape::Ntip(prewas::tree), ape::Ntip(test_results$tree))
  expect_identical(test_results$mat, dna)
  expect_identical(colnames(test_results$mat), test_results$tree$tip.label)

  # Same as above, but start with column names jumbled in dna object to check
  # that the matrix gets reordered
  temp_dna <- dna
  colnames(temp_dna) <- paste0("t", 1:14)
  test_results <- subset_tree_and_matrix(tree, temp_dna)
  expect_equal(length(test_results), 2)
  expect_identical(c("tree", "mat"), names(test_results))
  expect_equal(ape::Ntip(prewas::tree), ape::Ntip(test_results$tree))
  expect_identical(colnames(test_results$mat), test_results$tree$tip.label)

  # Some tree tips not found in matrix
  temp_dna <- dna[, 1:5, drop = FALSE]
  expect_warning(temp_results <- subset_tree_and_matrix(tree, temp_dna))
  expect_equal(length(temp_results), 2)
  expect_identical(c("tree", "mat"), names(temp_results))
  expect_equal(5, ape::Ntip(temp_results$tree))
  expect_equal(nrow(temp_results$mat), nrow(dna))
  expect_equal(ncol(temp_results$mat), 5)

  # Some matrix columns not found in tree
  temp_tree <- ape::drop.tip(tree, c("t1", "t2", "t3", "t4"))
  expect_warning(temp_results <- subset_tree_and_matrix(temp_tree, dna))
  expect_equal(length(temp_results), 2)
  expect_identical(c("tree", "mat"), names(temp_results))
  expect_equal(10, ape::Ntip(temp_results$tree))
  expect_equal(nrow(temp_results$mat), nrow(dna))
  expect_equal(ncol(temp_results$mat), 10)

  # Only 1 sample shared between tree and matrix
  temp_dna <- dna[, 1, drop = FALSE]
  expect_error(temp_results <- subset_tree_and_matrix(tree, temp_dna))

  # No samples shared between tree and matrix
  colnames(dna) <- letters[1:ncol(dna)]
  expect_error(subset_tree_and_matrix(tree, dna))

})

test_that("subset_tree_and_matrix() gives error when given invalid inputs", {
  # No tree
  expect_error(subset_tree_and_matrix("foo", prewas::vcf))
  expect_error(subset_tree_and_matrix(0, prewas::vcf))
  expect_error(subset_tree_and_matrix("tests/testthat/test-preprocess_tree_and_fasta.R", prewas::vcf))

  # No matrix
  expect_error(subset_tree_and_matrix(prewas::tree, "foo"))
  expect_error(subset_tree_and_matrix(prewas::tree, 0))
})
