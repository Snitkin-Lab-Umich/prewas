# reference_alleles -----------------------------------------------------------#
# make_all_tree_edges_positive ------------------------------------------------#
test_that("make_all_tree_edges_positive haves as expected when given valid inputs", {
  # Create tree with zero and negative edge lengths
  set.seed(1)
  temp_tree <- ape::rcoal(n = 10)
  temp_tree$edge.length[1] <- -.0001
  temp_tree$edge.length[2] <- 0

  expect_warning(only_positive_edges <- make_all_tree_edges_positive(temp_tree))
  expect_true(sum(only_positive_edges$edge.length > 0) == ape::Nedge(only_positive_edges))

  # Create a tree with only positive edge lengths
  set.seed(1)
  temp_tree <- ape::rcoal(10)
  # Expect no warnings
  only_positive_edges <- make_all_tree_edges_positive(temp_tree)
  expect_true(
    sum(only_positive_edges$edge.length > 0) == ape::Nedge(only_positive_edges))
})

# get_major_alleles -----------------------------------------------------------#
test_that("get_major_alleles behaves as expected when given valid inputs", {
  allele_mat <- load_vcf_file(prewas::vcf)
  major_allele_df <- get_major_alleles(allele_mat)
  expect_equal(class(major_allele_df), "character")
  expect_equal(length(major_allele_df), nrow(allele_mat))
})

test_that("get_major_alleles gives error when given invalid input", {
  expect_error(get_major_alleles("foo"))
  expect_error(get_major_alleles(15))
})

# get_ancestral_alleles -------------------------------------------------------#
test_that("get_ancestral_alleles behaves as expected when given valid inputs", {
  temp_tree <- prewas::tree
  temp_tree <- root_tree(temp_tree, "t1")
  temp_dna <- load_vcf_file(prewas::vcf)
  temp_dna <- temp_dna[1:10, , drop = FALSE]
  expect_warning(subsetted_data <- subset_tree_and_matrix(temp_tree, temp_dna))
  temp_tree <- subsetted_data$tree
  temp_dna <- keep_only_variant_sites(subsetted_data$mat)
  temp_ar_results <- get_ancestral_alleles(temp_tree, temp_dna)

  expect_equal(length(temp_ar_results), 2)
  expect_identical(names(temp_ar_results), c("ar_results", "tree"))
  expect_true(ape::is.rooted(temp_ar_results$tree))
  expect_equal(class(temp_ar_results$ar_results), "data.frame")
  expect_equal(class(temp_ar_results$tree), "phylo")
  expect_equal(nrow(temp_ar_results$ar_results), nrow(temp_dna))
  expect_equal(ncol(temp_ar_results$ar_results), 2)
})

test_that("get_ancestral_alleles gives error when given invalid input", {
  # DNA matrix not in correct format
  expect_error(get_ancestral_alleles(prewas::tree, prewas::vcf))

  # Unrooted tree
  temp_dna <- load_vcf_file(prewas::vcf)
  temp_dna <- temp_dna[1:10, , drop = FALSE]
  expect_error(get_ancestral_alleles(prewas::tree, temp_dna))

})

# remove_unknown_alleles ------------------------------------------------------#
test_that("remove_unknown_alleles correctly removes Ns when given valid input", {
  temp_tree <- prewas::tree
  temp_tree <- root_tree(temp_tree, "t1")
  temp_dna <- load_vcf_file(prewas::vcf)
  temp_dna <- temp_dna[1:10, , drop = FALSE]
  expect_warning(subsetted_data <- subset_tree_and_matrix(temp_tree, temp_dna))
  temp_tree <- subsetted_data$tree
  temp_dna <- keep_only_variant_sites(subsetted_data$mat)
  temp_ar_results <- get_ancestral_alleles(temp_tree, temp_dna)
  temp_ar_results$ar_results$ancestral_allele <- as.character(temp_ar_results$ar_results$ancestral_allele)
  temp_ar_results$ar_results$ancestral_allele[1] <- "N"
  temp_ar_results$ar_results$ancestral_allele[2] <- "-"
  temp_ar_results$ar_results$ancestral_allele <- as.factor(temp_ar_results$ar_results$ancestral_allele)


  expect_warning(temp_remove <-
    remove_unknown_alleles(temp_dna,
                           temp_ar_results$ar_results$ancestral_allele,
                           temp_ar_results$ar_results))
})

test_that("remove_unknown_alleles gives error when given invalid input", {
  expect_error(remove_unknown_alleles("foo", "foo", "foo"))
})

# make_binary_matrix ----------------------------------------------------------#
test_that("make_binary_matrix performs as expected when given valid input", {
  temp_tree <- prewas::tree
  temp_tree <- root_tree(temp_tree, "t1")
  temp_dna <- load_vcf_file(prewas::vcf)
  temp_dna <- temp_dna[1:10, , drop = FALSE]
  expect_warning(subsetted_data <- subset_tree_and_matrix(temp_tree, temp_dna))
  temp_tree <- subsetted_data$tree
  temp_dna <- keep_only_variant_sites(subsetted_data$mat)
  temp_ar_results <- get_ancestral_alleles(temp_tree, temp_dna)

  split <-
    split_multi_to_biallelic_snps(mat = temp_dna,
                                  ar_results = temp_ar_results$ar_results)

  allele_mat_split <- split$mat_split
  allele_results_split <- split$ar_results_split
  split_rows_flag <- split$split_rows_flag

  alleles <- allele_results_split$ancestral_allele
  names(alleles) <- rownames(allele_results_split)

  temp_bin <- make_binary_matrix(allele_mat_split, alleles)

  expect_silent(check_if_binary_matrix(temp_bin))
})

test_that("make_binary_matrix gives error when given invalid input", {
  expect_error(make_binary_matrix("foo", "foo"))
})
