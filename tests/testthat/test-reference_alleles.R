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
  major_allele_df <- get_major_alleles(allele_mat$vcf_geno_mat)
  expect_true(methods::is(major_allele_df, "character"))
  expect_equal(length(major_allele_df), nrow(allele_mat$vcf_geno_mat))
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
  temp_dna_mat <- temp_dna$vcf_geno_mat
  temp_dna_mat <- temp_dna_mat[1:10, , drop = FALSE]
  expect_warning(subsetted_data <- subset_tree_and_matrix(temp_tree, temp_dna_mat))
  temp_tree <- subsetted_data$tree
  temp_dna_mat <- keep_only_variant_sites(subsetted_data$mat,
                                          o_ref = temp_dna$vcf_ref_allele[1:10],
                                          o_alt = temp_dna$vcf_alt_allele[1:10],
                                          snpeff = temp_dna$snpeff_pred)
  temp_ar_results <- get_ancestral_alleles(temp_tree, temp_dna_mat$variant_only_dna_mat)

  expect_equal(length(temp_ar_results), 2)
  expect_identical(names(temp_ar_results), c("ar_results", "tree"))
  expect_true(ape::is.rooted(temp_ar_results$tree))
  expect_true(methods::is(temp_ar_results$ar_results, "data.frame"))
  expect_true(methods::is(temp_ar_results$tree, "phylo"))
  expect_equal(nrow(temp_ar_results$ar_results), nrow(temp_dna_mat$variant_only_dna_mat))
  expect_equal(ncol(temp_ar_results$ar_results), 2)
})

test_that("get_ancestral_alleles gives error when given invalid input", {
  # DNA matrix not in correct format
  expect_error(get_ancestral_alleles(prewas::tree, prewas::vcf))

  # Unrooted tree
  temp_dna <- load_vcf_file(prewas::vcf)
  temp_dna <- temp_dna$vcf_geno_mat[1:10, , drop = FALSE]
  expect_error(get_ancestral_alleles(prewas::tree, temp_dna))

})

# remove_unknown_alleles ------------------------------------------------------#
test_that("remove_unknown_alleles correctly removes Ns when given valid input", {
  temp_tree <- prewas::tree
  temp_tree <- root_tree(temp_tree, "t1")
  temp_dna <- load_vcf_file(prewas::vcf)
  temp_dna_mat <- temp_dna$vcf_geno_mat[1:10, , drop = FALSE]
  expect_warning(subsetted_data <- subset_tree_and_matrix(temp_tree,
                                                          temp_dna_mat))
  temp_tree <- subsetted_data$tree
  temp_dna_list <- keep_only_variant_sites(dna_mat = subsetted_data$mat,
                                           o_ref = temp_dna$vcf_ref_allele[1:10],
                                           o_alt = temp_dna$vcf_alt_allele[1:10],
                                           snpeff = temp_dna$snpeff_pred[1:10])
  temp_dna_mat <- temp_dna_list$variant_only_dna_mat
  temp_ar_results <- get_ancestral_alleles(temp_tree, temp_dna_mat)
  temp_ar_results$ar_results$ancestral_allele <-
    as.character(temp_ar_results$ar_results$ancestral_allele)
  temp_ar_results$ar_results$ancestral_allele[1] <- "N"
  temp_ar_results$ar_results$ancestral_allele[2] <- "-"
  temp_ar_results$ar_results$ancestral_allele <-
    as.factor(temp_ar_results$ar_results$ancestral_allele)


  # Data with Ns and dashes
  expect_true(sum(temp_ar_results$ar_results$ancestral_allele %in% c("N", "-")) == 2)

  expect_warning(temp_remove <-
    remove_unknown_alleles(allele_mat = temp_dna_mat,
                           alleles = temp_ar_results$ar_results$ancestral_allele,
                           ar_results = temp_ar_results$ar_results,
                           o_ref = temp_dna_list$o_ref_var_pos,
                           o_alt = temp_dna_list$o_alt_var_pos,
                           snpeff = temp_dna_list$snpeff_var_pos))

  # Now data no longer has Ns and dashes
  expect_true(sum(temp_remove$ar_results$ancestral_allele %in% c("N", "-")) == 0)
})

test_that("remove_unknown_alleles gives error when given invalid input", {
  expect_error(remove_unknown_alleles("foo", "foo", "foo", "foo", "foo", "foo"))
})

# make_binary_matrix ----------------------------------------------------------#
test_that("make_binary_matrix performs as expected when given valid input", {
  temp_tree <- prewas::tree
  temp_tree <- root_tree(temp_tree, "t1")
  temp_dna <- load_vcf_file(prewas::vcf)
  temp_dna_mat <- temp_dna$vcf_geno_mat
  temp_dna_mat <- temp_dna_mat[1:10, , drop = FALSE]
  expect_warning(subsetted_data <- subset_tree_and_matrix(temp_tree, temp_dna_mat))
  temp_tree <- subsetted_data$tree
  ref <- unname(subsetted_data$mat[, 1, drop = TRUE])
  alt <- c("G", "A", "A",
           "T", "A,G", "T,C",
           "T", "T", "T", "C")
  temp_dna_list <- keep_only_variant_sites(subsetted_data$mat,
                                           o_ref = ref,
                                           o_alt = alt,
                                           snpeff = NULL)
  temp_ar_results <-
    get_ancestral_alleles(tree = temp_tree,
                          mat = temp_dna_list$variant_only_dna_mat)

  split <-
    split_multi_to_biallelic_snps(mat = temp_dna_list$variant_only_dna_mat,
                                  ar_results = temp_ar_results$ar_results,
                                  o_ref = temp_dna_list$o_ref_var_pos,
                                  o_alt = temp_dna_list$o_alt_var_pos,
                                  snpeff = NULL)

  allele_mat_split <- split$mat_split
  allele_results_split <- split$ar_results_split

  alleles <- allele_results_split$ancestral_allele
  names(alleles) <- rownames(allele_results_split)

  temp_bin <- make_binary_matrix(allele_mat = allele_mat_split,
                                 n_ref = alleles,
                                 o_ref = split$o_ref_split,
                                 o_alt = split$o_alt_split)

  expect_silent(check_if_binary_matrix(temp_bin[[1]]))
})

test_that("make_binary_matrix gives error when given invalid input", {
  expect_error(make_binary_matrix("foo", "foo", "foo", "foo"))
})


# parse_snpeff ----------------------------------------------------------------#
test_that("parse_snpeff performs as expected when given valid input", {
  vcf_with_snpeff <- load_vcf_file(prewas::snpeff_vcf)
  temp_snpeff <- vcf_with_snpeff$snpeff_pred[1:10]
  temp_snpeff[[5]] <- temp_snpeff[[6]] <- temp_snpeff[[3]]
  temp_snpeff[[3]] <- temp_snpeff[[1]]
  temp_tree <- prewas::tree
  temp_tree <- root_tree(temp_tree, "t1")
  temp_dna <- load_vcf_file(prewas::vcf)
  temp_dna_mat <- temp_dna$vcf_geno_mat
  temp_dna_mat <- temp_dna_mat[1:10, , drop = FALSE]
  expect_warning(subsetted_data <- subset_tree_and_matrix(temp_tree,
                                                          temp_dna_mat))
  temp_tree <- subsetted_data$tree
  ref <- unname(subsetted_data$mat[, 1, drop = TRUE])
  alt <- c("G", "A", "A",
           "T", "A,G", "T,C",
           "T", "T", "T", "C")

  for (i in c(1:4, 7:10)) {
    temp_snpeff[[i]] <- gsub("^[A-Z]", alt[i], temp_snpeff[[i]])
  }

  for (i in c(5:6)) {
    temp_snpeff[[i]][[1]] <-
      gsub("^[A-Z]", substr(alt[i], 1, 1), temp_snpeff[[i]][[1]])
    temp_snpeff[[i]][[2]] <-
      gsub("^[A-Z]", substr(alt[i], 3, 3), temp_snpeff[[i]][[2]])
  }

  temp_dna_list <- keep_only_variant_sites(subsetted_data$mat,
                                           o_ref = ref,
                                           o_alt = alt,
                                           snpeff = temp_snpeff)
  temp_ar_results <-
    get_ancestral_alleles(tree = temp_tree,
                          mat = temp_dna_list$variant_only_dna_mat)

  split <-
    split_multi_to_biallelic_snps(mat = temp_dna_list$variant_only_dna_mat,
                                  ar_results = temp_ar_results$ar_results,
                                  o_ref = temp_dna_list$o_ref_var_pos,
                                  o_alt = temp_dna_list$o_alt_var_pos,
                                  snpeff = temp_dna_list$snpeff_var_pos)

  parse_output <- parse_snpeff(alt_allele = split$o_alt_split,
                               snpeff_split = split$snpeff_split)

  expected_impacts <- c(rep("MODIFIER", 8), "LOW", "MODERATE", "LOW", "HIGH")
  expected_genes <- c(rep("CHR_START-dnaA", 8),
                      rep("dnaA", 3),
                      "serS")
  # Tests
  expect_true(sum(parse_output$pred_impact != expected_impacts) == 0)
  expect_true(sum(parse_output$gene != expected_genes) == 0)
})

test_that("parse_snpeff gives error when given invalid input", {
  vcf_with_snpeff <- load_vcf_file(prewas::snpeff_vcf)
  temp_snpeff <- vcf_with_snpeff$snpeff_pred[1:10]
  temp_snpeff[[5]] <- temp_snpeff[[6]] <- temp_snpeff[[3]]
  temp_snpeff[[3]] <- temp_snpeff[[1]]
  temp_tree <- prewas::tree
  temp_tree <- root_tree(temp_tree, "t1")
  temp_dna <- load_vcf_file(prewas::vcf)
  temp_dna_mat <- temp_dna$vcf_geno_mat
  temp_dna_mat <- temp_dna_mat[1:10, , drop = FALSE]
  expect_warning(subsetted_data <- subset_tree_and_matrix(temp_tree,
                                                          temp_dna_mat))
  temp_tree <- subsetted_data$tree
  ref <- unname(subsetted_data$mat[, 1, drop = TRUE])
  alt <- c("G", "A", "A",
           "T", "A,G", "T,C",
           "T", "T", "T", "C")

  for (i in c(1:4, 7:10)) {
    temp_snpeff[[i]] <- gsub("^[A-Z]", alt[i], temp_snpeff[[i]])
  }

  for (i in c(5:6)) {
    temp_snpeff[[i]][[1]] <-
      gsub("^[A-Z]", substr(alt[i], 1, 1), temp_snpeff[[i]][[1]])
    temp_snpeff[[i]][[2]] <-
      gsub("^[A-Z]", substr(alt[i], 3, 3), temp_snpeff[[i]][[2]])
  }

  temp_dna_list <- keep_only_variant_sites(subsetted_data$mat,
                                           o_ref = ref,
                                           o_alt = alt,
                                           snpeff = temp_snpeff)
  temp_ar_results <-
    get_ancestral_alleles(tree = temp_tree,
                          mat = temp_dna_list$variant_only_dna_mat)

  split <-
    split_multi_to_biallelic_snps(mat = temp_dna_list$variant_only_dna_mat,
                                  ar_results = temp_ar_results$ar_results,
                                  o_ref = temp_dna_list$o_ref_var_pos,
                                  o_alt = temp_dna_list$o_alt_var_pos,
                                  snpeff = temp_dna_list$snpeff_var_pos)

  expect_error(parse_snpeff(alt_allele = NULL,
                            snpeff_split = split$snpeff_split))
  expect_error(parse_snpeff(alt_allele = split$o_alt_split,
                            snpeff_split = NULL))
  expect_error(parse_snpeff(alt_allele = split$o_alt_split,
                            snpeff_split = split$snpeff_split[[1]]))
  expect_error(parse_snpeff(alt_allele = split$o_alt_split[1:5],
                            snpeff_split = split$snpeff_split))
})
