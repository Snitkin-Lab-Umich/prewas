# prewas ----------------------------------------------------------------------#
test_that("prewas() gives expected output when given valid input", {
  # Warning given because the outgroup will be dropped from the DNA matrix.
  test_results <- prewas::prewas(dna = prewas::vcf,
                                 tree = prewas::tree,
                                 outgroup = prewas::outgroup,
                                 gff = prewas::gff,
                                 anc = FALSE)
  expect_identical(prewas::results$ar_results, test_results$ar_results)
  expect_identical(prewas::results$dup, test_results$dup)
  expect_identical(prewas::results$allele_mat, test_results$allele_mat)
  expect_identical(prewas::results$gene_mat, test_results$gene_mat)
  expect_identical(prewas::results$tree, test_results$tree)
  expect_identical(colSums(prewas::results$bin_mat), colSums(test_results$bin_mat))
  # Don't test that prewas::results$bin_mat and test_results$bin_mat are
  #   identical because rows can get flipped. Example:
  # > test_results$bin_mat[grepl("973", row.names(test_results$bin_mat)), ]
  # t7 t3 t10 t6 t8 t2 t5 t9 t11 t1 t13 t4 t14 t12
  # 973|gene_98    0  0   0  0  1  1  0  0   0  0   0  0   0   0
  # 973.1|gene_98  1  1   1  1  0  0  0  0   0  0   0  0   0   0
  # > prewas::results$bin_mat[grepl("973", row.names(prewas::results$bin_mat)),]
  # t7 t3 t10 t6 t8 t2 t5 t9 t11 t1 t13 t4 t14 t12
  # 973|gene_98    1  1   1  1  0  0  0  0   0  0   0  0   0   0
  # 973.1|gene_98  0  0   0  0  1  1  0  0   0  0   0  0   0   0
})

test_that("prewas() gives an error when given either a non-vcfR object or a non-file-path as the VCF data", {
  expect_error(prewas::prewas(dna = NULL))
  expect_error(prewas::prewas(dna = "foo"))
  expect_error(prewas::prewas(dna =  matrix(0, 10, 10)))
})

test_that("prewas() gives an error when given a non-tree object", {
  expect_error(prewas::prewas(dna = prewas::vcf, tree = 5))
  expect_error(prewas::prewas(dna = prewas::vcf, tree = "foo"))
  expect_error(prewas::prewas(dna = prewas::vcf, tree = matrix(0, 10, 10)))
})

test_that("prewas() gives an error when given an non-character outgroup object", {
  expect_error(prewas::prewas(dna = prewas::vcf,
                              tree = prewas::tree,
                              outgroup = 5))
  expect_error(prewas::prewas(dna = prewas::vcf,
                              tree = prewas::tree,
                              outgroup =  matrix(0, 10, 10)))
})

test_that("prewas() gives a warning when given an outgroup not found in the tree", {
  expect_warning(prewas::prewas(dna = prewas::vcf,
                                tree = prewas::tree,
                                outgroup = "foo"))
})

test_that("prewas() gives an error when given incorrect gff data", {
  expect_error(prewas::prewas(dna = prewas::vcf,
                              tree = prewas::tree,
                              outgroup = prewas::outgroup,
                              gff = 5))
  expect_error(prewas::prewas(dna = prewas::vcf,
                              tree = prewas::tree,
                              outgroup = prewas::outgroup,
                              gff = "foo"))
  expect_error(prewas::prewas(dna = prewas::vcf,
                              tree = prewas::tree,
                              outgroup = prewas::outgroup,
                              gff = matrix(0, 10, 10)))
})
