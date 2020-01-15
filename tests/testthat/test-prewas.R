# prewas ----------------------------------------------------------------------#
test_that("prewas() gives expected output when given valid input", {
  # Warning given because the outgroup will be dropped from the DNA matrix.
  test_results <- prewas::prewas(dna = prewas::vcf,
                                 tree = prewas::tree,
                                 outgroup = prewas::outgroup,
                                 gff = prewas::gff,
                                 anc = FALSE)
  expect_identical(prewas::results, test_results)
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
