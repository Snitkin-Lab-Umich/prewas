# format_inputs ---------------------------------------------------------------#
# format_inputs ---------------------------------------------------------------#
test_that("format_inputs gives expected results prewas test data", {
  vcf <- prewas::vcf
  tree <- prewas::tree
  outgroup <- prewas::outgroup
  gff <- prewas::gff
  results <- format_inputs(dna = vcf,
                           tree = tree,
                           outgroup = outgroup,
                           gff = gff)
  expect_equal(length(results), 4)
  expect_true(methods::is(results$dna, "matrix"))
  expect_true(methods::is(results$tree, "phylo"))
  expect_true(methods::is(results$outgroup, "character"))
  expect_true(methods::is(results$gff, "matrix"))


})

test_that("format_inputs gives errors when given no data", {
  expect_error(format_inputs(dna = NULL))
})
