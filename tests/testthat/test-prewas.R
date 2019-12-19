# prewas ----------------------------------------------------------------------#
test_that("Check that prewas() gives expected output when given valid input", {
  test_results <- prewas::prewas(dna = prewas::dna,
                                 tree = prewas::tree,
                                 outgroup = prewas::outgroup,
                                 gff = prewas::gff,
                                 anc = FALSE)
  recorded_results <- prewas::prewas_results
  expect_identical(recorded_results, test_results)
})

test_that("Check that prewas() gives an error when given incorrect VCF data", {
  expect_error(prewas::prewas(dna = NULL))
  expect_error(prewas::prewas(dna = "foo"))
  expect_error(prewas::prewas(dna =  matrix(0, 10, 10)))
})

test_that("Check that prewas() gives an error when given incorrect tree data", {
  expect_error(prewas::prewas(dna = prewas::dna, tree = 5))
  expect_error(prewas::prewas(dna = prewas::dna, tree = "foo"))
  expect_error(prewas::prewas(dna = prewas::dna, tree = matrix(0, 10, 10)))
})

test_that("Check that prewas() gives an error when given incorrect outgroup data", {
  expect_error(prewas::prewas(dna = prewas::dna,
                              tree = prewas:tree,
                              outgroup = 5))
  expect_error(prewas::prewas(dna = prewas::dna,
                              tree = prewas:tree,
                              outgroup = "foo"))
  expect_error(prewas::prewas(dna = prewas::dna,
                              tree = prewas:tree,
                              outgroup =  matrix(0, 10, 10)))
})

test_that("Check that prewas() gives an error when given incorrect gff data", {
  expect_error(prewas::prewas(dna = prewas::dna,
                              tree = prewas:tree,
                              outgroup = prewas:outgroup,
                              gff = 5))
  expect_error(prewas::prewas(dna = prewas::dna,
                              tree = prewas:tree,
                              outgroup = prewas:outgroup,
                              gff = "foo"))
  expect_error(prewas::prewas(dna = prewas::dna,
                              tree = prewas:tree,
                              outgroup = prewas:outgroup,
                              gff = matrix(0, 10, 10)))
})
