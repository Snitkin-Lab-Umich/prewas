# format_inputs ---------------------------------------------------------------#
# format_inputs ---------------------------------------------------------------#
test_that("format_inputs gives no errors when reading in test data included with prewas", {
  expect_silent(format_inputs(dna = prewas::vcf,
                              tree = prewas::tree,
                              outgroup = prewas::outgroup,
                              gff = prewas::gff))
})

test_that("format_inputs gives errors when given no data", {
  expect_error(format_inputs(dna = NULL))
})
