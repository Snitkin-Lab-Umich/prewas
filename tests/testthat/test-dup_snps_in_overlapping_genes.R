# dup_snps_in_overlapping_genes -----------------------------------------------#
# dup_snps_in_overlapping_genes -----------------------------------------------#
test_that("dup_snps_in_overlapping_genes() works as expected when given valid inputs", {
  temp_gff <- read_gff(prewas::gff)
  temp_gff <- temp_gff[1:4, , drop = FALSE]

  temp_bin_mat <- matrix(c(0, 0, 1), nrow = 7, ncol = 3)
  colnames(temp_bin_mat) <- c("t1", "t2", "t3")
  row.names(temp_bin_mat) <- c(25:31)

  temp_dup_mat <- dup_snps_in_overlapping_genes(temp_bin_mat, temp_gff)

  expect_equal(nrow(temp_dup_mat), 22)
  expect_equal(ncol(temp_dup_mat), 3)
})

test_that("dup_snps_in_overlapping_genes() gives error when given invalid inputs", {
  expect_error(dup_snps_in_overlapping_genes("foo", matrix(0, 10, 10)))
  expect_error(dup_snps_in_overlapping_genes(matrix(0, 10, 10), "foo"))
  expect_error(dup_snps_in_overlapping_genes(matrix(0, 10, 10), matrix(0, 10, 10)))
})

# dup_snps_in_overlapping_genes_snpeff ----------------------------------------#
test_that("dup_snps_in_overlapping_genes_snpeff() works as expected when given valid inputs", {

  temp_bin_mat <- matrix(c(0, 0, 1), nrow = 7, ncol = 3)
  colnames(temp_bin_mat) <- c("t1", "t2", "t3")
  row.names(temp_bin_mat) <- c(25:31)

  temp_predicted_impact <- c('LOW', 'LOW|HIGH', 'MODERATE', 'HIGH|HIGH', 'MODIFIER', 'LOW', 'LOW|HIGH')
  temp_gene <- c('gene1', 'gene2|gene3', 'gene4', 'gene5|gene6', 'gene7', 'gene8', 'gene9|gene10')

  temp_dup_mat <- dup_snps_in_overlapping_genes_snpeff(temp_bin_mat, temp_predicted_impact, temp_gene)

  expect_equal(nrow(temp_dup_mat$bin_mat_dup), 10)
  expect_equal(length(temp_dup_mat$predicted_impact_split), 10)
  expect_equal(ncol(temp_dup_mat$bin_mat_dup), 3)
})

test_that("dup_snps_in_overlapping_genes() gives error when given invalid inputs", {
  temp_bin_mat <- matrix(c(0, 0, 1), nrow = 7, ncol = 3)
  colnames(temp_bin_mat) <- c("t1", "t2", "t3")
  row.names(temp_bin_mat) <- c(25:31)

  temp_predicted_impact <- c('LOW', 'LOW|HIGH', 'MODERATE', 'HIGH|HIGH', 'MODIFIER', 'LOW', 'LOW|HIGH')
  temp_gene <- c('gene1', 'gene2|gene3', 'gene4', 'gene5|gene6', 'gene7', 'gene8', 'gene9|gene10')

  temp_dup_mat <- dup_snps_in_overlapping_genes_snpeff(temp_bin_mat, temp_predicted_impact, temp_gene)

  expect_error(dup_snps_in_overlapping_genes_snpeff("foo", "foo", "foo"))
  expect_error(dup_snps_in_overlapping_genes_snpeff(matrix(0, 10, 10), "foo", "foo"))
})
