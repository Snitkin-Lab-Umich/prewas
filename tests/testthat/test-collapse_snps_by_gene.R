collapse_snps_by_gene -------------------------------------------------------#
get_gene_names --------------------------------------------------------------#
test_that("Check that get_gene_names() extracts correct gene names from a valid bin_mat", {
  temp_bin_mat <- matrix(c(0, 1, 1), nrow = 10, ncol = 3)
  row.names(temp_bin_mat) <-
    paste0(1:10, "|gene_", c(rep("1", 3), rep("2", 5), rep("3", 2)))
  colnames(temp_bin_mat) <- c("t1", "t2", "t3")
  gene_names <- get_gene_names(temp_bin_mat)
  expect_equal(class(gene_names), "character")
  expect_equal(length(gene_names), nrow(temp_bin_mat))
})

test_that("Check that get_gene_names() gives error when not given a binary matrix", {
  expect_error(get_gene_names(0))
  expect_error(get_gene_names(""))
  expect_error(get_gene_names(as.data.frame(matrix(0, 10, 10))))
  expect_error(get_gene_names(matrix(1:10, 2, 5)))
})

# collapse_snps_into_genes ----------------------------------------------------#
test_that("Check that collpase_snps_into_genes() behaves as expected when given valid input", {
  # When collapsing is possible
  temp_bin_mat <- matrix(c(0, 0, 1), nrow = 10, ncol = 3)
  row.names(temp_bin_mat) <-
    paste0(1:10, "|gene_", c(rep("1", 3), rep("2", 5), rep("3", 2)))
  colnames(temp_bin_mat) <- c("t1", "t2", "t3")
  gene_names <- get_gene_names(temp_bin_mat)

  temp_gene_mat <- collapse_snps_into_genes(temp_bin_mat, gene_names)
  expect_equal(unname(temp_gene_mat[3, , drop = TRUE]), c(1, 0, 1))
  expect_equal(unname(temp_gene_mat[2, , drop = TRUE]), c(1, 1, 1))

  # When there is no collapsing possible
  temp_bin_mat <- matrix(c(0, 0, 1), nrow = 10, ncol = 3)
  row.names(temp_bin_mat) <-
    paste0(1:10, "|gene_", 1:10)
  colnames(temp_bin_mat) <- c("t1", "t2", "t3")
  gene_names <- get_gene_names(temp_bin_mat)

  temp_gene_mat <- collapse_snps_into_genes(temp_bin_mat, gene_names)
  expect_equal(dim(temp_gene_mat), dim(temp_bin_mat))
})

test_that("Check that collapse_snps_into_genes() give error when given invalid matrix", {
  # Non-binary matrix
  temp_bin_mat <- matrix(c(0, 0, 1), nrow = 10, ncol = 3)
  row.names(temp_bin_mat) <-
    paste0(1:10, "|gene_", c(rep("1", 3), rep("2", 5), rep("3", 2)))
  colnames(temp_bin_mat) <- c("t1", "t2", "t3")
  gene_names <- get_gene_names(temp_bin_mat)
  temp_bin_mat[temp_bin_mat == 1] <- 2
  expect_error(collapse_snps_into_genes(temp_bin_mat, gene_names))
})
