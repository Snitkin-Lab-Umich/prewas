# collapse_snps_by_gene -------------------------------------------------------#
# get_gene_names --------------------------------------------------------------#
test_that("get_gene_names() extracts correct gene names from a valid bin_mat", {
  temp_bin_mat <- matrix(c(0, 1, 1), nrow = 10, ncol = 3)
  row.names(temp_bin_mat) <-
    paste0(1:10, "|gene_", c(rep("1", 3), rep("2", 5), rep("3", 2)))
  colnames(temp_bin_mat) <- c("t1", "t2", "t3")
  gene_names <- get_gene_names(temp_bin_mat)
  expect_true(methods::is(gene_names, "character"))
  expect_equal(length(gene_names), nrow(temp_bin_mat))
})

test_that("get_gene_names() gives error when not given a binary matrix", {
  expect_error(get_gene_names(0))
  expect_error(get_gene_names(""))
  expect_error(get_gene_names(as.data.frame(matrix(0, 10, 10))))
  expect_error(get_gene_names(matrix(1:10, 2, 5)))
})

# collapse_snps_into_genes ----------------------------------------------------#
test_that("collpase_snps_into_genes() behaves as expected when given valid input", {
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

test_that("collapse_snps_into_genes() give error when given invalid matrix", {
  # Non-binary matrix
  temp_bin_mat <- matrix(c(0, 0, 1), nrow = 10, ncol = 3)
  row.names(temp_bin_mat) <-
    paste0(1:10, "|gene_", c(rep("1", 3), rep("2", 5), rep("3", 2)))
  colnames(temp_bin_mat) <- c("t1", "t2", "t3")
  gene_names <- get_gene_names(temp_bin_mat)
  temp_bin_mat[temp_bin_mat == 1] <- 2
  expect_error(collapse_snps_into_genes(temp_bin_mat, gene_names))
})

test_that("collapse_snps_into_genes() behaves as expected when given valid allele grouping", {
  # When multiallelic sites
  temp_bin_mat <- matrix(c(0, 0, 1), nrow = 10, ncol = 3)
  row.names(temp_bin_mat) <-
    paste0(1:10, "|gene_", c(rep("1", 3), rep("2", 5), rep("3", 2)))
  colnames(temp_bin_mat) <- c("t1", "t2", "t3")
  row.names(temp_bin_mat)[10] <- "9.1|gene_3"
  allele_names <- get_allele_names(temp_bin_mat)

  temp_new_bin_mat <- collapse_snps_into_genes(temp_bin_mat, allele_names)
  expect_equal(unname(temp_new_bin_mat[9, , drop = TRUE]), c(1, 0, 1))

  # When there is no collapsing possible
  temp_bin_mat <- matrix(c(0, 0, 1), nrow = 10, ncol = 3)
  row.names(temp_bin_mat) <-
    paste0(1:10, "|gene_", c(rep("1", 3), rep("2", 5), rep("3", 2)))
  colnames(temp_bin_mat) <- c("t1", "t2", "t3")
  allele_names <- get_allele_names(temp_bin_mat)

  temp_new_bin_mat <- collapse_snps_into_genes(temp_bin_mat, allele_names)
  expect_equal(dim(temp_new_bin_mat), dim(temp_bin_mat))
})

  # collapse_snps_into_genes_by_impact ----------------------------------------#
  test_that("collapse_snps_into_genes_by_impact() behaves as expected when given valid input", {
    # When collapsing is possible
    temp_bin_mat <- matrix(c(0, 0, 1), nrow = 10, ncol = 3)
    row.names(temp_bin_mat) <-
      paste0(1:10, "|gene_", c(rep("1", 3), rep("2", 5), rep("3", 2)))
    colnames(temp_bin_mat) <- c("t1", "t2", "t3")
    gene_names <- get_gene_names(temp_bin_mat)
    set.seed(5)
    temp_predicted_impact <- sample(c("MODIFIER", "MODERATE", "HIGH", "LOW"),
                                    length(gene_names),
                                    replace = TRUE)
    snpeff_grouping <- NULL
    temp_gene_mat_list <- collapse_snps_into_genes_by_impact(temp_bin_mat,
                                                        gene_names,
                                                        temp_predicted_impact,
                                                        snpeff_grouping)
    expect_true(length(temp_gene_mat_list) == 6)
    expect_true(length(unique(gene_names)) ==
                  nrow(temp_gene_mat_list$gene_mat_all))
    expect_true(sum(temp_predicted_impact == "HIGH") ==
                  sum(temp_gene_mat_list$gene_mat_high))
    expect_true(sum(temp_predicted_impact == "MODERATE") ==
                  sum(temp_gene_mat_list$gene_mat_moderate))
    expect_true(sum(temp_predicted_impact == "LOW") ==
                  sum(temp_gene_mat_list$gene_mat_low))
    expect_true(sum(temp_predicted_impact == "MODIFIER") ==
                  sum(temp_gene_mat_list$gene_mat_modifier))
    expect_true(sum(temp_gene_mat_list$gene_mat_custom) == 0)

    expect_equal(
      unname(temp_gene_mat_list$gene_mat_all[1, , drop = TRUE]), c(1, 1, 1))
    expect_equal(
      unname(temp_gene_mat_list$gene_mat_all[2, , drop = TRUE]), c(1, 1, 1))

    expect_equal(unname(
      temp_gene_mat_list$gene_mat_moderate[1, , drop = TRUE]), c(0, 0, 1))
    expect_equal(unname(
      temp_gene_mat_list$gene_mat_moderate[2, , drop = TRUE]), c(1, 0, 0))
    expect_equal(unname(
      temp_gene_mat_list$gene_mat_high[2, , drop = TRUE]), c(0, 1, 1))

    # if snpeff_grouping is not NULL
    snpeff_grouping <- c("HIGH", "MODERATE")
    temp_gene_mat_list <-
      collapse_snps_into_genes_by_impact(temp_bin_mat,
                                         gene_names,
                                         temp_predicted_impact,
                                         snpeff_grouping)

    expect_true(length(temp_gene_mat_list) == 6)
    expect_true(length(unique(gene_names)) ==
                  nrow(temp_gene_mat_list$gene_mat_all))
    expect_true(sum(temp_predicted_impact == "HIGH") ==
                  sum(temp_gene_mat_list$gene_mat_high))
    expect_true(sum(temp_predicted_impact == "MODERATE") ==
                  sum(temp_gene_mat_list$gene_mat_moderate))
    expect_true(sum(temp_predicted_impact == "LOW") ==
                  sum(temp_gene_mat_list$gene_mat_low))
    expect_true(sum(temp_predicted_impact == "MODIFIER") ==
                  sum(temp_gene_mat_list$gene_mat_modifier))
    expect_true(sum(temp_gene_mat_list$gene_mat_custom) ==
                  (sum(temp_gene_mat_list$gene_mat_moderate) +
                     sum(temp_gene_mat_list$gene_mat_high)))

    expect_equal(unname(
      temp_gene_mat_list$gene_mat_moderate[1, ,drop = TRUE]), c(0, 0, 1))
    expect_equal(unname(
      temp_gene_mat_list$gene_mat_moderate[2, ,drop = TRUE]), c(1, 0, 0))
    expect_equal(unname(
      temp_gene_mat_list$gene_mat_high[2, ,drop = TRUE]), c(0, 1, 1))

    expect_equal(unname(
      temp_gene_mat_list$gene_mat_custom[1, ,drop = TRUE]), c(0, 1, 1))
    expect_equal(unname(
      temp_gene_mat_list$gene_mat_custom[3, ,drop = TRUE]), c(1, 0, 1))


    # When there is no collapsing possible
    temp_bin_mat <- matrix(c(0, 0, 1), nrow = 10, ncol = 3)
    row.names(temp_bin_mat) <-
      paste0(1:10, "|gene_", 1:10)
    colnames(temp_bin_mat) <- c("t1", "t2", "t3")
    gene_names <- get_gene_names(temp_bin_mat)

    temp_gene_mat_list <- collapse_snps_into_genes_by_impact(temp_bin_mat,
                                                        gene_names,
                                                        temp_predicted_impact,
                                                        snpeff_grouping)
    expect_equal(dim(temp_gene_mat_list$gene_mat_all), dim(temp_bin_mat))
  })

  test_that("collapse_snps_into_genes_by_impact() give error when given invalid matrix", {
    # Non-binary matrix
    temp_bin_mat <- matrix(c(0, 0, 1), nrow = 10, ncol = 3)
    row.names(temp_bin_mat) <-
      paste0(1:10, "|gene_", c(rep("1", 3), rep("2", 5), rep("3", 2)))
    colnames(temp_bin_mat) <- c("t1", "t2", "t3")
    gene_names <- get_gene_names(temp_bin_mat)
    temp_bin_mat[temp_bin_mat == 1] <- 2

    set.seed(5)
    temp_predicted_impact <- sample(c("MODIFIER", "MODERATE", "HIGH", "LOW"),
                                    length(gene_names),
                                    replace = TRUE)
    snpeff_grouping <- NULL
    expect_error(
      collapse_snps_into_genes_by_impact(temp_bin_mat,
                                         gene_names,
                                         temp_predicted_impact,
                                         snpeff_grouping))
  })

  test_that("collapse_snps_into_genes_by_impact() behaves as expected when given valid allele grouping",{
    # When multiallelic sites

    temp_bin_mat <- matrix(c(0, 0, 1), nrow = 10, ncol = 3)
    row.names(temp_bin_mat) <-
      paste0(1:10, "|gene_", c(rep("1", 3), rep("2", 5), rep("3", 2)))
    colnames(temp_bin_mat) <- c("t1", "t2", "t3")
    row.names(temp_bin_mat)[10] <- "9.1|gene_3"
    allele_names <- get_allele_names(temp_bin_mat)

    set.seed(5)
    temp_predicted_impact <- sample(c("MODIFIER", "MODERATE", "HIGH", "LOW"),
                                    length(allele_names),
                                    replace = TRUE)
    snpeff_grouping <- NULL

    temp_new_bin_mat_list <-
      collapse_snps_into_genes_by_impact(temp_bin_mat,
                                         allele_names,
                                         temp_predicted_impact,
                                         snpeff_grouping)
    expect_equal(
      unname(temp_new_bin_mat_list$gene_mat_all[9, , drop = TRUE]), c(1, 0, 1))

    # When there is no collapsing possible
    temp_bin_mat <- matrix(c(0, 0, 1), nrow = 10, ncol = 3)
    row.names(temp_bin_mat) <-
      paste0(1:10, "|gene_", c(rep("1", 3), rep("2", 5), rep("3", 2)))
    colnames(temp_bin_mat) <- c("t1", "t2", "t3")
    allele_names <- get_allele_names(temp_bin_mat)

    temp_new_bin_mat_list <-
      collapse_snps_into_genes_by_impact(temp_bin_mat,
                                         allele_names,
                                         temp_predicted_impact,
                                         snpeff_grouping)
    expect_equal(dim(temp_new_bin_mat_list$gene_mat_all), dim(temp_bin_mat))
  })

  # get_gene_mat_by_impact ----------------------------------------------------#
  test_that("get_gene_mat_by_impact() behaves as expected when given valid input",{

    temp_bin_mat <- matrix(c(0, 0, 1), nrow = 10, ncol = 3)
    row.names(temp_bin_mat) <-
      paste0(1:10, "|gene_", c(rep("1", 3), rep("2", 5), rep("3", 2)))
    colnames(temp_bin_mat) <- c("t1", "t2", "t3")
    gene_names <- get_gene_names(temp_bin_mat)

    unique_gene_names <- unique(gene_names)
    num_unique_genes <- length(unique_gene_names)

    set.seed(5)
    temp_predicted_impact <- sample(c("MODIFIER", "MODERATE", "HIGH", "LOW"),
                                    length(gene_names),
                                    replace = TRUE)

    gene_mat <- get_gene_mat_by_impact(num_unique_genes,
                         unique_gene_names,
                         gene_names,
                         temp_bin_mat,
                         temp_predicted_impact,
                         "ALL")

    expect_equal(unname(gene_mat[1, , drop = TRUE]), c(1,1,1))

    gene_mat <- get_gene_mat_by_impact(num_unique_genes,
                                       unique_gene_names,
                                       gene_names,
                                       temp_bin_mat,
                                       temp_predicted_impact,
                                       "MODERATE")

    expect_equal(unname(gene_mat[2, , drop = TRUE]), c(1,0,0))
  })
