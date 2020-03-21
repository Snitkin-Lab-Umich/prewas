#' Check that function input is a numeric.
#'
#' Doesn't return anything. Gives an error message if input is not a number.
#'
#' @param num Number
#' @noRd
check_is_number <- function(num){
  if (!methods::is(num, "numeric")) {
    stop("Input must be a numeric")
  }
}

#' Tests if input object is a file.
#'
#' @param obj Any R object.
#' @noRd
#' @return is_file Logical
is_file <- function(obj){
  is_file <- FALSE
  if (methods::is(obj, "character")) {
    if (file.exists(obj)) {
      is_file <- TRUE
    }
  }
  return(is_file)
}

#' Tests if an input object has the specified class.
#'
#' @param obj Any R object.
#' @param current_class String. Name of the expected class of the R object.
#'
#' @return is_this_class: Logical.
#' @noRd
is_this_class <- function(obj, current_class){
  if (length(current_class) != 1) {
    stop("Current_class must have a length of 1")
  }
  if (!methods::is(current_class, "character")) {
    stop("Current_class is expected to be a string describing a class")
  }
  r_classes <- c("character",
                 "numeric",
                 "integer",
                 "logical",
                 "complex",
                 "phylo",
                 "DNAbin",
                 "phyDat",
                 "matrix",
                 "data.frame",
                 "factor",
                 "vcfR",
                 "dist",
                 "list")
  if (!(current_class %in% r_classes)) {
    stop("current_class is expected to be a R class")
  }

  is_this_class <- methods::is(obj, current_class)
  return(is_this_class)
}

#' Checks if an object is of the expected R class.
#'
#' Doesn't return anything. Gives error if the object is not of the expected R
#' class.
#'
#' @param obj Any R object.
#' @param current_class Character string. Name of R class
#' @noRd
check_is_this_class <- function(obj, current_class){
  class_log <- is_this_class(obj, current_class)
  if (class_log != TRUE) {
    stop(paste("Input must be a", current_class))
  }
}

#' Check that the input tree is actually a 'phylo' object.
#'
#' Doesn't return anything. Gives an error message if the object is not a
#' 'phylo' object.
#'
#' @param tree Phylogenetic tree.
#' @noRd
check_is_tree <- function(tree){
  if (!is_this_class(tree, "phylo")) {
    stop("Input requires either a path to a tree file or an ape phylo object")
  }
}

#' Check that the tree has a root.
#'
#' Doesn't return anything. Gives an error message if the object is not a
#' rooted tree.
#'
#' @param tree Phylogenetic tree.
#' @noRd
check_tree_is_rooted <- function(tree){
  check_is_tree(tree)
  if (!ape::is.rooted(tree)) {
    stop("Tree must be rooted.")
  }
}

#' Read in a GFF file given a file path
#'
#' @param gff Character or matrix. If character: path to GFF file.
#' @noRd
#' @return gff: matrix
read_gff <- function(gff){
  if (is_file(gff)) {
    gff <- utils::read.table(gff,
                             sep = "\t",
                             header = FALSE,
                             stringsAsFactors = FALSE,
                             quote = "",
                             comment.char = "#")
  }
  if (is_this_class(gff, "data.frame")) {
    gff <- as.matrix(gff)
  }

  check_is_this_class(gff, "matrix")

  if (ncol(gff) != 9) {
    stop("GFF file must have 9 columns")
  }
  if (sum(c(".", "+", "-") %in% unique(gff[, 7, drop = TRUE])) == 0) {
    stop("GFF file must only have strand information in column 7")
  }

  return(gff)
}

#' Keep only GFF features that are in CDS regions.
#'
#' @param gff Matrix. Rows = genomic features.
#' @noRd
#' @return gff: Matrix.
subset_gff <- function(gff){
  check_is_this_class(gff, "matrix")

  gff <- gff[gff[, 3] == "CDS", ] # subset gff on CDS
  if (nrow(gff) == 0) {
    stop(paste("GFF has no CDS regions. Remove the gff from prewas inputs and",
               "rerun or use a different GFF with CDS regions."))
  }
  return(gff)
}

#' Remove ID prefix from CDS feature descriptions (e.g. gene names).
#'
#' @param gff Data.frame. Rows = genomic features.
#' @noRd
#' @return gff: Data.frame.
clean_up_cds_name_from_gff <- function(gff){
  check_is_this_class(gff, "matrix")

  cds_name <- apply(gff, 1, function(row) {
    temp <- gsub("^ID=", "", row[9])
    temp <- gsub(";.*$", "", temp)
  })
  gff[, 9] <- cds_name
  return(gff)
}

#' Read in variant data from a VCF file or vcfR object prepared
#'
#' @description This function can now handle the multiVCF as prepared by
#'   bcftools / snpeff
#'
#' @param vcf Either character or vcfR. If character, it is a path to a VCF
#'   file.
#' @noRd
#' @return vcf_geno_mat: Matrix.
load_vcf_file <- function(vcf) {
  if (is_file(vcf)) {
    vcf <- vcfR::read.vcfR(file = vcf, convertNA = FALSE)
  }
  check_is_this_class(vcf, "vcfR")

  vcf_geno_mat <- vcf@gt
  vcf_geno_mat <- vcf_geno_mat[, colnames(vcf_geno_mat) != "FORMAT"]

  row.names(vcf_geno_mat) <- vcf@fix[, colnames(vcf@fix) == "POS", drop = TRUE]
  vcf_ref_allele <- vcf@fix[, colnames(vcf@fix) == "REF", drop = TRUE]
  vcf_alt_allele <- vcf@fix[, colnames(vcf@fix) == "ALT", drop = TRUE]

  geno_to_remove <- NULL
  for (i in 1:nrow(vcf_geno_mat)) {
    alt_alleles <- strsplit(vcf_alt_allele[i], split = ",")
    num_alt_alleles <- length(alt_alleles[[1]])

    # Simplify indel encoding from bcftools/snpeff for reference allele
    vcf_geno_mat[i, grepl(pattern = "^[.][/][.]", x = vcf_geno_mat[i, ])] <- "0"

    # Recode reference allele to character
    vcf_geno_mat[i, vcf_geno_mat[i, ] == "0"] <- vcf_ref_allele[i]

    for (j in 1:num_alt_alleles) {
      # Simplify indel encoding from bcftools/snpeff
      vcf_geno_mat[i,
                   grepl(pattern = paste0("^", j, "[/]", j),
                         x = vcf_geno_mat[i, ])] <- j

      # Recode numbers to allele characters
      vcf_geno_mat[i, vcf_geno_mat[i, ] == j] <- alt_alleles[[1]][j]
    }

    if (sum(
      vcf_geno_mat[i, grepl(pattern = "[/]", x = vcf_geno_mat[i, ])] > 0)) {
      geno_to_remove <- c(geno_to_remove, i)
    }
  }

  if (!is.null(geno_to_remove)) {
    warning(
      paste0("Removing ",
             length(geno_to_remove),
             " rows with heterozygous calls"))
    vcf_geno_mat <- vcf_geno_mat[-geno_to_remove, , drop = FALSE]
    vcf_alt_allele <- vcf_alt_allele[-geno_to_remove]
    vcf_ref_allele <- vcf_ref_allele[-geno_to_remove]
  }

  # Get SNPeff annotations, if they exist
  if (length(grep("SnpEff", vcf@meta)) > 0) {
    # Extract annotations
    snpeff_pred <- vcfR::extract_info_tidy(vcf, info_fields = "ANN")
    # Get unique annotations
    snpeff_pred <- sapply(snpeff_pred$ANN,
                          function(s){unique(unlist((strsplit(s, ","))))})
    snpeff_pred <- unname(snpeff_pred)

    if (!is.null(geno_to_remove)) {
      snpeff_pred <- snpeff_pred[-geno_to_remove]
    }
  } else {
    snpeff_pred <- NULL
  }
  return(list("vcf_geno_mat" = vcf_geno_mat,
              "vcf_ref_allele" = vcf_ref_allele,
              "vcf_alt_allele" = vcf_alt_allele,
              "snpeff_pred" = snpeff_pred))
}

#' Confirm that the tree and variant matrix contain exactly the same samples
#'
#' Doesn't return anything. Gives error if the two inputs do not match.
#'
#' @param tip_labels Character. Vector of tree$tip.labels.
#' @param colnames_mat Character. Vector of column names from variant matrix.
#' @noRd
check_setequal_tree_mat <- function(tip_labels, colnames_mat){
  if (!setequal(tip_labels, colnames_mat)) {
    stop("Tree and variant matrix sample names do not match.")
  }
}

#' Check that the matrix contains only 1s and/or 0s.
#'
#' Function gives error if matrix is not binary; otherwise, returns nothing.
#'
#' @param mat Matrix.
#' @noRd
check_if_binary_matrix <- function(mat) {
  check_is_this_class(mat, "matrix")

  if (sum(!(mat %in% c(0, 1))) > 0) {
    stop("Matrix should be only 1s and 0s")
  }
}

#' Check that the user supplied valid snpeff annotation groups
#'
#' @description Doesn't return anything. Stops prewas() from running if the
#'   user provided inputs are invalid.
#' @param snpeff_grouping NULL or vector of characters:
#'   c('HIGH','MODERATE', LOW', 'MODIFER')
#'
#' @noRd
check_snpeff_user_input <- function(snpeff_grouping){
  valid_annots <- c("MODERATE", "MODIFIER", "LOW", "HIGH")
  if (!is.null(snpeff_grouping)) {
    check_is_this_class(snpeff_grouping, "character")
    if (!sum(snpeff_grouping %in% valid_annots) > 0) {
      stop("If a user supplies asnpeff_grouping input, create a vector combination of 'HIGH','MODERATE', LOW', and/or 'MODIFER'. Must write the impact combinations in all caps (e.g. c('HIGH', 'MODERATE')).")
    }
  }
}

#' Check that the number of annotations match the number of genes at a position
#'
#' @description Doesn't return anything. Stops prewas() from running if the
#'   user provided inputs are invalid.
#' @param predicted_impact Character vector of snpeff predicted impact. Pipe
#' delimited when overlapping genes.
#' @param gene Character vector of genes/locus tag. Pipe
#' delimited when overlapping genes.
#' @noRd
check_num_overlap_genes_match_num_impact <- function(predicted_impact, gene) {
  num_annots <- sapply(predicted_impact, function(p) {
    length(unlist(strsplit(p, "[|]")))
  })
  num_genes <- sapply(gene, function(g) {
    length(unlist(strsplit(g, "[|]")))  })
  if (sum(unname(num_annots) != unname(num_genes)) != 0) {
    stop("Number of annotations do not match the number of genes at at least 1 position")
  }
}
