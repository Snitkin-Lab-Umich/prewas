#' Check that function input is a numeric.
#'
#' Doesn't return anything. Gives an error message if input is not a number.
#'
#' @param num Number
#' @examples
#' check_is_number(100)
check_is_number <- function(num){
  if (!class(num) %in% c("numeric")) {
    stop("Input must be a numeric")
  }
}

#' Tests if input object is a file.
#'
#' @param obj Any R object.
#'
#' @return is_file Logical
is_file <- function(obj){
  is_file <- FALSE
  if (class(obj) == "character") {
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
#' @examples
#' object <- "example"
#' is_this_class(object, "character")
is_this_class <- function(obj, current_class){
  if (class(current_class) != "character") {
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
                 "vcfR")
  if (!(current_class %in% r_classes)) {
    stop("current_class is expected to be a R class")
  }

  is_this_class <- FALSE
  if (class(obj) == current_class) {
    is_this_class <- TRUE
  }
  return(is_this_class)
}

#' Checks if an object is of the expected R class.
#'
#' Doesn't return anything. Gives error if the object is not of the expected R
#' class.
#'
#' @param obj Any R object.
#' @param current_class Character string. Name of R class
#'
#' @export
#'
#' @examples
#' object <- "example"
#' check_is_this_class(object, "character")
check_is_this_class <- function(obj, current_class){
  class_log <- is_this_class(obj, current_class)
  if (class_log != TRUE) {
    stop(paste('Input must be a', current_class))
  }
}

#' Check that the input tree is actually a 'phylo' object.
#'
#' Doesn't return anything. Gives an error message if the object is not a
#' 'phylo' object.
#'
#' @param tree Phylogenetic tree.
#'
#' @examples
#' tree <- ape::rtree(10)
#' check_is_tree(tree)
check_is_tree <- function(tree){
  if (!is_this_class(tree, 'phylo')) {
    stop('Input requires either a path to a tree file or an ape phylo object')
  }
}

#' Check that the tree has a root.
#'
#' Doesn't return anything. Gives an error message if the object is not a
#' rooted tree.
#'
#' @param tree Phylogenetic tree.
#'
#' @examples
#' tree <- ape::rtree(10)
#' check_tree_is_rooted(tree)
check_tree_is_rooted <- function(tree){
  check_is_tree(tree)
  if (!ape::is.rooted(tree)) {
    stop('Tree must be rooted.')
  }
}

#' Format user-provided inputs to prewas().
#'
#' @param dna Character or vcfR. Path to VCF file or vcfR object.
#' @param tree NULL, character, or phylo. If character it should be a path to
#'   the tree file.
#' @param outgroup NULL or character. If character it should be either a string
#'   naming the outgroup in the dataset or a path to a file containing the
#'   outgroup name.
#' @param gff NULL or character. If character it should be a path to the GFF
#'   file.
#'
#' @return A list with the following four elements:
#'   \describe{
#'     \item{dna}{vcfR}
#'     \item{tree}{phylo}
#'     \item{outgroup}{character}
#'     \item{gff}{data.frame}
#'   }
format_inputs <- function(dna, tree = NULL, outgroup = NULL, gff = NULL){
  # DNA
  if (is.null(dna)) {
    stop("User must provide a vcfR object or path to a VCF 4.1 file")
  }

  # Convert vcfR object or VCF file to allele matrix
  dna <- load_vcf_file(dna)

  # Tree
  if (is.null(tree)) {
    tree <- build_tree(dna)
  } else {
      if (is_file(tree)) {
        # If tree stored as .tree file, read in
        tree <- ape::read.tree(tree)
      }
      if (!is_this_class(tree, 'phylo')) {
        stop("Tree should either be a .tree or a phylo object")
      }
  }

 # Outgroup
  if (!is.null(outgroup)) {
    # If file, read it in
    if (is_file(outgroup)) {
      outgroup <- readLines(con = outgroup)
      outgroup <- as.character(outgroup)
    } else {
    # If object, check is a character string
      if (!is_this_class(outgroup, "character")) {
        stop("Outgroup must be a character string")
      }
    }

    # Check if outgroup is in the tree
    if (!outgroup %in% tree$tip.label) {
      warning("Outgroup not found in tree so ignoring it.")
      outgroup <- NULL
    }
  }

  # GFF
  if (!is.null(gff)) {
    # Read in file or data.frame and then format GFF data.frame
    gff <- read_gff(gff)

    # Subset GFF to only CDS
    gff <- subset_gff(gff)

    # Remove ID= prefix from gene name in column 9
    gff <- clean_up_cds_name_from_gff(gff)
  }

  results <- list("dna" = dna,
                  "tree" = tree,
                  "outgroup" = outgroup,
                  "gff" = gff)
  return(results)
}

#' Read in a GFF file given a file path
#'
#' @param gff Character or data.frame. If character: path to GFF file.
#'
#' @return gff: Data.frame.
#' @export
#'
read_gff <- function(gff){
  if (is_file(gff)) {
    gff <- utils::read.table(gff_path,
                             sep = '\t',
                             header = FALSE,
                             stringsAsFactors = FALSE,
                             quote = "",
                             comment.char = '#')
  }

  check_is_this_class(gff, "data.frame")

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
#' @param gff Data.frame. Rows = genomic features.
#'
#' @return gff: Data.frame.
#' @export
#'
subset_gff <- function(gff){
  check_is_this_class(gff, "data.frame")
  # TODO add a check_dimenions() call and write the check_dimensions function

  gff <- gff[gff[, 3] == 'CDS', ] # subset gff on CDS
  if (nrow(gff) == 0) {
    stop("GFF has no CDS regions. Remove the gff from prewas inputs and rerun or use a different GFF with CDS regions.")
  }
  return(gff)
}

#' Remove ID prefix from CDS feature descriptions (e.g. gene names).
#'
#' @param gff Data.frame. Rows = genomic features.
#'
#' @return gff: Data.frame.
#' @export
#'
clean_up_cds_name_from_gff <- function(gff){
  check_is_this_class(gff, "data.frame")
  # TODO add a check_dimenions() call and write the check_dimensions function

  cds_name <- apply(gff, 1, function(row){
    gsub('^ID=', '', row[9]) %>% gsub(';.*$', '', .)
  })
  gff[, 9] <- cds_name
  return(gff)
}

#' Read in variant data from a VCF file or vcfR object
#'
#' @param vcf Either character or vcfR. If character, it is a path to a VCF
#'   file.
#'
#' @return vcf_geno_mat: Matrix.
#' @export
#'
load_vcf_file <- function(vcf) {
  if (is_file(vcf)) {
    vcf <- vcfR::read.vcfR(file = vcf)
  }
  check_is_this_class(vcf, "vcfR")

  vcf_geno_mat <- vcf@gt[, 2:ncol(vcf@gt), drop = FALSE]
  row.names(vcf_geno_mat) <- vcf@fix[, colnames(vcf@fix) == "POS", drop = TRUE]
  vcf_ref_allele <- vcf@fix[, colnames(vcf@fix) == "REF", drop = TRUE]
  vcf_alt_allele <- vcf@fix[, colnames(vcf@fix) == "ALT", drop = TRUE]

  for (i in 1:nrow(vcf_geno_mat)) {
    alt_alleles <- strsplit(vcf_alt_allele[i], split = ",")
    vcf_alt_allele_1 <- alt_alleles[[1]][1]
    vcf_alt_allele_2 <- alt_alleles[[1]][2]
    vcf_alt_allele_3 <- alt_alleles[[1]][3]

    vcf_geno_mat[i, vcf_geno_mat[i, ] == "0"] <- vcf_ref_allele[i]
    vcf_geno_mat[i, vcf_geno_mat[i, ] == "1"] <- vcf_alt_allele_1
    vcf_geno_mat[i, vcf_geno_mat[i, ] == "2"] <- vcf_alt_allele_2
    vcf_geno_mat[i, vcf_geno_mat[i, ] == "3"] <- vcf_alt_allele_3
  }
  return(vcf_geno_mat)
}

#' Confirm that the tree and variant matrix contain exactly the same samples
#'
#' Doesn't return anything. Gives error if the two inputs do not match.
#'
#' @param tip_labels Character. Vector of tree$tip.labels.
#' @param colnames_mat Character. Vector of column names from variant matrix.
#'
#' @export
#'
check_setequal_tree_mat <- function(tip_labels, colnames_mat){
  if (!setequal(tip_labels, colnames_mat)) {
    stop('Tree and variant matrix sample names do not match.')
  }
}
