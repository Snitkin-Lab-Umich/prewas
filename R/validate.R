#' Check that function input is a numeric.
#'
#' Doesn't return anything. Gives an error message if input is not a number.
#'
#' @param num Number
#'
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
#' @return is_this_class Logical
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
                 "DNAbin")
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
#' @return
#' @export
#'
#' @examples
check_is_this_class <- function(obj, current_class){
  if (class(current_class) != "character") {
    stop("Current_class is expected to be a string describing a class")
  }
  r_classes <- c("character",
                 "numeric",
                 "integer",
                 "logical",
                 "complex",
                 "phylo",
                 "DNAbin")
  if (!(current_class %in% r_classes)) {
    stop("current_class is expected to be a R class")
  }

  if (class(obj) != current_class) {
    stop(paste('Input must be a', current_class))
  }
}

#' Check that the input tree is actually a 'phylo' object.
#'
#' Doesn't return anything. Gives an error message if the object is not a
#' 'phylo' object.
#' @param tree Phylogenetic tree.
check_is_tree <- function(tree){
  if (!is_this_class(tree, 'phylo')) {
    stop('Input requires either a path to a tree file or an ape phylo object')
  }
}

check_tree_is_rooted <- function(tree){
  check_is_tree(tree)
  if (!ape::is.rooted(tree)) {
    stop('Tree must be rooted.')
  }
}

format_inputs <- function(dna, tree, outgroup = NULL, gff = NULL){
  # DNA
  if (is.null(dna)) {
    stop("Must include a VCF file")
  } else {
    if (is_file(dna)) {
      # If DNA stored in VCF file, read in and convert to allele matrix
      dna <- load_vcf_file(dna)
    } else {
      stop("DNA input should be a path to a VCF 4.1 file")
    }
  }

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
    if (is_file(gff)) {
      gff <- read_gff(gff)

      #subset gff on CDS
      gff <- subset_gff(gff)

      #get gene name from column 9
      gff <- clean_up_cds_name_from_gff(gff)
    } else {
      gff <- NULL
      print("Prewas only accepts a gff file path as the GFF input.")
    }

  }


  results <- list("dna" = dna,
                  "tree" = tree,
                  "outgroup" = outgroup,
                  "gff" = gff)
  return(results)
}

read_gff <- function(gff_path){
  # TODO add ability to load in a gff object, not just take in a path
  gff <- utils::read.table(gff_path,
                           sep = '\t',
                           header = FALSE,
                           stringsAsFactors = FALSE,
                           quote = "",
                           comment.char = '#'
  )

  if (ncol(gff) != 9) {
    stop("GFF file must have 9 columns")
  }
  if (sum(c(".", "+", "-") %in% unique(gff[, 7, drop = TRUE])) == 0) {
    stop("GFF file must only have strand information in column 7")
  }

  return(gff)
}

subset_gff <- function(gff){
  check_is_this_class(gff, "data.frame")
  # TODO add a check_dimenions() call and write the check_dimensions function

  gff <- gff[gff[, 3] == 'CDS', ] # subset gff on CDS
  if (nrow(gff) == 0) {
    stop("GFF has no CDS regions. Remove the gff from prewas inputs and rerun or use a different GFF with CDS regions.")
  }
  return(gff)
}

clean_up_cds_name_from_gff <- function(gff){
  check_is_this_class(gff, "data.frame")
  # TODO add a check_dimenions() call and write the check_dimensions function

  cds_name <- apply(gff, 1, function(row){
    gsub('^ID=', '', row[9]) %>% gsub(';.*$', '', .)
  })
  gff[, 9] <- cds_name
  return(gff)
}

read_dna <- function(fasta_path){
  # TODO -- do we still need this function?
  dna <- ape::read.dna(file = fasta_path, as.character = TRUE, format = "fasta")
  colnames(dna) <- 1:ncol(dna)
  dna <- t(dna)
  return(dna)
}

convert_dnabin_to_matrix <- function(dnabin){
  # TODO -- do we still need this function?
  dna <- ape::as.alignment(dnabin)
  dna_mat <- matrix(NA, nrow = nrow(dnabin), ncol = nchar(dna$seq[1]))
  for (i in 1:nrow(dna_mat)) {
    dna_mat[i, ] <- strsplit(dna$seq[[i]], "")[[1]]
  }
  row.names(dna_mat) <- row.names(dnabin)
  colnames(dna_mat) <- 1:ncol(dna_mat)
  dna_mat <- t(dna_mat)
  return(dna_mat)
}

load_vcf_file <- function(vcf_path) {
  vcf <- vcfR::read.vcfR(file = vcf_path)
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
