#' Format user-provided inputs to prewas().
#'
#' @param dna Character or vcfR. Path to VCF file or vcfR object.
#' @param tree NULL, character, or phylo. If character it should be a path to
#'   the tree file.
#' @param outgroup NULL or character. If character it should be either a string
#'   naming the outgroup in the dataset or a path to a file containing the
#'   outgroup name.
#' @param gff NULL, character, or matrix. If character it should be a path to
#'   the GFF file. If NULL ignored. If a matrix or data.frame it should be the
#'   GFF information stored in 9 columns with the genes as rows.
#' @param anc logical. If TRUE ancestral reconstruction is performed and the
#'   ancestral allele is used for referencing. If FALSE no ancestral
#'   reconstruction is performed and the major allele is used for referencing.
#'
#' @return A list with the following four elements:
#'   \describe{
#'     \item{dna}{vcfR}
#'     \item{tree}{phylo}
#'     \item{outgroup}{character}
#'     \item{gff}{data.frame}
#'   }
format_inputs <- function(dna, tree = NULL, outgroup = NULL, gff = NULL, anc = TRUE){
  # DNA
  if (is.null(dna)) {
    stop("User must provide a vcfR object or path to a VCF 4.1 file")
  }

  # Convert vcfR object or VCF file to allele matrix
  dna <- load_vcf_file(dna)

  # Tree
  if (is.null(tree)) {
    if(anc) tree <- build_tree(dna)
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
    # Read in file or load matrix and then format GFF matrix
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
