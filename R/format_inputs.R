#' Format and check validity of user-provided prewas() inputs.
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
#' @param snpeff_grouping `NULL`, `character`. Optional input. Only used when a
#'   snpeff annotated multivcf is inputted. Use when you want to group SNPs by
#'   gene and snpeff impact. If `NULL` no custom-grouped gene matrix will be
#'   generated. Options for input are a vector combination of 'HIGH',
#'   'MODERATE', LOW', 'MODIFER'. Must write the impact combinations in all caps
#'   (e.g. c('HIGH', 'MODERATE')). Defaults to `NULL`.
#' @param grp_nonref `Logical`. Optional input. When `TRUE` prewas collapses all
#'   non-reference alleles for multi-allelic sites. When `FALSE` prewas keeps
#'   multi-allelic sites separate. Defaults to `FALSE`.
#'
#' @noRd
#' @return A list with the following four elements:
#'   \describe{
#'     \item{dna}{vcfR}
#'     \item{tree}{phylo}
#'     \item{outgroup}{character}
#'     \item{gff}{data.frame}
#'     \item{o_ref}{Original reference alleles. Character vector.}
#'     \item{o_alt}{Original alternative alleles from VCF. Character vector.}
#'     \item{snpeff}{VCF provided snpeff annotations.}
#'   }
format_inputs <- function(dna,
                          tree = NULL,
                          outgroup = NULL,
                          gff = NULL,
                          anc = TRUE,
                          snpeff_grouping = NULL,
                          grp_nonref = FALSE){
  # DNA
  if (is.null(dna)) {
    stop("User must provide a vcfR object or path to a VCF 4.1 file")
  }

  # Convert vcfR object or VCF file to allele matrix
  vcf_info <- load_vcf_file(dna)
  dna <- vcf_info$vcf_geno_mat
  o_ref <- vcf_info$vcf_ref_allele
  o_alt <- vcf_info$vcf_alt_allele
  snpeff <- vcf_info$snpeff_pred


  # Tree
  if (is.null(tree)) {
    if (anc) tree <- build_tree(dna)
  } else {
    if (is_file(tree)) {
      # If tree stored as .tree file, read in
      tree <- ape::read.tree(tree)
    }
    if (!is_this_class(tree, "phylo")) {
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

  # Check that group non-reference is logical
  if (!is.logical(grp_nonref)) {
    stop("User must provide a logical grp_nonref value")
  }

  check_snpeff_user_input(snpeff_grouping)

  results <- list("dna" = dna,
                  "tree" = tree,
                  "outgroup" = outgroup,
                  "gff" = gff,
                  "o_ref" = o_ref,
                  "o_alt" = o_alt,
                  "snpeff" = snpeff)
  return(results)
}
