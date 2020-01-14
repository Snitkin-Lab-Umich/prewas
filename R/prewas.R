#' Preprocess SNPs before GWAS
#'
#' @description Prepocess multiallelic sites, variants in overlapping genes, and
#'   reference all variants to either the ancestral or major allele.
#'
#' @param dna Character. Required input. A path to to a VCF4.1 file.
#' @param tree Phylo or Character. Optional input. Either a phylogenetic tree or
#'   a path to a file containing a phylogenetic tree. Defaults to NULL.
#' @param outgroup Character. Optional input. Either a path to a file containing
#'   only the outgroup name or a string of the outgroup name. Must be a tip in
#'   the phylogenetic tree. Defaults to NULL.
#' @param gff Character. Optional input. Path to GFF3 file location. Defaults to
#'   NULL.
#' @param anc Logical. Optional input. When TRUE prewas performs ancestral
#'   reconstruction. When FALSE prewas calculates the major allele. Defaults to
#'   TRUE.
#'
#' @return An list with the following items:
#'   \describe{
#'     \item{allele_mat}{Allele matrix.}
#'     \item{bin_mat}{Binary matrix.}
#'     \item{ar_results}{TODO}
#'     \item{dup}{TODO}
#'     \item{gene_mat}{Gene matrix.}
#'   }
#' @export
#' @examples
#' \dontrun{
#' vcf = prewas::vcf
#' gff = prewas::gff
#' tree = prewas::tree
#' outgroup = prewas::outgroup
#' prewas(dna = vcf,
#'        tree = tree,
#'        outgroup = outgroup,
#'        gff = gff,
#'        anc = TRUE)
#' }
prewas <- function(dna,
                   tree = NULL,
                   outgroup = NULL,
                   gff = NULL,
                   anc = TRUE){
  # Check inputs ---------------------------------------------------------------
  inputs <- format_inputs(dna, tree, outgroup, gff, anc)
  # dna_mat: rows == variants & columns == isolates
  dna_mat <- inputs$dna
  tree <- inputs$tree
  outgroup_char <- inputs$outgroup
  gff_mat <- inputs$gff

  # preprocess tree and dna_mat ------------------------------------------------
  if(anc){
    tree <- root_tree(tree, outgroup_char)
    subsetted_data <- subset_tree_and_matrix(tree, dna_mat)
    tree <- subsetted_data$tree
    allele_mat_only_var <- keep_only_variant_sites(subsetted_data$mat)
  }else{
    check_is_this_class(dna_mat, 'matrix')
    allele_mat_only_var <- keep_only_variant_sites(dna_mat)
  }

  # ancestral reconstruction ---------------------------------------------------
  if (anc) {
    anc_alleles <- get_ancestral_alleles(tree, allele_mat_only_var)
    allele_results <- anc_alleles$ar_results
    tree <- anc_alleles$tree
    remove_unknown_anc_results <-
      remove_unknown_alleles(allele_mat_only_var,
                             allele_results$ancestral_allele,
                             allele_results)
    allele_mat_only_var <- remove_unknown_anc_results$allele_mat
    allele_results <- remove_unknown_anc_results$ar_results
  } else {
    alleles <- get_major_alleles(allele_mat_only_var)
    tree <- NULL
    allele_results <- data.frame(major_allele = alleles)
    rownames(allele_results) <- rownames(allele_mat_only_var)
    remove_unknown_anc_results <-
      remove_unknown_alleles(allele_mat_only_var,
                             allele_results$major_allele,
                             allele_results)
    allele_mat_only_var <- remove_unknown_anc_results$allele_mat
    allele_results <- remove_unknown_anc_results$ar_results
  }

  # split multiallelic snps ----------------------------------------------------
  split <- split_multi_to_biallelic_snps(mat = allele_mat_only_var,
                                        ar_results = allele_results)

  allele_mat_split <- split$mat_split
  allele_results_split <- split$ar_results_split
  split_rows_flag <- split$split_rows_flag

  if (anc) {
    alleles <- allele_results_split$ancestral_allele
    names(alleles) <- rownames(allele_results_split)
  } else {
    alleles <- allele_results_split$major_allele
    names(alleles) <- rownames(allele_results_split)
  }

  # reference to ancestral state -----------------------------------------------
  bin_mat <- make_binary_matrix(allele_mat_split, alleles)

  if (!is.null(gff_mat)) {
    # overlapping genes --------------------------------------------------------
    bin_mat <- dup_snps_in_overlapping_genes(bin_mat, gff_mat)

    # collapse snps by gene ----------------------------------------------------
    gene_names <- get_gene_names(bin_mat)
    gene_mat <- collapse_snps_into_genes(bin_mat, gene_names)
  } else {
    gene_mat = NULL
  }

  return(list(allele_mat = allele_mat_split,
              bin_mat = bin_mat,
              ar_results = allele_results_split,
              dup = split_rows_flag,
              gene_mat = gene_mat,
              tree = tree))
}
