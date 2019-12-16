#' Title
#'
#' @param dna
#' @param tree
#' @param outgroup
#' @param gff
#' @param out_prefix
#' @param anc Logical. Default to TRUE. When TRUE prewas performs ancestral reconstruction. When FALSE prewas calculates the major allele.
#'
#' @return
#' @export
#'
prewas <- function(dna,
                   tree = NULL,
                   outgroup = NULL,
                   gff = NULL,
                   anc = TRUE){

  # Check inputs ---------------------------------------------------------------
  inputs <- format_inputs(dna, tree, outgroup, gff)

  # dna_mat: rows == variants & columns == isolates
  dna_mat <- inputs$dna

  tree <- inputs$tree
  outgroup_char <- inputs$outgroup
  gff_mat <- inputs$gff

  # preprocess tree and dna_mat ------------------------------------------------
  tree = root_tree(tree, outgroup_char)
  subsetted_data = subset_tree_and_matrix(tree, dna_mat)
  tree = subsetted_data$tree
  allele_mat_only_var <- keep_only_variant_sites(subsetted_data$mat)

  # ancestral reconstruction ---------------------------------------------------
  if (anc) {
    anc_alleles = get_ancestral_alleles(tree, allele_mat_only_var)
    allele_results = anc_alleles$ar_results
    tree = anc_alleles$tree
    remove_unknown_anc_results = remove_unknown_alleles(allele_mat_only_var,allele_results$ancestral_allele,allele_results)
    allele_mat_only_var = remove_unknown_anc_results$allele_mat
    allele_results = remove_unknown_anc_results$ar_results
  }else{
    alleles = get_major_alleles(allele_mat_only_var)
    tree = NULL
    allele_results = data.frame(major_allele = alleles)
    rownames(allele_results) = rownames(allele_mat_only_var)
    remove_unknown_anc_results = remove_unknown_alleles(allele_mat_only_var,allele_results$major_allele,allele_results)
    allele_mat_only_var = remove_unknown_anc_results$allele_mat
    allele_results = remove_unknown_anc_results$ar_results
  }

  # split multiallelic snps ----------------------------------------------------
  split = split_multi_to_biallelic_snps(mat = allele_mat_only_var, ar_results = allele_results)

  allele_mat_split = split$mat_split
  allele_results_split = split$ar_results_split
  split_rows_flag = split$split_rows_flag

  if (anc) {
    alleles = allele_results_split$ancestral_allele
    names(alleles) = rownames(allele_results_split)
  } else {
    alleles = allele_results_split$major_allele
    names(alleles) = rownames(allele_results_split)
  }

  # reference to ancestral state ------------------------------------------------
  bin_mat = make_binary_matrix(allele_mat_split, alleles)

  if(!is.null(gff_mat)){
    # overlapping genes ----------------------------------------------------------
    bin_mat = dup_snps_in_overlapping_genes(bin_mat, gff_mat)

    # collapse snps by gene ------------------------------------------------------
    gene_names <- get_gene_names(bin_mat)
    gene_mat <- collapse_snps_into_genes(bin_mat, gene_names)
  }else{
    gene_mat = NULL
  }

  return(list(allele_mat = allele_mat_split,
              bin_mat = bin_mat,
              ar_results = allele_results_split,
              dup = split_rows_flag,
              gene_mat = gene_mat))
}
