#' Preprocess SNPs before bGWAS
#'
#' @description prewas is a tool to standardize the pre-processing of your
#'   genomic data before performing a bacterial genome-wide association study
#'   (bGWAS). prewas creates a variant matrix (where each row is a variant, each
#'   column is a sample, and the entries are presence - 1 - or absence - 0 - of
#'   the variant) that can be used as input for bGWAS tools. When creating the
#'   binary variant matrix, prewas can perform 3 pre-processing steps including:
#'   dealing with  multiallelic SNPs, (optional) dealing with SNPs in
#'   overlapping genes, and choosing a reference allele. prewas can output
#'   matrices for use with both SNP-based bGWAS and gene-based bGWAS.
#'
#' @param dna `Character` or `vcfR`. Required input. Path to VCF4.1 file or
#'   `vcfR` object.
#' @param tree `NULL`, `character`, or `phylo`. Optional input. Ignored if
#'   `NULL`. If `character` it should be a path to a .tree file. Defaults to
#'   `NULL`.
#' @param outgroup `NULL` or `character`. Optional input. If `character` it
#'   should be either a string naming the outgroup in the tree or a path to a
#'   file containing only the outgroup name. Ignored if `NULL`. Defaults to
#'   `NULL`.
#' @param gff `NULL`, `character`, `matrix`, or `data.frame`. Optional input. If
#'  `NULL` it is ignored. If `character` it should be a path to a GFF3 file. If
#'   a `matrix` or `data.frame` it should be the GFF information stored in 9
#'   columns with the genes as rows. Defaults to `NULL`.
#' @param anc `Logical`. Optional input. When `TRUE` prewas performs ancestral
#'   reconstruction. When `FALSE` prewas calculates the major allele. Defaults
#'   to `TRUE`.
#' @param grp_nonref `Logical`. Optional input. When `TRUE` prewas collapses all
#'   non-reference alleles for multi-allelic sites. When `FALSE` prewas keeps
#'   multi-allelic sites separate. Defaults to `FALSE`.
#'
#' @return A list with the following items:
#'   \describe{
#'     \item{allele_mat}{`matrix`. An allele matrix, created from the vcf where
#'     each multiallelic site will be on its own line. The rowname will be the
#'     position of the variant in the vcf file. If the position is triallelic,
#'     there will be two rows containing the same information. The rows will be
#'     labeled "pos" and "pos.1". If the position is quadallelic, there will be
#'     three rows containing the same information. The rows will be labeled
#'     "pos", "pos.1", and "pos.2"}
#'     \item{bin_mat}{`matrix`. A binary matrix, where 0 is the reference allele
#'     and 1 indicates a variant. The dimensions may not match the `allele_mat`
#'     if the gff file is provided, because SNPs in overlapping genes are
#'     represented on multiple lines in `bin_mat`; in that case both position
#'     and locus tag name are provided in the rowname.}
#'     \item{ar_results}{`data.frame`. This data.frame records the alleles used
#'     as the reference alleles. Rows correspond to variant loci. If `anc =
#'     TRUE` the data.frame has two columns which contain the ancestrally
#'     reconstructed allele and the probability of the reconstruction. If `anc =
#'     FALSE` there is only one column which contains the major allele.}
#'     \item{dup}{`integer`. A vector of integers. It's an index that identifies
#'     duplicated rows. If the index is unique (appears once), that means it is
#'     not a multiallelic site. If the index appears more than once, that means
#'     the row was replicated `x` times, where `x` is the number of alternative
#'     alleles. Note: the multiple indices indicates multiallelic site splits,
#'     not overlapping genes splits.}
#'     \item{gene_mat}{`NULL` or `matrix`. `NULL` if no gene information
#'     provided (`gff = NULL`). If gene information is provided, a gene matrix
#'     is generated where each row is a gene and each column is a sample.}
#'     \item{tree}{`NULL` or `phylo`. If `anc = FALSE` no tree is use or
#'     generated and the function returns `NULL`. If `anc = TRUE` and the user
#'     provides a tree but no outgroup: the function returns the tree after
#'     midpoint rooting. If `anc = TRUE` and the user provides both a tree and
#'     an outgroup: the function returns a tree rooted on the outgroup and the
#'     outgroup is removed from the tree. If the user does not provide a tree
#'     and `anc = TRUE` the function returns the midpoint rooted tree
#'     generated.}
#'   }
#' @export
#' @examples
#' vcf = prewas::vcf
#' gff = prewas::gff
#' tree = prewas::tree
#' outgroup = prewas::outgroup
#' output <- prewas(dna = vcf,
#'                  tree = tree,
#'                  outgroup = outgroup,
#'                  gff = gff,
#'                  anc = FALSE ,
#'                  grp_nonref = FALSE)
prewas <- function(dna,
                   tree = NULL,
                   outgroup = NULL,
                   gff = NULL,
                   anc = TRUE,
                   grp_nonref = FALSE){
  # Check inputs ---------------------------------------------------------------
  inputs <- format_inputs(dna, tree, outgroup, gff, anc)
  dna_mat <- inputs$dna
  tree <- inputs$tree
  outgroup_char <- inputs$outgroup
  gff_mat <- inputs$gff
  o_ref <- inputs$o_ref
  o_alt <- inputs$o_alt
  snpeff <- inputs$snpeff

  # preprocess tree and dna_mat ------------------------------------------------
  if (anc) {
    tree <- root_tree(tree, outgroup_char)
    subsetted_data <- subset_tree_and_matrix(tree, dna_mat)
    tree <- subsetted_data$tree
    allele_mat_snpeff_only_var <- keep_only_variant_sites(subsetted_data$mat, o_ref, o_alt, snpeff)
    allele_mat_only_var <- allele_mat_snpeff_only_var$variant_only_dna_mat
    o_ref_only_var <- allele_mat_snpeff_only_var$o_ref_var_pos
    o_alt_only_var <-  allele_mat_snpeff_only_var$o_alt_var_pos
    snpeff_only_var <- allele_mat_snpeff_only_var$snpeff_var_pos
  } else {
    check_is_this_class(dna_mat, "matrix")
    allele_mat_snpeff_only_var <- keep_only_variant_sites(dna_mat, o_ref, o_alt, snpeff)
    allele_mat_only_var <- allele_mat_snpeff_only_var$variant_only_dna_mat
    o_ref_only_var <- allele_mat_snpeff_only_var$o_ref_var_pos
    o_alt_only_var <-  allele_mat_snpeff_only_var$o_alt_var_pos
    snpeff_only_var <- allele_mat_snpeff_only_var$snpeff_var_pos
  }

  # ancestral reconstruction ---------------------------------------------------
  if (anc) {
    anc_alleles <- get_ancestral_alleles(tree, allele_mat_only_var)
    allele_results <- anc_alleles$ar_results
    tree <- anc_alleles$tree
    remove_unknown_anc_results <-
      remove_unknown_alleles(allele_mat_only_var,
                             allele_results$ancestral_allele,
                             allele_results,
                             o_ref_only_var,
                             o_alt_only_var,
                             snpeff)
    allele_mat_only_var <- remove_unknown_anc_results$allele_mat
    allele_results <- remove_unknown_anc_results$ar_results
    o_ref <- remove_unknown_anc_results$o_ref
    o_alt <- remove_unknown_anc_results$o_alt
    snpeff <- remove_unknown_anc_results$snpeff
  } else {
    alleles <- get_major_alleles(allele_mat_only_var)
    tree <- NULL
    allele_results <- data.frame(major_allele = alleles)
    rownames(allele_results) <- rownames(allele_mat_only_var)
    allele_results$major_allele <- as.factor(allele_results$major_allele)
    remove_unknown_anc_results <-
      remove_unknown_alleles(allele_mat_only_var,
                             allele_results$major_allele,
                             allele_results,
                             o_ref_only_var,
                             o_alt_only_var,
                             snpeff)
    allele_mat_only_var <- remove_unknown_anc_results$allele_mat
    allele_results <- remove_unknown_anc_results$ar_results
    o_ref <- remove_unknown_anc_results$o_ref
    o_alt <- remove_unknown_anc_results$o_alt
    snpeff <- remove_unknown_anc_results$snpeff
  }

  print(o_alt)
  # split multiallelic snps ----------------------------------------------------
  split <- split_multi_to_biallelic_snps(mat = allele_mat_only_var,
                                         o_ref = o_ref,
                                         o_alt = o_alt,
                                         ar_results = allele_results,
                                         snpeff = snpeff)

  allele_mat_split <- split$mat_split
  allele_results_split <- split$ar_results_split
  snpeff_split <- split$snpeff_split
  split_rows_flag <- split$split_rows_flag
  o_ref_split <- split$o_ref_split
  o_alt_split <- split$o_alt_split

  if (anc) {
    n_ref <- allele_results_split$ancestral_allele
    names(n_ref) <- rownames(allele_results_split)
  } else {
    n_ref <- allele_results_split$major_allele
    names(n_ref) <- rownames(allele_results_split)
  }
  print(o_ref_split)
  print(n_ref)
  print(o_alt_split)
  # reference to reference allele ----------------------------------------------
  bin_mat_results <- make_binary_matrix(allele_mat_split, o_ref_split, n_ref, o_alt_split)
  bin_mat <- bin_mat_results[[1]]
  n_alt <- bin_mat_results[[2]] # new alternative allele, after splitting

  # snpeff ---------------------------------------------------------------------
  snpeff_parsed <- parse_snpeff(o_alt_split, snpeff_split)
  predicted_impact <- snpeff_parsed$pred_impact
  gene <- snpeff_parsed$gene
  print(predicted_impact)
  print(gene)

  if (!is.null(snpeff_split)) {
    if (!is.null(gff_mat)) {
      print('Warning: GFF is not being used; annotations are coming from vcf file')
    }

    output <- dup_snps_in_overlapping_genes_snpeff(bin_mat, predicted_impact, gene)
    bin_mat <- output$bin_mat_dup
    predicted_impact_split <- output$predicted_impact_split

    # collapse snps by gene and impact -----------------------------------------
    gene_names <- get_gene_names(bin_mat)
    gene_mat <- collapse_snps_into_genes_by_impact(bin_mat, gene_names, predicted_impact_split)
  }else if (!is.null(gff_mat)) {
    # overlapping genes --------------------------------------------------------
    bin_mat <- dup_snps_in_overlapping_genes(bin_mat, gff_mat)

    # collapse snps by gene ----------------------------------------------------
    gene_names <- get_gene_names(bin_mat)
    gene_mat <- collapse_snps_into_genes(bin_mat, gene_names)
  } else {
    gene_mat <- NULL
  }

  if (grp_nonref) {
    # collapse variants by position --------------------------------------------
    allele_names <- get_allele_names(bin_mat)
    bin_mat <- collapse_snps_into_genes(bin_mat, allele_names)
  }

  return(list(allele_mat = allele_mat_split,
              bin_mat = bin_mat,
              ar_results = allele_results_split,
              dup = split_rows_flag,
              gene_mat = gene_mat,
              tree = tree))
}
