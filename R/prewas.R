prewas <- function(dna,
                   tree,
                   outgroup = NULL,
                   gff = NULL,
                   out_prefix = "prewas"){

  # Check inputs
  inputs <- check_inputs(dna, tree, outgroup, gff)

  # DNA: as a matrix with rows == variants & columns == isolates
  dna_mat <- inputs$dna

  tree <- inputs$tree
  outgroup_char <- inputs$outgroup
  gff_mat <- inputs$gff

  # preprocess tree and fasta

  # generate allele matrix
  allele_mat_only_var <- generate_allele_matrix(dna_mat)

  # ancestral reconstruction

  # split multiallelic snps

  # refernce to ancestral state

  # overlapping genes

  # collapse snps by gene

}
