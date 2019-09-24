#' Generate multiple phylogenetic trees to use as test data for package.
#'
#' @param seed Integer to act as seed for random tree generation process.
#'
#' @return
#' @export
#'
#' @examples
generate_test_trees <- function(num_samples, seed){
  # Make tree
  set.seed(seed)
  tree_with_og <- ape::rcoal(n = num_samples)

  # Unrooted with outgroup
  tree_with_og_unrooted <- ape::unroot(tree_with_og)

  # Root with outgroup
  tree_with_og_rooted <- ape::root(phy = tree_with_og_unrooted, outgroup = "t1")

  # Keep root but remove outgroup
  tree_no_og_rooted <- ape::drop.tip(tree_with_og_rooted, tip = "t1")

  # No root or outgroup
  tree_no_og_unrooted <- ape::unroot(tree_no_og_rooted)

  # Tree tip labels are wrong
  tree_bad_labels <-  tree_no_og_rooted
  tree_bad_labels$tip.label <- letters[1:ape::Ntip(tree_bad_labels)]

  # Highly clonal
  clonal_tree <- tree_no_og_rooted
  clonal_tree$edge.length <- rnorm(ape::Nedge(clonal_tree), mean = 0.005, sd = 0.001)

  # Very diverse
  diverse_tree <- tree_no_og_rooted
  diverse_tree$edge.length <- rnorm(ape::Nedge(diverse_tree), mean = 1.25, sd = 0.25)


  return(list("tree_with_og_unrooted" = tree_with_og_unrooted,
              "tree_with_og_rooted" = tree_with_og_rooted,
              "tree_no_og_rooted" = tree_no_og_rooted,
              "tree_no_og_unrooted" = tree_no_og_unrooted,
              "tree_bad_labels" = tree_bad_labels,
              "clonal_tree" = clonal_tree,
              "diverse_tree" = diverse_tree))
}

#' Given a phylogenetic tree generate plausible DNA sequences for each tip.
#'
#' @param tree Phylogenetic tree.
#' @param seed Number. Seed for set.seed() function.
#' @param seq_length Number. Length of the generated DNA sequence.
#'
#' @return
#' @export
generate_test_dna <- function(tree, seq_length, seed){
  set.seed(seed)
  dna_seq <- phytools::genSeq(tree = tree, l = seq_length, format = "phyDat")
  return(dna_seq)
}

#' Save a fasta file given a phyDat sequence object.
#'
#' @param phydat_obj DNA sequence. phyDat object.
#' @param file_prefix String. File name (without .fna).
#'
#' @export
save_fasta <- function(phydat_obj, file_prefix){
  file_name <- paste0("data/", file_prefix, ".fna")
  phangorn::write.phyDat(x = phydat_obj, file = file_name, format = "fasta")
}

#' Save a tree to a .tree file.
#'
#' @param tree Phylogenetic tree.
#' @param file_prefix String. File name (without .fna).
#'
#' @export
save_tree <- function(tree, file_prefix){
  file_name <- paste0("data/", file_prefix, ".tree")
  ape::write.tree(phy = tree, file = file_name)
}

#' Given a dna sequence, generate a dummy GFF3 formatted file.
#'
#' @param phydat_obj phyDat object. DNA alignment.
#'
#' @return
#' @export
generate_test_gff <- function(phydat_obj, seq_length){
  num_rows <- 100
  avg_gene_length <- seq_length/num_rows

  gff <- matrix(".", ncol = 9, nrow = num_rows)
  # seqid
  gff[, 1] <- "chr1"

  # source
  gff[, 2] <- "prewas"

  # type
  gff[, 3] <- "CDS"

  # strand
  gff[seq(1, nrow(gff), 2), 7] <- "+"
  gff[seq(2, nrow(gff), 2), 7] <- "-"

  # start
  gff[, 4:5] <- seq(from = 1, to = seq_length, by = avg_gene_length)

  # end
  gff[, 5] <- as.numeric(gff[, 5]) + (avg_gene_length - 1)

  # ID
  gff[, 9] <- paste0("ID=gene_", 1:nrow(gff))

  # Now make some overlapping genes:
  gff[1:4, 4] <- 1
  gff[1:4, 5] <- 4 * avg_gene_length

  return(gff)
}

save_gff3 <- function(gff, file_prefix){
  first_row <- "##gff-version 3"
  file_name <- paste0("data/", file_prefix, ".gff")

  write(first_row, file = file_name)

  write.table(x = gff,
              file_name,
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE,
              append = TRUE)
}

append_fasta_to_gff <- function(gff_path) {
  fake_fasta <- ">ATGATGATGATGATGATG"
  write(fake_fasta, file = gff_path, append = TRUE)
}

generate_test_data <- function(num_samples, seq_length, seed){
  # Generate & save trees
  trees <- generate_test_trees(num_samples, seed)
  save_tree(trees$tree_with_og_unrooted, "og_unrooted")
  save_tree(trees$tree_with_og_rooted, "og_rooted")
  save_tree(trees$tree_no_og_unrooted, "no_og_unrooted")
  save_tree(trees$tree_no_og_rooted, "no_og_rooted")
  save_tree(trees$tree_bad_labels, "bad_labels")
  save_tree(trees$clonal_tree, "clonal")
  save_tree(trees$diverse_tree, "diverse")

  # Generate & save DNA alignments
  dna <- generate_test_dna(trees$tree_no_og_rooted, seq_length, seed)
  bad_labels_dna <- dna
  names(bad_labels_dna) <- letters[1:length(dna)]
  clonal_dna <- generate_test_dna(trees$clonal_tree, seq_length, seed)
  diverse_dna <- generate_test_dna(trees$diverse_tree, seq_length, seed)
  save_fasta(dna, "rooted_no_og")
  save_fasta(bad_labels_dna, "bad_labels")
  save_fasta(clonal_dna, "clonal")
  save_fasta(diverse_dna, "diverse")

  # Generate & save GFF files
  gff <- generate_test_gff(dna, seq_length)
  save_gff3(gff, "no_fasta_appended")
  save_gff3(gff, "fasta_appended")
  append_fasta_to_gff("data/fasta_appended.gff")
}

# generate_test_data(15, 1000, 1)

