#' Generate multiple phylogenetic trees to use as test data for package.
#'
#' @param seed Integer to act as seed for random tree generation process.
#'
#' @return
#' @export
#'
#' @examples
generate_test_trees <- function(seed){
  # TODO add input check

  # Make tree
  set.seed(seed)
  num_samples <- 11
  tree_with_og <- ape::rcoal(n = num_samples)

  # Unrooted with outgroup
  tree_with_og_unrooted <- ape::unroot(tree_with_og)

  # Root with outgroup
  tree_with_og_rooted <- ape::root(phy = tree_with_og_unrooted, outgroup = "t1")

  # Keep root but remove outgroup
  tree_no_og_rooted <- ape::drop.tip(tree_with_og_rooted, tip = "t1")

  # No root or outgroup
  tree_no_og_unrooted <- ape::unroot(tree_no_og_rooted)

  return(list("tree_with_og_unrooted" = tree_with_og_unrooted,
              "tree_with_og_rooted" = tree_with_og_rooted,
              "tree_no_og_rooted" = tree_no_og_rooted,
              "tree_no_og_unrooted" = tree_no_og_unrooted))
}

#' Given a phylogenetic tree generate plausible DNA sequences for each tip.
#'
#' @param tree
#'
#' @return
#' @export
#'
#' @examples
generate_test_dna <- function(tree){
  num_nucleotides <- 4
  return_anc_state <- TRUE
  dna_seq <- ape::rTraitDisc(tree,
                             model = "ARD",
                             k = generate_test_dna,
                             states = c("A", "T", "G", "C"),
                             p = num_nucleotides,
                             ancestor = return_anc_state)
}
