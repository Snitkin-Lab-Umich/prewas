#' Make tree edges positive
#'
#' If tree edges are zero or negative, makes them a very small positive number
#' (1/1000th of the smallest edge length). This is to prevent
#' \code{\link[ape]{ace}} from breaking. A warning is thrown. Returns a tree
#' with all positive edge lengths.
#'
#' @param tree (\code{ape phylo}).
#'
#' @return tree: (\code{ape phylo}). Tree now has only positive edge lengths.
#' @noRd
make_all_tree_edges_positive <- function(tree){
  check_is_tree(tree)
  if (sum(tree$edge.length <= 0) > 0) {
    warning(paste0("All non-positive branch lengths changed to small positive ",
                   "number to be able to perform ancestral reconstruction."))
    # Change any edge lengths that are zero to a very small number (so ancestral
    # reconstruction doesn't break)
    tree$edge.length[tree$edge.length <= 0] <-
      min(tree$edge.length[tree$edge.length > 0]) / 1000
  }
  return(tree)
}

#' Get major alleles
#'
#' Finds the major allele (most common allele) for each variant position in an
#' allele matrix.
#'
#' @param allele_mat Matrix. Allele matrix. Rows are variants. Columns are
#'   samples.
#'
#' @return major_allele: Character vector. Gives the major allele for each
#'   position. Names are the genomic loci. Values are the nucleotides.
#' @noRd
get_major_alleles <- function(allele_mat){
  major_allele <- apply(allele_mat, 1, function(x) {
    names(which.max(table(x)))
  })
  names(major_allele) <- rownames(allele_mat)
  return(major_allele)
}

#' Get ancestral state of alleles
#'
#' Finds the most likely ancestral allele for each variant position in an allele
#' matrix using \code{\link[ape]{ace}}, given a rooted tree. Returns the most
#' likely ancestral allele, its probability, and the tree used to perform
#' ancestral reconstruction.
#'
#' @param tree (\code{ape phylo}). Rooted tree.
#' @param mat Matrix. Allele matrix. Rows are variants, columns are samples.
#'
#' @return list of two elements:
#'   \describe{
#'     \item{ar_results}{Data.frame. Data.frame of most likely ancestral allele
#'     for each row in allele matrix and probability that that is the ancestral
#'     state. Rows are variants. First column is the allele. Second column is
#'     the probability.}
#'     \item{tree}{phylo. Tree used for ancestral state reconstruction.}
#'   }
#' @noRd
get_ancestral_alleles <- function(tree, mat){
  check_is_tree(tree)
  check_is_this_class(mat, "matrix")
  check_setequal_tree_mat(tree$tip.label, colnames(mat))
  check_tree_is_rooted(tree)

  future::plan(future::multiprocess)

  # ORDER MATRIX TO MATCH TREE TIP LABELS
  mat <- mat[, tree$tip.label, drop = FALSE]

  tree <- make_all_tree_edges_positive(tree)

  # Get ancestral state of root
  ar_all <- t(future.apply::future_apply(mat, 1, function(tip_states) {
    tip_state <- unique(tip_states)
    if (length(tip_state) > 1) {
      ar <- ape::ace(x = tip_states, phy = tree, type = "discrete")
      states <- ar$lik.anc[1, ]
      tip_state <- names(states)[which.max(states)]
      prob <- states[which.max(states)]
      c(tip_state, prob)
    } else {
      c(tip_states, 1)
    }
  }))
  ar_all <- data.frame(ar_all)
  colnames(ar_all) <- c("ancestral_allele", "probability")
 ar_all$ancestral_allele <- as.factor(ar_all$ancestral_allele)

  return(list(ar_results = ar_all, tree = tree))
}

#' Remove unknown ancestral states
#'
#' Removes rows from variant matrix where the reference allele (ancestral allele
#' or major allele) is unknown (- or N).
#'
#' @param allele_mat Matrix. Allele matrix. Rows are variants, columns are
#'   samples.
#' @param alleles Factor. Vector of reference alleles (ancestral alleles or
#'   major alleles)
#' @param ar_results Data.frame. Results from ancestral reconstruction
#'
#' @return A list of three objects:
#'   \describe{
#'     \item{allele_mat}{Matrix. Rows are variants. Columns are samples.}
#'     \item{ar_results}{Data.frame. Variants with unknown ancestral states
#'     removed. Rows are variants. If ancestral reconstruction was performed:
#'     first column in ancestral allele & second column is probability. If no
#'     ancestral reconstruction performed: the only column is the major allele.}
#'     \item{removed}{Character. Vector with names of removed samples.}
#'   }
#' @noRd
remove_unknown_alleles <- function(allele_mat, alleles, ar_results, o_ref, o_alt, snpeff){
  check_is_this_class(allele_mat, "matrix")
  check_is_this_class(alleles, "factor")
  check_is_this_class(ar_results, "data.frame")

  unknown <- alleles %in% c("-", "N")
  removed <- rownames(allele_mat)[unknown]
  if (length(removed) > 0) {
    warning(paste(length(removed),
                  "positions removed because ancestral allele is unknown"))
  }
  return(list(allele_mat = allele_mat[!unknown, ],
              ar_results = ar_results[!unknown, , drop = FALSE],
              snpeff = snpeff[!unknown],
              o_ref = o_ref[!unknown],
              o_alt = o_alt[!unknown],
              removed = removed))
}

#' Make binary matrix from allele matrix
#'
#' Returns a binary matrix of variant presence/absence created from an allele
#' matrix and a reference allele (ancestral allele or major allele) vector that
#' can be used in bGWAS. The reference allele is 0 and the non-reference allele
#' is 1.
#'
#' @param allele_mat Matrix. Allele matrix (split by multi-allelic site). Rows
#'   are variants. Columns are samples.
#' @param reference_allele Factor. Vector of alleles that are 0 in binary matrix
#'   (ancestral allele or major allele). Named vector. Names are genetic loci.
#'
#' @return bin_mat. Matrix. Binary matrix of variant presence/absence.
#' @noRd
make_binary_matrix <- function(allele_mat, o_ref, n_ref, o_alt){
  #check_is_this_class(reference_allele, "factor")
  check_is_this_class(allele_mat, "matrix")

  # make matrix of reference allele that's the same size as the allele matrix
  ref_allele_mat <- replicate(ncol(allele_mat), n_ref)
  # initialize binary matrix
  bin_mat <- allele_mat
  # if allele is the reference allele, code it as 0
  bin_mat[bin_mat == ref_allele_mat] <- 0

  # initialize n_alt (new alternative allele)
  n_alt <-  rep(NA, length(o_alt))

  # assign new alternative (after new reference)
  for (i in 1:length(o_alt)) {
    if (n_ref[i] == o_alt[i]) {
      n_alt[i] <- o_ref[i]
    }else{
      n_alt[i] <- o_alt[i]
    }
  }

  # iterate over each row in bin_mat
  for (j in 1:nrow(bin_mat)) {
    bin_mat[j,bin_mat[j,] != n_alt[j]] <- 0
    bin_mat[j,bin_mat[j,] == n_alt[j]] <- 1
  }

  # update rownames and colnames of  binary matrix
  rownames(bin_mat) <- rownames(allele_mat)
  colnames(bin_mat) <- colnames(allele_mat)

  # make binary matrix numeric
  class(bin_mat) <- "numeric"
  return(list(bin_mat, n_alt))
}

#' Title
#'
#' @param alt_allele Character vector of alternative allele that matches the
#' snpeff annotation (original allele, not re-referenced alternative allele)
#' @param snpeff_split Character vector of snpeff annotations extracted from vcf
#'
#' @return list of predicted functional impact and gene names
#' (locus tag or symbols in the snpeff annotation)
#' @noRd
parse_snpeff <- function(alt_allele, snpeff_split){
  pred_impact = rep(NA, length(alt_allele))
  gene = rep(NA, length(alt_allele))

  for (i in 1:length(alt_allele)) {
    annot = snpeff_split[[i]][grep(paste0('^',alt_allele[i]), snpeff_split[[i]])]

    if (length(annot) == 1) {
      pred_impact[i] <- unlist(strsplit(annot, '[|]'))[3]
      gene[i] <- unlist(strsplit(annot, '[|]'))[4]
    }else if (length(annot) == 0){
      pred_impact[i] <- ""
      gene[i] <- ""
    }else{
      multi_annots <- unlist(strsplit(annot, ',')) # would occur if overlapping genes
      pred_impact[i] <- paste(sapply(multi_annots, function(s) {
        unlist(strsplit(s, '[|]'))[3]
      }), collapse = '|')
      gene[i] <- paste(sapply(multi_annots, function(s) {
        unlist(strsplit(s, '[|]'))[4]
      }), collapse = '|')
    }
  }

  return(list(pred_impact = pred_impact,
              gene = gene))

}
