
check_setequal_tree_mat = function(tip_labels, colnames_mat){
  if(!setequal(tip_labels,colnames_mat)){
    stop('Tree and variant matrix sample names do not match.')
  }
}

#' Make tree edges positive
#'
#' If tree edges are zero or negative, makes them a very small positive number
#' (1/1000th of the smallest edge length). This is to prevent
#' \code{\link[ape]{ace}} from breaking. A warning is thrown. Returns a tree
#' with all positive edge lengths.
#'
#' @param tree (\code{ape phylo}) object
#'
#' @return tree (\code{ape phylo}) object with all positive edge lengths
#' @export
#'
#' @examples
make_all_tree_edges_positive = function(tree){
  if(sum(tree$edge.length <= 0) > 0){
    warning('All non-positive branch lengths changed to small positive number to be able to perform ancestral reconstruction.')
    # Change any edge lengths that are zero to a very small number (so ancestral reconstruction doesn't break)
    tree$edge.length[tree$edge.length <= 0] = min(tree$edge.length[tree$edge.length > 0])/1000
  }
  return(tree)
}

#' Get major alleles
#'
#' Finds the major allele (most common allele) for each variant position in an allele matrix.
#'
#' @param allele_mat allele matrix
#'
#' @return major allele for each position
#' @export
#'
#' @examples
get_major_alleles = function(allele_mat){
  major_allele = apply(allele_mat,1,function(x){
    names(which.max(table(x)))
  })
  names(major_allele) = rownames(allele_mat)
  return(major_allele)
}

#' Get ancestral state of alleles
#'
#' Finds the most likely ancestral allele for each variant position in an allele
#' matrix using \code{\link[ape]{ace}}, given a rooted tree. Returns the most
#' likely ancestral allele, its probability, and the tree used to perform
#' ancestral reconstruction.
#'
#' @param tree rooted tree (\code{ape phylo} object)
#' @param mat allele matrix (rows are variants, columns are samples)
#'
#' @return list of (1) ar_results: matrix of most likely ancestral allele for
#'   each row in allele matrix and probability that that is the ancestral state
#'   and (2) tree: tree used for ancestral state reconstruction
#' @export
#'
#' @examples
get_ancestral_alleles = function(tree,mat){
  future::plan(future::multiprocess)

  check_setequal_tree_mat(tree$tip.label, colnames(mat))
  check_tree_is_rooted(tree)

  # ORDER MATRIX TO MATCH TREE TIP LABELS
  mat = mat[,tree$tip.label]

  tree = make_all_tree_edges_positive(tree)

  # Get ancestral state of root
  ar_all = t(future.apply::future_apply(mat,1,function(tip_states){
    tip_state = unique(tip_states)
    if(length(tip_state) > 1){
      ar = ape::ace(x = tip_states,phy = tree,type = 'discrete')
      states = ar$lik.anc[1,]
      tip_state = names(states)[which.max(states)]
      prob = states[which.max(states)]
      c(tip_state,prob)
    }else{
      c(tip_states,1)
    }
  }))
  ar_all = data.frame(ar_all)
  colnames(ar_all) = c('ancestral_allele','probability')

  return(list(ar_results=ar_all,tree=tree))

}

#' Remove unknown ancestral states
#'
#' Removes rows from variant matrix where the reference allele (ancestral allele or major allele) is unknown (- or N)
#'
#' @param allele_mat  allele matrix (rows are variants, columns are samples)
#' @param alleles reference alleles (ancestral alleles or major alleles)
#' @param ar_results results from ancestral reconstruction
#'
#' @return
#' @export
#'
#' @examples
remove_unknown_alleles = function(allele_mat, alleles, ar_results){
  unknown = alleles %in% c('-','N')
  removed = rownames(allele_mat)[unknown]
  if (length(removed) > 0) {
    warning(paste(length(removed),'positions removed because ancestral allele is unknown'))
  }
  return(list(allele_mat = allele_mat[!unknown, ],
              ar_results = ar_results[!unknown, , drop = FALSE],
              removed = removed))
}

#' Make binary matrix from allele matrix
#'
#' Returns a binary matrix of variant presence/absence created from an allele
#' matrix and a reference allele (ancestral allele or major allele) vector that
#' can be used in bGWAS. The reference allele is 0 and the non-reference allele
#' is 1.
#'
#' @param allele_mat allele matrix (split by multi-allelic site)
#' @param reference_allele vector of alleles that are 0 in binary matrix
#'   (ancestral allele or major allele)
#'
#' @return binary matrix of variant presence/absence
#' @export
#'
#' @examples
make_binary_matrix = function(allele_mat,reference_allele){
  # make matrix of reference allele that's the same size as the allele matrix
  ref_allele_mat = replicate(ncol(allele_mat),reference_allele)
  # initialize binary matrix
  bin_mat = allele_mat
  # if allele is the reference allele, code it as 0
  bin_mat[bin_mat == ref_allele_mat] = 0
  # get variant positions
  sites = unique(gsub('\\..*','',rownames(bin_mat)))
  # iterate over each site (to handle multiallelic sites)
  bin_mat = sapply(sites, function(x){
    site = bin_mat[rownames(bin_mat) == x,] # get 1st site
    # get all bases at that position
    bases = unique(site)
    # remove reference (0) and unknown (N) bases
    bases = bases[bases != "0" & bases != "N"]
    # create mini-matrix for this variant position where rows are bases and columns are samples
    binsplit = matrix(NA,nrow=length(bases),ncol=ncol(bin_mat))
    rownames(binsplit) = bases
    # for each base, code that base as 1 and all others as 0
    for(b in bases){
      binsite = site
      binsite[binsite != b] = 0
      binsite[binsite == b] = 1
      binsite = as.numeric(binsite)
      binsplit[b,] = binsite
    }
    return(binsplit)
  })
  # change from list to matrix
  bin_mat = do.call(rbind,bin_mat)
  # update rownames and colnames of  binary matrix
  rownames(bin_mat) = rownames(allele_mat)
  colnames(bin_mat) = colnames(allele_mat)
  # make binary matrix numeric
  class(bin_mat) = 'numeric'
  return(bin_mat)
}

