
check_setequal_tree_mat = function(tip_labels, colnames_mat){
  if(!setequal(tip_labels,colnames_mat)){
    stop('Tree and variant matrix sample names do not match.')
  }
}

make_all_tree_edges_positive = function(tree){
  if(sum(tree$edge.length <= 0) > 0){
    warning('All non-positive branch lengths changed to small positive number to be able to perform ancestral reconstruction.')
    # Change any edge lengths that are zero to a very small number (so ancestral reconstruction doesn't break)
    tree$edge.length[tree$edge.length <= 0] = min(tree$edge.length[tree$edge.length > 0])/1000
  }
  return(tree)
}

#' Get_major_alleles
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
#' @description Rereference alleles based on rooted tree
#'
#' @param tree rooted tree
#' @param mat allele matrix (rows are variants, columns are samples)
#'
#' @return list of (1) ar_results: matrix of most likely ancestral allele for each row in allele matrix and probability that that is the ancestral state and (2) tree: tree used for ancestral state reconstruction
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
#' @description Remove rows from variant matrix where the ancestral state is unknown (- or N)
#'
#' @param allele_mat
#' @param alleles
#' @param ar_results
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
#' @param allele_mat allele matrix (split by multi-allelic site)
#' @param reference_allele vector of alleles that are 0 in binary matrix (ancestral allele or major allele)
#'
#' @return
#' @export
#'
#' @examples
make_binary_matrix = function(allele_mat,reference_allele){
  ref_allele_mat = replicate(ncol(allele_mat),reference_allele)
  bin_mat = allele_mat
  bin_mat[bin_mat == ref_allele_mat] = 0
  sites = unique(gsub('\\..*','',rownames(bin_mat)))
  #rownames(bin_mat) = gsub('\\..*','',rownames(bin_mat))
  bin_mat = sapply(sites, function(x){
    site = bin_mat[rownames(bin_mat) == x,] # get 1st site
    bases = unique(site)
    bases = bases[bases != "0" & bases != "N"]
    binsplit = matrix(NA,nrow=length(bases),ncol=ncol(bin_mat))
    rownames(binsplit) = bases
    for(b in bases){
      binsite = site
      binsite[binsite != b] = 0
      binsite[binsite == b] = 1
      binsite = as.numeric(binsite)
      binsplit[b,] = binsite
    }
    binsplit
  })
  bin_mat = do.call(rbind,bin_mat)
  rownames(bin_mat) = rownames(allele_mat)
  colnames(bin_mat) = colnames(allele_mat)


  class(bin_mat) = 'numeric'
  return(bin_mat)
}

#' Reference to ancestral state
#'
#' @param allele_mat allele matrix
#' @param anc_alleles ancestral alleles
#'
#' @return binary matrix referenced to ancestral state
#' @export
#'
#' @examples
# reference_alleles = function(allele_mat,alleles){
#     # reference to ancestral state
#     alleles = ar_results$ancestral_allele
#     names(alleles) = rownames(ar_results)
#     # check that columns in allele_mat match with rows in alleles
#     allele_mat = allele_mat[rownames(allele_mat),]
#     alleles = alleles[rownames(allele_mat)]
#     if(!all(rownames(allele_mat) == names(alleles))){
#       stop('Allele matrix columns and ancestral allele state names must match.')
#     }
#
#   # MAKE BINARY MATRIX
#   make_binary_matrix(allele_mat,alleles)
#
#   return(bin_mat)
# }



