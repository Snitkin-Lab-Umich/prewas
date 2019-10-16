
#' Read in tree
#'
#' @param tree tree path or ape phylo object
#'
#' @return ape phylo object
#' @export
#'
#' @examples
read_in_tree = function(tree){
  if(is_file(tree)){
    # LOAD IN TREE
    tree = read.tree(tree)
    return(tree)
  }else{
    check_is_tree(tree)
  }
  return(tree)
}

#' Root tree on outgroup
#' @description Root tree based on outgroup. If outgroup is provided, root based on outgroup. If outgroup is null and the tree isn't rooted, midpoint root tree. If tree is rooted and no outgroup is provided, return tree as is.
#'
#' @param tree phylogenetic tree or file path of tree
#' @param outgroup tip name of outgroup in phylogeny. If NULL, midpoint root if not rooted
#'
#' @return rooted tree without outgroup
#' @export
#'
#' @examples
#' tree = rcoal(100)
#' is.rooted(tree)
#' root_tree_og(tree)
#' is.rooted(tree)
root_tree = function(tree,outgroup=NULL){
  # READ IN TREE
  tree = read_in_tree(tree)
  # IF OUTGROUP SUPPLIED
  if(!is.null(outgroup)){
    # ROOT TREE ON OUTGROUP
    tree = ape::root(tree,outgroup)
    tree = ape::drop.tip(tree,outgroup)
    return(tree)
  # IF NO OUTGROUP AND TREE IS UNROOTED
  }else if(is.null(outgroup) & !ape::is.rooted(tree)){
    # MIDPOINT ROOT TREE
    tree = phytools::midpoint.root(tree)
    return(tree)
  }
  return(tree)
}

#' Subset tree and matrix
#'
#' @param tree path to tree file or ape phylo object
#' @param mat variant matrix (columns are samples, rows are variants)
#'
#' @return
#' @export
#'
#' @examples
subset_tree_and_matrix = function(tree,mat){
  tree = read_in_tree(tree)
  check_is_this_class(mat,'matrix')
  # return if tree tip labels and variant matrix column names match
  if(setequal(tree$tip.label,colnames(mat))){
    # order matrix same as tree tip labels
    mat = mat[,tree$tip.label]
    return(list(tree=tree,mat=mat))
  }
  in_both = intersect(tree$tip.label,colnames(mat))
  if(length(in_both) < 2){
    stop('Not enough samples match between the tree and variant matrix. Tree tip labels and variant matrix column names must match.')
  }
  drop_from_tree = setdiff(tree$tip.label,in_both)
  if(length(drop_from_tree) > 0){
    tree = ape::drop.tip(tree,drop_from_tree)
    warning(paste('These samples were dropped from the tree:',paste(drop_from_tree,collapse=', ')))
  }
  drop_from_mat = setdiff(colnames(mat),in_both)
  if(length(drop_from_mat) > 0){
    mat = mat[,in_both]
    warning(paste('These samples were dropped from the matrix:',paste(drop_from_mat,collapse=', ')))
  }
  mat = mat[,tree$tip.label]
  return(list(tree=tree,mat=mat))
}
