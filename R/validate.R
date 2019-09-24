check_is_number <- function(num){
  if (class(num) %in% c("numeric")) {
    stop("Input must be a numeric")
  }
}

is_file <- function(obj){
  is_file <- FALSE
  if (class(obj) == "character" & file.exists(obj)) {
      is_file <- TRUE
    }
  }
  return(is_file)
}

is_this_class <- function(obj, current_class){
  is_this_class <- FALSE
  if (class(obj) == current_class) {
    is_this_class <- TRUE
  }
  return(is_this_class)
}

check_inputs <- function(dna, tree, outgroup = NULL, gff = NULL){


  # DNA
  if (is.null(dna)) {
    stop("Must include a fasta file")
  } else {
    if (is_file(dna)) {
      # If DNA stored in FASTA file, read in and convert to matrix
      dna <- read_dna(dna)
    } else if (is_this_class(dna, 'DNAbin')) {
      dna <- convert_dnabin_to_matrix(dna)
    } else {
      stop("DNA input should be a path to a fasta file or a DNAbin object")
    }
  }

  # Tree
  if (is.null(tree)) {
    stop("Must include a tree")
  } else {
      if (is_file(tree)) {
        # If tree stored as .tree file, read in
        tree <- ape::read.tree(tree)
      }
      if (!is_this_class(tree, 'phylo')) {
        stop("Tree should either be a .tree or a phylo object")
      }
  }

 # Outgroup
  if (!is.null(outgroup)) {
    # If file, read it in
    if (is_file(outgroup)) {
      outgroup <- readLines(con = outgroup)
      outgroup <- as.character(outgroup)
    } else {
    # If object, check is a character string
      if (!is_this_class(outgroup, "character")) {
        stop("Outgroup must be a character string")
      }
    }

    # Check if outgroup is in the tree
    if (!outgroup %in% tree$tip.label) {
      warning("Outgroup not found in tree so ignoring it.")
      outgroup <- NULL
    }
  }

  # GFF
  if (!is.null(gff)) {
    gff <- read_gff(gff)
    if (ncol(gff) != 9) {
      stop("GFF file must have 9 columns")
    }
    if (sum(c(".", "+", "-") %in% unique(gff[, 7, drop = TRUE])) == 0) {
      stop("GFF file must only have strand information in column 7")
    }
  }


  results <- list("dna" = dna,
                    "tree" = tree,
                    "outgroup" = outgroup,
                    "gff" = gff)
  return(results)
}

read_gff <- function(gff_path){
  gff <- read.table(gff_path,
                    sep = '\t',
                    header = FALSE,
                    stringsAsFactors = FALSE,
                    quote = "",
                    comment.char = '#'
  )
  return(gff)
}

read_dna <- function(fasta_path){
  dna <- ape::read.dna(file = fasta_path, as.character = TRUE, format = "fasta")
  colnames(dna) <- 1:ncol(dna)
  dna <- t(dna)
  return(dna)
}

convert_dnabin_to_matrix <- function(dnabin){
  dna <- ape::as.alignment(dnabin)
  dna_mat <- matrix(NA, nrow = nrow(dnabin), ncol = nchar(dna$seq[1]))
  for (i in 1:nrow(dna_mat)) {
    dna_mat[i, ] <- strsplit(dna$seq[[i]], "")[[1]]
  }
  row.names(dna_mat) <- row.names(dnabin)
  colnames(dna_mat) <- 1:ncol(dna_mat)
  dna_mat <- t(dna_mat)
  return(dna_mat)
}
