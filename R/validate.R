check_is_number <- function(num){
  if (class(num) %in% c("numeric")) {
    stop("Input must be a numeric")
  }
}

is_file <- function(obj){
  is_file <- FALSE
  if (class(obj) == "character") {
    if (file.exists(obj)) {
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

check_is_this_class <- function(obj, current_class){
  if (class(obj) != current_class) {
    stop(paste('Input must be a',current_class))
  }

}

check_is_tree = function(tree){
  if (!is_this_class(tree,'phylo')) {
    stop('Input requires either a path to a tree file or an ape phylo object')
  }
}

check_tree_is_rooted = function(tree){
  if (!ape::is.rooted(tree)) {
    stop('Tree must be rooted.')
  }
}

check_inputs <- function(dna, tree, outgroup = NULL, gff = NULL){

  # DNA
  if (is.null(dna)) {
    stop("Must include a VCF file")
  } else {
    if (is_file(dna)) {
      # If DNA stored in VCF file, read in and convert to allele matrix
      dna <- load_vcf_file(dna)
    } else {
      stop("DNA input should be a path to a VCF 4.1 file")
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

    #subset gff on CDS
    gff <- subset_gff(gff)

    #get gene name from column 9
    gff <- clean_up_cds_name_from_gff(gff)
  }


  results <- list("dna" = dna,
                  "tree" = tree,
                  "outgroup" = outgroup,
                  "gff" = gff)
  return(results)
}

read_gff <- function(gff_path){
  gff <- utils::read.table(gff_path,
                           sep = '\t',
                           header = FALSE,
                           stringsAsFactors = FALSE,
                           quote = "",
                           comment.char = '#'
  )
  return(gff)
}

subset_gff <- function(gff){
  gff <- gff[gff[,3] == 'CDS',] # subset gff on CDS
  return(gff)
}

clean_up_cds_name_from_gff <- function(gff){
  cds_name <- apply(gff, 1, function(row){
    gsub('^ID=', '', row[9]) %>% gsub(';.*$', '', .)
  })
  gff[, 9] = cds_name
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

load_vcf_file <- function(vcf_path) {
  vcf <- vcfR::read.vcfR(file = vcf_path)
  vcf_geno_mat <- vcf@gt[, 2:ncol(vcf@gt), drop = FALSE]
  row.names(vcf_geno_mat) <- vcf@fix[, colnames(vcf@fix) == "POS", drop = TRUE]
  vcf_ref_allele <- vcf@fix[, colnames(vcf@fix) == "REF", drop = TRUE]
  vcf_alt_allele <- vcf@fix[, colnames(vcf@fix) == "ALT", drop = TRUE]

  for (i in 1:nrow(vcf_geno_mat)) {
    alt_alleles <- strsplit(vcf_alt_allele[i], split = ",")
    vcf_alt_allele_1 <- alt_alleles[[1]][1]
    vcf_alt_allele_2 <- alt_alleles[[1]][2]
    vcf_alt_allele_3 <- alt_alleles[[1]][3]

    vcf_geno_mat[i, vcf_geno_mat[i, ] == "0"] <- vcf_ref_allele[i]
    vcf_geno_mat[i, vcf_geno_mat[i, ] == "1"] <- vcf_alt_allele_1
    vcf_geno_mat[i, vcf_geno_mat[i, ] == "2"] <- vcf_alt_allele_2
    vcf_geno_mat[i, vcf_geno_mat[i, ] == "3"] <- vcf_alt_allele_3
  }
  return(vcf_geno_mat)
}

convert_vcf_to_allele_matrix <- function(vcf){

}
