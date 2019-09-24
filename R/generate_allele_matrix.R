convert_fasta_to_matrix <- function(fasta){
  # TODO add code to check if fasta is a path to a fasta file, or is a phyDat object, or a DNAbin object
  # TODO add code to convertn phyDat and DNAbins to matrices
  if (class(fasta) == "character") {
    dna <- ape::read.dna(file = fasta, as.character = TRUE, format = "fasta")
  }

  # Convert object matrix
  # ADD check that dna object is a matrix now.
  colnames(dna) <- 1:ncol(dna)

  return(dna)
}

replace_non_ATGC_with_N <- function(mat){
  mat <- apply(mat, 2, toupper)
  mat[is.na(mat)] <- "N"
  mat[!(mat %in% c("A", "T", "G", "C"))] <- "N"
  return(mat)
}

identify_variant_sites <- function(mat){
  cols_to_keep <- apply(mat, 2, function(row){
    sum(unique(row) %in% c("A", "T", "G", "C"))
  })
  cols_to_keep <- cols_to_keep > 1
  cols_to_keep <- unname(cols_to_keep)
  return(cols_to_keep)
}

remove_invariant_sites <- function(mat, cols_to_keep){
  mat <- mat[, cols_to_keep, drop = FALSE]
  return(mat)
}

generate_allele_matrix <- function(fasta_path){
  dna_mat <- convert_fasta_to_matrix("data/test_dna.fna")
  dna_mat <- replace_non_ATGC_with_N(dna_mat)
  variant_site_log <- identify_variant_sites(dna_mat)
  variant_only_dna_mat <- remove_invariant_sites(dna_mat, variant_site_log)
  # Put SNP positions in the rows, isolates in columns
  t_variant_only_dna_mat <- t(variant_only_dna_mat)
  return(t_variant_only_dna_mat)
}

# Running this code:
# foo <- generate_allele_matrix("data/test_dna.fna")

# temp <- convert_fasta_to_matrix("data/test_dna.fna")
# temp <- replace_non_ATGC_with_N(temp)
# temp[,5] <- "T"
# variant_site_log <- identify_variant_sites(temp)
# variant_only_temp <- remove_invariant_sites(temp, variant_site_log)
