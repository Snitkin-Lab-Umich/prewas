#' VCF
#'
#' Example VCF data.
#' @format TODO add Description
#' \describe{
#'   \item{temp}{longer description}
#'  }
"vcf"

#' Phylogenetic tree
#'
#' Example phylogenetic tree.
#' @format TODO add Description
#' \describe{
#'   \item{temp}{longer description}
#'   \item{temp2}{longer description}
#'  }
"tree"

#' Name of outgroup in the phylogenetic tree.
#'
#' Example VCF data.
#' @format Character string
#' \describe{
#'   \item{outgroup}{Name of the outgroup in the phylogenetic tree.}
#'  }
"outgroup"

#' FASTA file
#'
#' Example FASTA file.
#' @format TODO add Description
#' \describe{
#'   \item{temp}{longer description}
#'  }
"fasta"

#' GFF3 file
#'
#' Example GFF3 data.
#' @format TODO add Description
#' \describe{
#'   \item{temp}{longer description}
#'  }
"gff"


#' Results from running prewas() on the example data.
#'
#' Output from prewas().
#'
#' @format List of 6 objects.
#' \describe{
#'   \item{allele_mat}{Allele matrix (nucleotides stored as characters). Multiallelic sites represented on multiple lines in the matrix. TODO: what are rows? What are columns?}
#'   \item{bin_mat}{Binary matrix (nucleotides stored as 0 or 1). Multiallelic sites represented on multiple lines in the matrix. TODO: what are rows? What are columns?}
#'   \item{ar_results}{TODO - add description - is this a list? Vector? Length?}
#'   \item{dup}{TODO - add description}
#'   \item{gene_mat}{Gene-based matrix. Genes with any SNP stored as 1, Genes without SNPs stored as 0). TODO: What are rows? What are columns?}
#'  }
"prewas_results"
