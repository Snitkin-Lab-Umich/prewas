#' Nucleotide variants in example genome samples
#'
#' An example dataset containing 313 variants from 13 genome samples.
#'
#' @format vcfR class object with three sections:
#' \describe{
#'   \item{meta}{The metadata for the VCF file including the file format version
#'   number}
#'   \item{fix}{A character matrix with 313 rows and 8 columns. Contains
#'   information on chromosome (CHROM), genome position (POS), reference genome
#'   allele (REF), and alternative allele (ALT)}
#'   \item{gt}{A character matrix with 313 rows and 14 columns. Presence/absence
#'   for each variant defined in fix. Colnames are sample IDs.}
#'  }
"vcf"

#' Phylogenetic tree of example genomes
#'
#' Example rooted phylogenetic tree.
#'
#' @format An ape phylo object with 13 tips.
"tree"

#' Name of outgroup in the phylogenetic tree.
#'
#' @format Character string.
"outgroup"

#' Multiple sequence alignment of example genomes.
#'
#' Example DNA fasta data for example genomes.
#'
#' @format A phyDat object with 13 genomes and 216 nucleotides each
"fasta"

#' GFF3 file for example genomes
#'
#' An example of GFF3 formatted genome information.
#'
#' @format A character matrix with 110 rows and 9 columns:
#' \describe{
#'   \item{Column 1}{Chromosome}
#'   \item{Column 2}{Data source}
#'   \item{Column 3}{Feature type}
#'   \item{Column 4}{Feature start position}
#'   \item{Column 5}{Feature stop position}
#'   \item{Column 6}{Score}
#'   \item{Column 7}{Strand}
#'   \item{Column 8}{Phase}
#'   \item{Column 9}{Locus ID}
#'  }
"gff"


#' Results from running prewas() on the example data.
#'
#' Output from prewas().
#'
#' @format List of 6 objects.
#' \describe{
#'   \item{allele_mat}{Allele matrix (nucleotides stored as characters).
#'   Multiallelic sites represented on multiple lines in the matrix. TODO: what
#'   are rows? What are columns?}
#'   \item{bin_mat}{Binary matrix (nucleotides stored as 0 or 1). Multiallelic
#'   sites represented on multiple lines in the matrix. TODO: what are rows?
#'   What are columns?}
#'   \item{ar_results}{TODO - add description - is this a list? Vector? Length?}
#'   \item{dup}{TODO - add description}
#'   \item{gene_mat}{Gene-based matrix. Genes with any SNP stored as 1, Genes
#'   without SNPs stored as 0). TODO: What are rows? What are columns?}
#'  }
"prewas_results"
