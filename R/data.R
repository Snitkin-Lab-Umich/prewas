#' Nucleotide variants in example genome samples
#'
#' An example dataset containing 326 variants from 14 genome samples.
#'
#' @format vcfR class object with three sections:
#' \describe{
#'   \item{meta}{The metadata for the VCF file including the file format version
#'   number}
#'   \item{fix}{A character matrix with 326 rows and 8 columns. Contains
#'   information on chromosome (CHROM), genome position (POS), reference genome
#'   allele (REF), and alternative allele (ALT)}
#'   \item{gt}{A character matrix with 326 rows and 15 columns. Presence/absence
#'   for each variant defined in fix. Colnames are sample IDs.}
#'  }
"vcf"

#' Phylogenetic tree of example genomes
#'
#' Example unrooted phylogenetic tree.
#'
#' @format An ape phylo object with 14 tips.
"tree"

#' Name of outgroup in the phylogenetic tree.
#'
#' @format Character string.
"outgroup"

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
#' @format List of 5 objects.
#' \describe{
#'   \item{allele_mat}{Matrix. Character matrix of nucleotides (alleles).
#'   Multiallelic sites represented on multiple lines in the matrix. Dim: 360 x
#'   13. Rows are genomic loci. Columns are samples. Row names include only
#'   genomic position and do not have gene information.}
#'   \item{bin_mat}{Matrix. Binary matrix (nucleotides stored as 0 or 1).
#'   Multiallelic sites represented on multiple lines in the matrix. Alleles in
#'   overlapping genes are represened on multiple lines in the matrix. Rownames
#'   include genomic position and gene. Dim: 1016 x 13. Rows are genomic loci.
#'   Columns are samples.}
#'   \item{ar_results}{Data.frame. Dim: 360 x 1. Rows are genomic loci. The
#'   column is the major allele at that position. If anc=TRUE, then this object
#'   would be a 306 x 2 data.frame where the first column is the ancestral
#'   allele at that position inferred from ancestral reconstruction and the
#'   second column is the maximum likelihood probability.}
#'   \item{dup}{Integer vector. Length = 360. The number refers to the original
#'   genomic loci in the VCF file. The occurence count of the number is one less
#'   than the number of alleles. Ex: the 1st genomic locus (Position "1") occurs
#'   once in `dup` indicating that this is a biallelic site. In contrast, the
#'   5th genomic locus in the vcf (Position 18) occurs twice indicating that
#'   this is a triallelic site (represented in two rows: 18 and 18.1)}
#'   \item{gene_mat}{Matrix. Gene-based matrix. Genes with any SNP stored as 1,
#'   genes without SNPs stored as 0. Rows are genes. Columns are samples. Dim:
#'   96 x 13.}
#'  }
"results"
