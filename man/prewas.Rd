% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prewas.R
\name{prewas}
\alias{prewas}
\title{Preprocess SNPs before bGWAS}
\usage{
prewas(
  dna,
  tree = NULL,
  outgroup = NULL,
  gff = NULL,
  anc = TRUE,
  snpeff_grouping = NULL,
  grp_nonref = FALSE
)
}
\arguments{
\item{dna}{`Character` or `vcfR`. Required input. Path to VCF4.1 file or
`vcfR` object.}

\item{tree}{`NULL`, `character`, or `phylo`. Optional input. Ignored if
`NULL`. If `character` it should be a path to a .tree file. Defaults to
`NULL`.}

\item{outgroup}{`NULL` or `character`. Optional input. If `character` it
should be either a string naming the outgroup in the tree or a path to a
file containing only the outgroup name. Ignored if `NULL`. Defaults to
`NULL`.}

\item{gff}{`NULL`, `character`, `matrix`, or `data.frame`. Optional input. If
`NULL` it is ignored. If `character` it should be a path to a GFF3 file. If
 a `matrix` or `data.frame` it should be the GFF information stored in 9
 columns with the genes as rows. Defaults to `NULL`.}

\item{anc}{`Logical`. Optional input. When `TRUE` prewas performs ancestral
reconstruction. When `FALSE` prewas calculates the major allele. Defaults
to `TRUE`.}

\item{snpeff_grouping}{`NULL`, `character`. Optional input. Only used when a
SnpEff annotated multivcf is inputted. Use when you want to group SNPs by
gene and SnpEff impact. If `NULL` no custom-grouped gene matrix will be
generated. Options for input are a vector combination of 'HIGH',
'MODERATE', LOW', 'MODIFER'. Must write the impact combinations in all caps
(e.g. c('HIGH', 'MODERATE')). Defaults to `NULL`.}

\item{grp_nonref}{`Logical`. Optional input. When `TRUE` prewas collapses all
non-reference alleles for multi-allelic sites. When `FALSE` prewas keeps
multi-allelic sites separate. Defaults to `FALSE`.}
}
\value{
A list with the following items:
  \describe{
    \item{allele_mat}{`matrix`. An allele matrix, created from the vcf where
    each multiallelic site will be on its own line. The rowname will be the
    position of the variant in the vcf file. If the position is triallelic,
    there will be two rows containing the same information. The rows will be
    labeled "pos" and "pos.1". If the position is quadallelic, there will be
    three rows containing the same information. The rows will be labeled
    "pos", "pos.1", and "pos.2"}
    \item{bin_mat}{`matrix`. A binary matrix, where 0 is the reference allele
    and 1 indicates a variant. The dimensions may not match the `allele_mat`
    if the gff file is provided, because SNPs in overlapping genes are
    represented on multiple lines in `bin_mat`; in that case both position
    and locus tag name are provided in the rowname.}
    \item{ar_results}{`data.frame`. This data.frame records the alleles used
    as the reference alleles. Rows correspond to variant loci. If `anc =
    TRUE` the data.frame has two columns which contain the ancestrally
    reconstructed allele and the probability of the reconstruction. If `anc =
    FALSE` there is only one column which contains the major allele.}
    \item{dup}{`integer`. A vector of integers. It's an index that identifies
    duplicated rows. If the index is unique (appears once), that means it is
    not a multiallelic site. If the index appears more than once, that means
    the row was replicated `x` times, where `x` is the number of alternative
    alleles. Note: the multiple indices indicates multiallelic site splits,
    not overlapping genes splits.}
    \item{gene_mat}{`NULL` or `matrix` or `list`. `NULL` if no gene
    information provided (`gff = NULL` and no SnpEff annotation provided in
    VCF). If gene information is provided by a gff then a gene matrix is
    generated where each row is a gene and each column is a sample. If gene
    information is provided by a SnpEff annotated vcf then a list of up to
    six gene matrices are returned. The first matrix, `gene_mat_all` is a
    gene matrix for all variants.`gene_mat_modifier` is a gene matrix for
    only variants annoted as MODIFIER by SnpEff. Similarly there is a
    `gene_mat_low`, `gene_mat_moderate`, and `gene_mat_high.` If the user
    asks for a combination SnpEff annotations the final `gene_mat_custom`
    will contain that matrix.}
    \item{tree}{`NULL` or `phylo`. If `anc = FALSE` no tree is use or
    generated and the function returns `NULL`. If `anc = TRUE` and the user
    provides a tree but no outgroup: the function returns the tree after
    midpoint rooting. If `anc = TRUE` and the user provides both a tree and
    an outgroup: the function returns a tree rooted on the outgroup and the
    outgroup is removed from the tree. If the user does not provide a tree
    and `anc = TRUE` the function returns the midpoint rooted tree
    generated.}
  }
}
\description{
prewas is a tool to standardize the pre-processing of your
  genomic data before performing a bacterial genome-wide association study
  (bGWAS). prewas creates a variant matrix (where each row is a variant, each
  column is a sample, and the entries are presence - 1 - or absence - 0 - of
  the variant) that can be used as input for bGWAS tools. When creating the
  binary variant matrix, prewas can perform 3 pre-processing steps including:
  dealing with  multiallelic SNPs, (optional) dealing with SNPs in
  overlapping genes, and choosing a reference allele. prewas can output
  matrices for use with both SNP-based bGWAS and gene-based bGWAS. prewas can
  also provide gene matrices for variants with specific SnpEff annotations
  (Cingolani et al. 2012).
}
\examples{
vcf = prewas::vcf
gff = prewas::gff
tree = prewas::tree
outgroup = prewas::outgroup
output <- prewas(dna = vcf,
                 tree = tree,
                 outgroup = outgroup,
                 gff = gff,
                 anc = FALSE ,
                 grp_nonref = FALSE)
}
