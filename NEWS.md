# prewas 1.0.0

* Added a `NEWS.md` file to track changes to the package.
* Initial CRAN submission 

# prewas 1.0.0.9000

* Added snpeff functionality. 
* Added an optional flag, grp_nonref. Its value tells prewas whether keep all alleles at multiallelic sites separate or to collapse all non-reference alleles for multi-allelic sites. 
* Prewas can now take bcftools-prepared indel vcf file and prepare these variants in the same manner as SNPs.
* Given updates to as.factor default behavior in R 4.0 code was added to keep data correctly typed. 
