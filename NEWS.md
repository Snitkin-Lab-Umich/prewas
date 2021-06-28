# prewas 1.0.0

* Added a `NEWS.md` file to track changes to the package.
* Initial CRAN submission 

# prewas 1.1.0

* Given updates to as.factor default behavior in R 4.0 code was changed to keep data correctly typed. 
* New feature: prewas can now subset variants by SnpEff annotation.
* New feature: Added an optional flag, grp_nonref. Its value tells prewas whether keep all alleles at multiallelic sites separate or to collapse all non-reference alleles for multi-allelic sites. 
* New feature: prewas can now handle bcftools-prepared indel vcf file and prepare these variants in the same manner as SNPs.

# prewas 1.1.1

* Prewas was archived from CRAN because a dependency was archived from CRAN. The dependent package is now updated and therefore prewas can be on CRAN again. 
