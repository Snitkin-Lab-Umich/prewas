# [prewas: data pre-processing for more informative bacterial GWAS](https://www.biorxiv.org/content/10.1101/2019.12.20.873158v1)


[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/prewas)](https://cran.r-project.org/package=prewas)
[![](http://cranlogs.r-pkg.org/badges/grand-total/prewas)](https://cran.r-project.org/package=prewas)

## Introduction
The prewas R package allows users to create a binary SNP matrix from a whole genome alignment. The SNP matrix includes the following features: (1) multiple line representation of multiallelic sites, (2) multiple line representation for SNPs present in overlapping genes, and (3) choice over the reference allele. Additionally, users can collapse SNPs into genes so the output is a binary gene matrix. Output from the prewas package should be used as the input to bacterial GWAS tools such as [hogwash](https://github.com/katiesaund/hogwash).
  
## Installation  
To install prewas follow these commands in R:  
 
```
install.packages("prewas")
library(prewas)
```

Note: this package depends on R (>= 3.5.0).

## Documentation
prewas is described in the paper: ["prewas: data pre-processing for more informative bacterial GWAS"](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000368). The Rscripts and data for the paper's figures and analyses can be found [in the manuscript analysis repository](https://github.com/Snitkin-Lab-Umich/prewas_manuscript_analysis).

A tutorial explaining how to use the package can be found in the  [vignette](https://github.com/Snitkin-Lab-Umich/prewas/blob/master/vignettes/getting_started_with_prewas.Rmd). 

## Contributors
[Katie Saund](https://github.com/katiesaund), [Stephanie Thiede](https://github.com/sthiede), and [Zena Lapp](https://github.com/zenalapp) contributed to this code.

## Questions

Please open an [issue](https://github.com/Snitkin-Lab-Umich/prewas/issues) or contact Evan Snitkin (esnitkin@umich.edu) with any questions. 
