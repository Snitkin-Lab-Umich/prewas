# prewas: data pre-processing for more informative bacterial GWAS

## Introduction
The prewas R package allows users to create a binary SNP matrix from a whole genome alignment. The SNP matrix includes the following features: (1) multiple line representation of multiallelic sites, (2) multiple line represention for SNPs present in overlapping genes, and (3) choice over the reference allele. Additionally, users can collapse SNPs into genes so the output is a binary gene matrix. Output from the prewas package should be used as the input to bacterial GWAS tools such as [hogwash](https://github.com/katiesaund/hogwash). The package is currently under development by the [Snitkin Lab](http://thesnitkinlab.com/) - please check for updates in December 2019. 
  
## Installation  
To install prewas follow these commands in R:  
 
 ```
## install devtools 
install.packages("devtools", dep=TRUE)
library(devtools)

## install prewas from github:
install_github("Snitkin-Lab-Umich/prewas")
library(prewas)
```

## Documentation
Documentation (Github wiki + vignette) as well as a publication are forthcoming.

## Contributors
Katie Saund, Stephanie Theide, and Zena Lapp contributed to this code.
