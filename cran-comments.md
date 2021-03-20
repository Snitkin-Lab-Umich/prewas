## Package update
Prewas was removed from CRAN on 2020-07-05 because a dependency, vcfR, was archived on CRAN. Since then vcfR has been updated on CRAN and so prewas will now be able to install successfully. 

I tested the package on: 
* local OS X install, R 3.6.2
* local OS X install, R. 4.0.2
* win, R 3.6.3

## Check results for OS X:
No errors, warnings, or notes.

## Check results for win:
0 ERROR | 0 WARNING | 1 NOTE

NOTE: 
```
Possibly mis-spelled words in DESCRIPTION:
  Cingolani (29:61)
  SnpEff (29:41)
```

The NOTE is spurious because these words are all correctly spelled technical terms and names.
