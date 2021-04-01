## Resubmission

* Placed package names and software names in single quotes in DESCRIPTION
* Updated CITATION with current journal article information
* Explained that the word 'SnpEff' in the descriptoin text is not an acronym, but the name of a piece of software

## Package update
Prewas was archived from CRAN on 2020-07-05 because a dependency, vcfR, was archived on CRAN. Since then vcfR has been updated on CRAN and so prewas will now be able to install successfully. 

I tested the package on: 
* mac, R 3.6.2
* mac, R 4.0.2
* win, R dev

## Check results for mac:
No errors, warnings, or notes.

## Check results for win:
0 ERROR | 0 WARNING | 1 NOTE

NOTE: 
```
Possibly mis-spelled words in DESCRIPTION:
    Cingolani (29:66)
    Pre (3:13)
    Saund (27:71)
    SnpEff (29:46)
    al (28:5, 30:5)
    bGWAS (20:59, 23:23, 27:15, 27:36)
    et (27:77, 29:76)
    multiallelic (25:5)
    pre (19:30, 24:26)
    prewas (20:67, 24:5, 26:34, 28:57)
```

The NOTE is spurious because these words are all correctly spelled technical terms (bGWAS, multiallelic), Latin (et al), software (SnpEff, prewas), or last names (Cingolani, Saund).
