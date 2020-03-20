## Package update
In this package version I have:

* Given change in as.factor=FALSE default behavior in R 4.0 code was changed to keep data correctly typed. 
* Added three new user-facing features. 
* Added more data.
* Updated documentation.

I tested the package on: 
* local OS X install, R 3.6.2 & devel
* win, R 3.5.3, 3.6.2 & devel

## Check results for OS X:
No errors, warnings, or notes.

## Check results for win:
0 ERROR | 0 WARNING | 1 NOTE

NOTE: 
```
Possibly mis-spelled words in DESCRIPTION:
  Cingolani (29:61)
  SnpEff (29:41)
  Pre (3:13)
  Saund (27:71)
  al (28:5)
  bGWAS (20:59, 23:23, 27:15, 27:36)
  et (27:77)
  multiallelic (25:5)
  pre (19:30, 24:26)
  prewas (20:67, 24:5, 26:34)

```

The NOTE is spurious because these words are all correctly spelled technical 
terms, names, or Latin.
